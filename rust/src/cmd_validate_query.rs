//! `pipp_util validate-query` — Rust port of script/01-1a-A.validate_query.rb.
//!
//! Reads a (optionally gzipped) query FASTA in a single pass and:
//!   * writes a cleaned FASTA (`>{id} ` header — the trailing space is
//!     required so gappa chunkify doesn't parse a trailing `_<digit>` as an
//!     abundance — followed by the whitespace-stripped sequence)
//!   * writes the list of all sequence IDs (one per line, in input order)
//!   * filters out sequences shorter than `--min-seq-len`
//!   * errors on duplicate IDs
//!   * writes a JSON metadata file (name, counts, paths, MD5s)
//!
//! MD5s of both the input file (raw bytes, matching Ruby's
//! `Digest::MD5.file`) and the output FASTA are computed *while streaming*,
//! so the whole job is one pass over the input plus one pass writing the
//! output — versus the Ruby version's three passes (parse + 2× MD5 reread).

use anyhow::{anyhow, Context, Result};
use flate2::read::MultiGzDecoder;
use md5::{Digest, Md5};
use serde::Serialize;
use std::cell::RefCell;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::rc::Rc;

#[derive(clap::Args)]
pub struct Args {
    /// Input query FASTA (.gz / .gzip auto-detected)
    #[arg(short = 'i', long = "query")]
    pub query: PathBuf,

    /// Display name recorded in the JSON metadata
    #[arg(long)]
    pub name: String,

    /// Cleaned FASTA output path
    #[arg(long = "fa")]
    pub fa: PathBuf,

    /// JSON metadata output path
    #[arg(long = "json")]
    pub json: PathBuf,

    /// Sequence-ID list output path
    #[arg(long = "list")]
    pub list: PathBuf,

    /// Minimum amino-acid length; shorter sequences are dropped from the
    /// cleaned FASTA (but still listed in the ID list)
    #[arg(long = "min-seq-len", default_value_t = 0)]
    pub min_seq_len: usize,

    /// Maximum amino-acid length; longer sequences are dropped from the
    /// cleaned FASTA (but still listed in the ID list). Default 100000 is
    /// HMMER's hard protein-pipeline limit: a single sequence > 100000 aa
    /// makes hmmsearch abort with a fatal exception (p7_pipeline.c), taking
    /// the whole run down, so such sequences are skipped up front.
    #[arg(long = "max-seq-len", default_value_t = 100_000)]
    pub max_seq_len: usize,
}

#[derive(Serialize)]
struct Meta {
    name: String,
    numseq: usize,
    numtooshortseq: usize,
    numtoolongseq: usize,
    fasta: String,
    fjson: String,
    fasta_ori: String,
    #[serde(rename = "fasta_MD5")]
    fasta_md5: String,
    #[serde(rename = "fasta_ori_MD5")]
    fasta_ori_md5: String,
}

/// Reader wrapper that MD5-hashes every byte that passes through, before
/// any decompression — so the digest matches `Digest::MD5.file` on the raw
/// on-disk file (compressed bytes for a .gz input). The hasher is shared
/// via Rc/RefCell so it can be finalized after the decode chain (which
/// takes ownership of this reader) is dropped.
struct HashingReader<R> {
    inner: R,
    hasher: Rc<RefCell<Md5>>,
}
impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        if n > 0 {
            self.hasher.borrow_mut().update(&buf[..n]);
        }
        Ok(n)
    }
}

/// Writer wrapper that MD5-hashes everything written to it, giving the
/// output FASTA's digest without a re-read.
struct HashingWriter<W> {
    inner: W,
    hasher: Md5,
}
impl<W: Write> Write for HashingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let n = self.inner.write(buf)?;
        if n > 0 {
            self.hasher.update(&buf[..n]);
        }
        Ok(n)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

fn is_gzip(p: &Path) -> bool {
    let s = p.to_string_lossy();
    s.ends_with(".gz") || s.ends_with(".gzip")
}

fn absolute_path(p: &Path) -> PathBuf {
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|cwd| cwd.join(p))
            .unwrap_or_else(|_| p.to_path_buf())
    }
}

pub fn run(args: Args) -> Result<()> {
    // ---- input: raw-byte hasher -> (gunzip?) -> line reader ----
    let in_hasher = Rc::new(RefCell::new(Md5::new()));
    let raw = File::open(&args.query)
        .with_context(|| format!("opening query {}", args.query.display()))?;
    let hashing = HashingReader {
        inner: raw,
        hasher: Rc::clone(&in_hasher),
    };
    // The decode chain takes ownership of `hashing`; we recover the digest
    // after dropping the reader (only `in_hasher` clone remains).
    let mut reader: Box<dyn BufRead> = if is_gzip(&args.query) {
        Box::new(BufReader::new(MultiGzDecoder::new(hashing)))
    } else {
        Box::new(BufReader::new(hashing))
    };

    // ---- output: cleaned FASTA via hashing writer ----
    let fa_file =
        File::create(&args.fa).with_context(|| format!("creating fa {}", args.fa.display()))?;
    let mut fa_w = HashingWriter {
        inner: BufWriter::new(fa_file),
        hasher: Md5::new(),
    };

    let list_file = File::create(&args.list)
        .with_context(|| format!("creating list {}", args.list.display()))?;
    let mut list_w = BufWriter::new(list_file);

    let mut numseq = 0usize;
    let mut numtooshort = 0usize;
    let mut numtoolong = 0usize;
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();

    let mut cur_id: Option<String> = None;
    let mut seq = String::new();

    let mut line = Vec::<u8>::new();
    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .context("reading query line")?;
        if n == 0 {
            break;
        }
        // Strip the trailing newline for header detection; sequence lines
        // get all whitespace stripped anyway.
        let is_header = line.first() == Some(&b'>');
        if is_header {
            // flush previous record
            if let Some(id) = cur_id.take() {
                flush_record(
                    &mut fa_w,
                    &id,
                    &seq,
                    args.min_seq_len,
                    args.max_seq_len,
                    &mut numseq,
                    &mut numtooshort,
                    &mut numtoolong,
                )?;
            }
            seq.clear();
            // parse id: first whitespace-delimited token after '>'
            let header = String::from_utf8_lossy(&line[1..]);
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            if !seen.insert(id.clone()) {
                return Err(anyhow!(
                    "sequence id {id} is found twice. Ensure all query sequence ids are unique."
                ));
            }
            writeln!(list_w, "{id}").context("writing id list")?;
            cur_id = Some(id);
        } else {
            // append sequence with all whitespace removed
            for &b in &line {
                if !b.is_ascii_whitespace() {
                    seq.push(b as char);
                }
            }
        }
    }
    if let Some(id) = cur_id.take() {
        flush_record(
            &mut fa_w,
            &id,
            &seq,
            args.min_seq_len,
            args.max_seq_len,
            &mut numseq,
            &mut numtooshort,
            &mut numtoolong,
        )?;
    }

    fa_w.flush().context("flushing fa")?;
    list_w.flush().context("flushing list")?;

    // Finalize output MD5.
    let fasta_md5 = hex(fa_w.hasher.finalize().as_slice());

    // Finalize input MD5 from the streaming hash: drop the reader (and thus
    // the HashingReader's Rc clone), then take the hasher back out. No extra
    // pass over the input.
    drop(reader);
    let in_md5 = Rc::try_unwrap(in_hasher)
        .map_err(|_| anyhow!("input hasher still referenced"))?
        .into_inner();
    let fasta_ori_md5 = hex(in_md5.finalize().as_slice());

    let meta = Meta {
        name: args.name.clone(),
        numseq,
        numtooshortseq: numtooshort,
        numtoolongseq: numtoolong,
        fasta: args.fa.to_string_lossy().into_owned(),
        fjson: args.json.to_string_lossy().into_owned(),
        fasta_ori: absolute_path(&args.query).to_string_lossy().into_owned(),
        fasta_md5,
        fasta_ori_md5,
    };
    let json = serde_json::to_string(&meta).context("serializing metadata")?;
    let mut jf = File::create(&args.json)
        .with_context(|| format!("creating json {}", args.json.display()))?;
    writeln!(jf, "{json}").context("writing json")?;

    eprintln!(
        "validate-query: {numseq} kept, {numtooshort} too short (< {} aa), {numtoolong} too long (> {} aa)",
        args.min_seq_len, args.max_seq_len
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn flush_record<W: Write>(
    fa_w: &mut W,
    id: &str,
    seq: &str,
    min_len: usize,
    max_len: usize,
    numseq: &mut usize,
    numtooshort: &mut usize,
    numtoolong: &mut usize,
) -> Result<()> {
    let len = seq.len();
    if len < min_len {
        *numtooshort += 1;
    } else if len > max_len {
        // HMMER's protein pipeline aborts on target sequences > 100000 aa
        // (the default max_len), so these are dropped before hmmsearch.
        *numtoolong += 1;
    } else {
        *numseq += 1;
        // Trailing space after id is intentional (gappa chunkify abundance
        // parsing). Matches Ruby `fw.puts [">"+id+" ", seq]`.
        writeln!(fa_w, ">{id} ").context("writing fa header")?;
        writeln!(fa_w, "{seq}").context("writing fa seq")?;
    }
    Ok(())
}

fn hex(bytes: &[u8]) -> String {
    let mut s = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_tmp(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".fa").tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn cleans_and_counts() {
        let q = write_tmp(">s1 desc\nMKLA\nPCDE\n>s2\nAC\n>s3\nMMMMMM\n");
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("out.fa");
        let json = dir.path().join("out.json");
        let list = dir.path().join("out.list");
        run(Args {
            query: q.path().to_path_buf(),
            name: "q".into(),
            fa: fa.clone(),
            json: json.clone(),
            list: list.clone(),
            min_seq_len: 4,
            max_seq_len: 100_000,
        })
        .unwrap();

        // s1 = MKLAPCDE (8), s2 = AC (2, dropped), s3 = MMMMMM (6)
        let fa_s = std::fs::read_to_string(&fa).unwrap();
        assert_eq!(fa_s, ">s1 \nMKLAPCDE\n>s3 \nMMMMMM\n");
        // list has all three ids
        let list_s = std::fs::read_to_string(&list).unwrap();
        assert_eq!(list_s, "s1\ns2\ns3\n");
        // json counts
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
        assert_eq!(v["numseq"], 2);
        assert_eq!(v["numtooshortseq"], 1);
        assert_eq!(v["name"], "q");
        assert!(v["fasta_MD5"].as_str().unwrap().len() == 32);
        assert!(v["fasta_ori_MD5"].as_str().unwrap().len() == 32);
    }

    #[test]
    fn rejects_duplicate_ids() {
        let q = write_tmp(">dup\nAAAA\n>dup\nCCCC\n");
        let dir = tempfile::tempdir().unwrap();
        let r = run(Args {
            query: q.path().to_path_buf(),
            name: "q".into(),
            fa: dir.path().join("o.fa"),
            json: dir.path().join("o.json"),
            list: dir.path().join("o.list"),
            min_seq_len: 0,
            max_seq_len: 100_000,
        });
        assert!(r.is_err());
        assert!(r.unwrap_err().to_string().contains("twice"));
    }

    #[test]
    fn drops_too_long_sequences() {
        // s1 ok (10), s2 too long (15 > max 12), s3 ok (5)
        let long = "A".repeat(15);
        let q = write_tmp(&format!(">s1\nMKLAPCDEFG\n>s2\n{long}\n>s3\nMMMMM\n"));
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("o.fa");
        let json = dir.path().join("o.json");
        run(Args {
            query: q.path().to_path_buf(),
            name: "q".into(),
            fa: fa.clone(),
            json: json.clone(),
            list: dir.path().join("o.list"),
            min_seq_len: 0,
            max_seq_len: 12,
        })
        .unwrap();

        let fa_s = std::fs::read_to_string(&fa).unwrap();
        assert_eq!(fa_s, ">s1 \nMKLAPCDEFG\n>s3 \nMMMMM\n"); // s2 dropped
        let v: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json).unwrap()).unwrap();
        assert_eq!(v["numseq"], 2);
        assert_eq!(v["numtoolongseq"], 1);
        assert_eq!(v["numtooshortseq"], 0);
    }
}
