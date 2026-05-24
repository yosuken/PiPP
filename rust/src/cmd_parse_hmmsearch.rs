//! `pipp_util parse-hmmsearch` — Rust port of script/parse_hmmsearch.rb.
//!
//! Parses an hmmsearch text output and produces:
//!   - all-hit.tsv         (every domain hit that passes domain-level filters)
//!   - best-hit.tsv        (non-overlapping best hits, with linked-hit merging)
//!   - best-hit.whole.fa   (one full-length sequence per query that had any best-hit)
//!   - best-hit.fa         (sliced sequence per best-hit region)
//!   - evalues.tsv         (when --create-evalue-table is passed; gene × hmm matrix)
//!
//! Two passes over the hmmsearch file are used: pass 1 collects gene IDs that
//! beat the gene-level evalue threshold, so pass 2 (fasta) only keeps those
//! sequences in memory. Pass 3 re-reads hmmsearch for per-domain detail.

use anyhow::{anyhow, Context, Result};
use std::collections::HashMap;
use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(clap::Args)]
pub struct Args {
    /// hmmsearch output file (concatenated, optionally .gz is NOT supported)
    #[arg(short = 'i', long = "hmmsearchout")]
    pub hmmsearch_out: PathBuf,

    /// Protein FASTA
    #[arg(short = 'f', long = "fasta")]
    pub fasta: PathBuf,

    /// Output directory (created if missing)
    #[arg(short = 'o', long = "outdir")]
    pub outdir: PathBuf,

    /// Gene-level i-Evalue threshold (default: 1e-2)
    #[arg(long = "gene-evalue", short = 'g', default_value_t = 1e-2)]
    pub gene_evalue: f64,

    /// Domain-level i-Evalue threshold (default: 10.0)
    #[arg(short = 'e', long = "evalue", default_value_t = 10.0)]
    pub evalue: f64,

    /// Emit evalues.tsv (matrix of gene × hmm best i-Evalues)
    #[arg(long = "create-evalue-table")]
    pub create_evalue_table: bool,

    // Linked-hit filters
    #[arg(long = "min-hmm-len", default_value_t = 0)]
    pub min_hmm_len: usize,
    #[arg(long = "min-hmm-cov", default_value_t = 0.0)]
    pub min_hmm_cov: f64,
    #[arg(long = "min-ali-len", default_value_t = 0)]
    pub min_ali_len: usize,
    #[arg(long = "min-ali-cov", default_value_t = 0.0)]
    pub min_ali_cov: f64,

    // Domain-level filters
    #[arg(long = "min-hmm-len-dom", default_value_t = 0)]
    pub min_hmm_len_dom: usize,
    #[arg(long = "min-hmm-cov-dom", default_value_t = 0.0)]
    pub min_hmm_cov_dom: f64,
    #[arg(long = "min-ali-len-dom", default_value_t = 0)]
    pub min_ali_len_dom: usize,
    #[arg(long = "min-ali-cov-dom", default_value_t = 0.0)]
    pub min_ali_cov_dom: f64,

    // Best-hit overlap thresholds (alignment side)
    #[arg(long = "max-ali-ovp-frc", default_value_t = 0.2)]
    pub max_ali_ovp_frc: f64,
    #[arg(long = "max-ali-ovp-len", default_value_t = 100_000)]
    pub max_ali_ovp_len: usize,

    // Link detection thresholds (hmm side)
    #[arg(long = "max-hmm-ovp-frc", default_value_t = 0.2)]
    pub max_hmm_ovp_frc: f64,
    #[arg(long = "max-hmm-ovp-len", default_value_t = 100_000)]
    pub max_hmm_ovp_len: usize,
}

// ---- types ----------------------------------------------------------------

#[derive(Default)]
struct HmmMeta {
    acc: String,
    desc: String,
    len: usize,
}

#[derive(Clone)]
struct DomainHit {
    hmm: String,
    // Verbatim hmmsearch text for output-only columns. Ruby keeps these
    // as strings (`values_at` on a split line) so we preserve them byte
    // for byte to match. Ruby's `Float#to_s` happens to roundtrip these
    // exactly, so emitting the original string also matches Ruby's
    // formatted output byte for byte.
    score: String,
    bias: String,
    c_evalue: String,
    i_evalue_str: String,
    env_fm: String,
    env_to: String,
    acc: String,
    full_evalue_str: String,
    full_score: String,
    // Parsed copies used for filtering, sorting, and range arithmetic.
    i_evalue: f64,
    hmm_fm: usize,
    hmm_to: usize,
    ali_fm: usize,
    ali_to: usize,
    full_evalue: f64,
}

// ---- formatting helper ----------------------------------------------------

/// Apply Ruby's `Float#to_s` cosmetic normalization to a numeric string,
/// without parsing it to f64:
///   * scientific notation (contains `e`): ensure mantissa has at least
///     one fractional digit ("4e-41" → "4.0e-41") and the exponent is
///     signed and zero-padded to >= 2 digits ("1e-6" → "1.0e-06").
///   * decimal notation (no `e`): ensure a `.` is present ("100" → "100.0").
///
/// hmmsearch's summary table emits some values without the `.0` Ruby
/// expects (e.g. `4e-41` for full-Evalue), so the raw text doesn't match
/// Ruby's `"4e-41".to_f.to_s == "4.0e-41"`. This pure-string normalizer
/// closes that gap for the two summary-table columns Ruby reformats
/// (full-Evalue and full-score), without doing a lossy parse/format
/// roundtrip.
fn ruby_format_normalize(s: &str) -> String {
    if s.is_empty() {
        return s.to_string();
    }
    match s.find('e') {
        None => {
            if s.contains('.') {
                s.to_string()
            } else {
                format!("{s}.0")
            }
        }
        Some(e_pos) => {
            let mantissa_raw = &s[..e_pos];
            let exp_str = &s[e_pos + 1..];
            let mantissa = if mantissa_raw.contains('.') {
                mantissa_raw.to_string()
            } else {
                format!("{mantissa_raw}.0")
            };
            let (sign, digits) = match exp_str.as_bytes().first().copied() {
                Some(b'-') => ('-', &exp_str[1..]),
                Some(b'+') => ('+', &exp_str[1..]),
                _ => ('+', exp_str),
            };
            let padded = if digits.len() < 2 {
                format!("{digits:0>2}")
            } else {
                digits.to_string()
            };
            format!("{mantissa}e{sign}{padded}")
        }
    }
}

// ---- range helpers --------------------------------------------------------

/// Returns overlap length of two inclusive `[a..=b]` integer ranges, or 0 if
/// they don't overlap.
fn overlap_len(a: (usize, usize), b: (usize, usize)) -> usize {
    let lo = a.0.max(b.0);
    let hi = a.1.min(b.1);
    if hi >= lo {
        hi - lo + 1
    } else {
        0
    }
}

/// Merge pairs (i, j) interpreted as inclusive ranges into the smallest set of
/// non-overlapping ranges. Mirrors `merge_ranges` in the Ruby version.
fn merge_pair_ranges(mut pairs: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if pairs.len() < 2 {
        return pairs;
    }
    pairs.sort_by_key(|p| p.0);
    let mut out: Vec<(usize, usize)> = Vec::with_capacity(pairs.len());
    for p in pairs {
        if let Some(last) = out.last_mut() {
            // overlap or contiguous touching is considered "overlap" for our
            // "pair indices" use (Ruby treats touching as no overlap, but pairs
            // generated by consecutive (i..i+1) overlap on the shared index).
            if overlap_len(*last, p) > 0 {
                last.1 = last.1.max(p.1);
                continue;
            }
        }
        out.push(p);
    }
    out
}

// ---- pass 1: scan summary -------------------------------------------------

/// Pass 1: read the per-query "Scores for complete sequences" tables and
/// collect every gene ID that beats `gene_evalue`.
fn pass1_collect_gids(path: &Path, gene_evalue: f64) -> Result<HashMap<String, ()>> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let rdr = BufReader::new(f);
    let mut gids: HashMap<String, ()> = HashMap::new();
    let mut flag = "";
    for line in rdr.lines() {
        let line = line?;
        let trimmed = line.trim_start();
        if trimmed.starts_with("Scores for complete sequences") {
            flag = "full";
        } else if trimmed.contains("--- inclusion threshold ---") {
            flag = "above_inc";
        } else if line.trim().is_empty() {
            flag = "";
        } else if flag == "full" && starts_with_indent_digit(&line) {
            // columns: full-evalue, score, bias, best-evalue, score, bias, exp, N, sequence, [desc...]
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() < 9 {
                continue;
            }
            let full_e: f64 = match cols[0].parse() {
                Ok(v) => v,
                Err(_) => continue,
            };
            if full_e >= gene_evalue {
                continue;
            }
            let gid = cols[8].to_string();
            gids.insert(gid, ());
        }
    }
    Ok(gids)
}

/// A line in the summary table starts with whitespace then a digit/decimal,
/// e.g. "    4.6e-12   51.2 ...". This guard mirrors the Ruby regex `^\s+\d`.
fn starts_with_indent_digit(line: &str) -> bool {
    let mut chars = line.chars();
    let mut saw_ws = false;
    for c in chars.by_ref() {
        if c == ' ' || c == '\t' {
            saw_ws = true;
            continue;
        }
        return saw_ws && c.is_ascii_digit();
    }
    false
}

// ---- pass 2: fasta --------------------------------------------------------

#[derive(Default)]
struct FastaStore {
    /// Preserves the file order so output mirrors the input fasta.
    order: Vec<String>,
    info: HashMap<String, String>,
    seq: HashMap<String, String>,
}

impl FastaStore {
    fn len_of(&self, gid: &str) -> usize {
        self.seq.get(gid).map(String::len).unwrap_or(0)
    }
}

fn pass2_read_fasta(path: &Path, want: &HashMap<String, ()>) -> Result<FastaStore> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let rdr = BufReader::new(f);
    let mut store = FastaStore::default();
    let mut cur: Option<String> = None;
    let mut buf = String::new();
    for line in rdr.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            // flush previous
            if let Some(gid) = cur.take() {
                if want.contains_key(&gid) {
                    store.seq.insert(gid, std::mem::take(&mut buf));
                } else {
                    buf.clear();
                }
            }
            let mut split = rest.splitn(2, char::is_whitespace);
            let gid = split.next().unwrap_or("").to_string();
            let info = split.next().unwrap_or("").trim().to_string();
            if want.contains_key(&gid) {
                store.order.push(gid.clone());
                store.info.insert(gid.clone(), info);
            }
            cur = Some(gid);
        } else {
            // append sequence line (strip whitespace)
            for ch in line.chars() {
                if !ch.is_whitespace() {
                    buf.push(ch);
                }
            }
        }
    }
    if let Some(gid) = cur.take() {
        if want.contains_key(&gid) {
            store.seq.insert(gid, buf);
        }
    }
    Ok(store)
}

// ---- pass 3: full hmmsearch parse ----------------------------------------

struct ParsedHmmsearch {
    /// HMM names in encounter order.
    hmm_order: Vec<String>,
    hmm_meta: HashMap<String, HmmMeta>,
    /// hits[gid][hmm] = vec of DomainHit (already filtered by per-domain rules).
    hits: HashMap<String, HashMap<String, Vec<DomainHit>>>,
}

fn pass3_parse_hmmsearch(path: &Path, args: &Args, fa: &FastaStore) -> Result<ParsedHmmsearch> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let rdr = BufReader::new(f);

    let mut out = ParsedHmmsearch {
        hmm_order: Vec::new(),
        hmm_meta: HashMap::new(),
        hits: HashMap::new(),
    };

    // inc_gids[gid][hmm] = (full_evalue_f64, full_evalue_str, full_score_str)
    // for genes above the gene-level inclusion threshold under the current
    // Query. The f64 is needed to filter against `--gene-evalue`; the
    // strings are kept verbatim from hmmsearch and used in output, which
    // matches Ruby exactly (Ruby parses to Float then Float#to_s, which
    // roundtrips for these well-formed values).
    let mut inc_gids: HashMap<String, HashMap<String, (f64, String, String)>> = HashMap::new();

    let mut flag = "";
    let mut cur_hmm = String::new();
    let mut cur_gid = String::new();
    let mut skipped_rows: usize = 0;

    for line in rdr.lines() {
        let line = line?;

        if let Some(rest) = line.strip_prefix("Query:") {
            // "Query:       <name>  [M=<len>]"
            let rest = rest.trim_start();
            let mut tokens = rest.split_whitespace();
            let name = tokens.next().unwrap_or("").to_string();
            // find [M=<int>]
            let mut hlen: usize = 0;
            for tok in tokens {
                if let Some(s) = tok.strip_prefix("[M=") {
                    if let Some(s) = s.strip_suffix(']') {
                        hlen = s.parse().unwrap_or(0);
                    }
                }
            }
            cur_hmm = name.clone();
            out.hmm_order.push(name.clone());
            out.hmm_meta.entry(name).or_insert_with(|| HmmMeta {
                len: hlen,
                ..HmmMeta::default()
            });
            // reset per-query state
            inc_gids.clear();
            flag = "";
        } else if let Some(rest) = line.strip_prefix("Accession:") {
            let v = rest.trim().to_string();
            if let Some(meta) = out.hmm_meta.get_mut(&cur_hmm) {
                meta.acc = v;
            }
        } else if let Some(rest) = line.strip_prefix("Description:") {
            let v = rest.trim().to_string();
            if let Some(meta) = out.hmm_meta.get_mut(&cur_hmm) {
                meta.desc = v;
            }
        } else if line
            .trim_start()
            .starts_with("Scores for complete sequences")
        {
            flag = "parse_full";
        } else if line.contains("--- inclusion threshold ---") {
            flag = "above_inclusion_threshold";
        } else if line.trim().is_empty() {
            flag = "";
        } else if let Some(rest) = line.strip_prefix(">>") {
            // "geneX" — start parsing per-domain table for this gene under cur_hmm
            cur_gid = rest.split_whitespace().next().unwrap_or("").to_string();
            let pass = inc_gids
                .get(&cur_gid)
                .and_then(|h| h.get(&cur_hmm))
                .map(|(e, _, _)| *e < args.gene_evalue)
                .unwrap_or(false);
            flag = if pass { "parse_each" } else { "" };
        } else if flag == "parse_full" && starts_with_indent_digit(&line) {
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() < 9 {
                continue;
            }
            let Ok(full_e) = cols[0].parse::<f64>() else {
                skipped_rows += 1;
                debug_assert!(false, "summary row full-evalue parse failed: {line:?}");
                continue;
            };
            if full_e >= args.gene_evalue {
                continue;
            }
            // Normalize summary-table strings to Ruby's Float#to_s form
            // (e.g. "4e-41" → "4.0e-41"). Ruby reformats these columns;
            // pass-through verbatim alone is not enough.
            let full_e_str = ruby_format_normalize(cols[0]);
            let score_str = ruby_format_normalize(cols[1]);
            let gid = cols[8].to_string();
            inc_gids
                .entry(gid)
                .or_default()
                .insert(cur_hmm.clone(), (full_e, full_e_str, score_str));
        } else if flag == "parse_each" && starts_with_indent_digit(&line) {
            // domain table row:
            // " #  ?  score  bias  c-Evalue  i-Evalue hmmfrom  hmmto    alifrom  alito    envfrom  envto    acc"
            //   0  1     2     3        4         5       6      7         8      9         10     11      12
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() < 16 {
                continue;
            }
            // Ruby code uses values_at(2..7, 9..10, 12..13, 15) → score, bias, c-Eval, i-Eval, hmm.fm, hmm.to, ali.fm, ali.to, env.fm, env.to, acc
            let Some(hit) = parse_domain_row(&cols, &cur_hmm) else {
                skipped_rows += 1;
                debug_assert!(false, "domain row parse failed: {line:?}");
                continue;
            };

            // domain-level filters
            if hit.i_evalue >= args.evalue {
                continue;
            }
            let hmm_len = out.hmm_meta.get(&cur_hmm).map(|m| m.len).unwrap_or(0);
            let g_len = fa.len_of(&cur_gid);
            if g_len == 0 {
                continue;
            }
            let hmmlen = hit.hmm_to.saturating_sub(hit.hmm_fm) + 1;
            let alilen = hit.ali_to.saturating_sub(hit.ali_fm) + 1;
            if hmmlen < args.min_hmm_len_dom {
                continue;
            }
            if alilen < args.min_ali_len_dom {
                continue;
            }
            if hmm_len == 0 || (hmmlen as f64) / (hmm_len as f64) < args.min_hmm_cov_dom {
                continue;
            }
            if (alilen as f64) / (g_len as f64) < args.min_ali_cov_dom {
                continue;
            }

            // Attach full-evalue/full-score (parsed + verbatim) for this
            // (gid, cur_hmm). At this point the gid passed the gene-level
            // filter so an inc_gids entry must exist; if it somehow
            // doesn't we fall back to zeroes (matches Ruby's `nil*"\t"`).
            let (fe_f64, fe_str, fs_str) = inc_gids
                .get(&cur_gid)
                .and_then(|h| h.get(&cur_hmm))
                .cloned()
                .unwrap_or((0.0, "0.0".to_string(), "0.0".to_string()));
            let mut hit = hit;
            hit.full_evalue = fe_f64;
            hit.full_evalue_str = fe_str;
            hit.full_score = fs_str;

            out.hits
                .entry(cur_gid.clone())
                .or_default()
                .entry(cur_hmm.clone())
                .or_default()
                .push(hit);
        }
    }

    // de-duplicate hmm_order (Query may repeat across concatenated outputs)
    let mut seen: HashMap<String, ()> = HashMap::new();
    out.hmm_order
        .retain(|h| seen.insert(h.clone(), ()).is_none());

    if skipped_rows > 0 {
        eprintln!(
            "  warning: skipped {skipped_rows} unparseable row(s) in {}",
            path.display()
        );
    }

    Ok(out)
}

/// Parse a domain table row. Numeric columns we actually use (i-Evalue
/// and the four hmm/ali coordinates) must parse; everything else is kept
/// verbatim. Returns None if any required parse fails; the caller bumps a
/// skip counter and continues.
fn parse_domain_row(cols: &[&str], cur_hmm: &str) -> Option<DomainHit> {
    Some(DomainHit {
        hmm: cur_hmm.to_string(),
        score: cols.get(2)?.to_string(),
        bias: cols.get(3)?.to_string(),
        c_evalue: cols.get(4)?.to_string(),
        i_evalue_str: cols.get(5)?.to_string(),
        i_evalue: cols.get(5)?.parse().ok()?,
        hmm_fm: cols.get(6)?.parse().ok()?,
        hmm_to: cols.get(7)?.parse().ok()?,
        ali_fm: cols.get(9)?.parse().ok()?,
        ali_to: cols.get(10)?.parse().ok()?,
        env_fm: cols.get(12)?.to_string(),
        env_to: cols.get(13)?.to_string(),
        acc: cols.get(15)?.to_string(),
        full_evalue: 0.0,
        full_evalue_str: String::new(),
        full_score: String::new(),
    })
}

// ---- output -------------------------------------------------------------

const ALL_HEADER: &[&str] = &[
    "protein",
    "length(aa)",
    "protein_info",
    "hmm_name",
    "hmm_acc",
    "hmm_desc",
    "hmm_len",
    "score",
    "bias",
    "c-Evalue",
    "i-Evalue",
    "hmm.fm",
    "hmm.to",
    "ali.fm",
    "ali.to",
    "env.fm",
    "env.to",
    "acc",
    "full-Evalue",
    "full-score",
];

fn write_outputs(args: &Args, fa: &FastaStore, ph: &ParsedHmmsearch) -> Result<()> {
    create_dir_all(&args.outdir)
        .with_context(|| format!("creating outdir {}", args.outdir.display()))?;

    let mut fwf = BufWriter::new(File::create(args.outdir.join("all-hit.tsv"))?);
    let mut fwb = BufWriter::new(File::create(args.outdir.join("best-hit.tsv"))?);
    let mut fwbw = BufWriter::new(File::create(args.outdir.join("best-hit.whole.fa"))?);
    let mut fwbf = BufWriter::new(File::create(args.outdir.join("best-hit.fa"))?);

    writeln!(fwf, "{}", ALL_HEADER.join("\t"))?;
    let mut best_header: Vec<&str> = ALL_HEADER.to_vec();
    best_header.push("link");
    best_header.push("region_name");
    writeln!(fwb, "{}", best_header.join("\t"))?;

    for gid in &fa.order {
        let Some(per_hmm) = ph.hits.get(gid) else {
            continue;
        };
        let glen = fa.len_of(gid);
        let ginfo = fa.info.get(gid).cloned().unwrap_or_default();

        // Flatten all per-domain hits.
        let mut infos: Vec<DomainHit> = per_hmm.values().flat_map(|v| v.iter().cloned()).collect();

        // ----- best-hit pass: greedy by i-Evalue ascending, filter ali overlap.
        infos.sort_by(|a, b| {
            a.i_evalue
                .partial_cmp(&b.i_evalue)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut binfos: Vec<DomainHit> = Vec::new();
        let mut alis: Vec<(usize, usize)> = Vec::new();
        for info in &infos {
            let r = (info.ali_fm, info.ali_to);
            let r_len = r.1 - r.0 + 1;
            let mut clash = false;
            for prev in &alis {
                let ovp = overlap_len(*prev, r);
                if ovp == 0 {
                    continue;
                }
                let prev_len = prev.1 - prev.0 + 1;
                let min_len = prev_len.min(r_len) as f64;
                let frc = ovp as f64 / min_len;
                if ovp > args.max_ali_ovp_len || frc > args.max_ali_ovp_frc {
                    clash = true;
                    break;
                }
            }
            if !clash {
                alis.push(r);
                binfos.push(info.clone());
            }
        }

        // ----- all-hit output: sort by [ali_fm, i-Evalue]
        let mut all_sorted = infos.clone();
        all_sorted.sort_by(|a, b| match a.ali_fm.cmp(&b.ali_fm) {
            std::cmp::Ordering::Equal => a
                .i_evalue
                .partial_cmp(&b.i_evalue)
                .unwrap_or(std::cmp::Ordering::Equal),
            o => o,
        });
        for info in &all_sorted {
            write_all_row(&mut fwf, gid, glen, &ginfo, info, ph)?;
        }

        // ----- best output: sort by [ali_fm, i-Evalue]
        binfos.sort_by(|a, b| match a.ali_fm.cmp(&b.ali_fm) {
            std::cmp::Ordering::Equal => a
                .i_evalue
                .partial_cmp(&b.i_evalue)
                .unwrap_or(std::cmp::Ordering::Equal),
            o => o,
        });

        // Link detection: per hmm, find consecutive non-reversed, non-overlapping
        // (in hmm coords) binfos and merge into one linked range.
        let mut link_cands: HashMap<String, Vec<usize>> = HashMap::new();
        for (idx, info) in binfos.iter().enumerate() {
            link_cands.entry(info.hmm.clone()).or_default().push(idx);
        }

        // idx -> (link_no, merged_ali_fm, merged_ali_to, merged_hmm_fm, merged_hmm_to)
        let mut idx_to_link: HashMap<usize, (usize, usize, usize, usize, usize)> = HashMap::new();

        for cand_idxs in link_cands.values() {
            if cand_idxs.len() <= 1 {
                continue;
            }
            let mut pairs: Vec<(usize, usize)> = Vec::new();
            for w in cand_idxs.windows(2) {
                let i = w[0];
                let j = w[1];
                let bi = &binfos[i];
                let bj = &binfos[j];
                if bi.hmm_fm > bj.hmm_fm {
                    continue;
                } // reversed
                let ir = (bi.hmm_fm, bi.hmm_to);
                let jr = (bj.hmm_fm, bj.hmm_to);
                let ovp = overlap_len(ir, jr);
                let mut overlapping = false;
                if ovp > 0 {
                    let i_len = ir.1 - ir.0 + 1;
                    let j_len = jr.1 - jr.0 + 1;
                    let frc = ovp as f64 / (i_len.min(j_len) as f64);
                    if ovp > args.max_hmm_ovp_len || frc > args.max_hmm_ovp_frc {
                        overlapping = true;
                    }
                }
                if !overlapping {
                    // store as index range over cand_idxs
                    let pi = cand_idxs.iter().position(|&x| x == i).unwrap();
                    let pj = cand_idxs.iter().position(|&x| x == j).unwrap();
                    pairs.push((pi, pj));
                }
            }
            let links = merge_pair_ranges(pairs);
            for (link_no, (lo, hi)) in links.iter().enumerate() {
                let members: Vec<usize> = (*lo..=*hi).map(|p| cand_idxs[p]).collect();
                let hmm_fm = members.iter().map(|&j| binfos[j].hmm_fm).min().unwrap();
                let hmm_to = members.iter().map(|&j| binfos[j].hmm_to).max().unwrap();
                let ali_fm = members.iter().map(|&j| binfos[j].ali_fm).min().unwrap();
                let ali_to = members.iter().map(|&j| binfos[j].ali_to).max().unwrap();
                for j in members {
                    idx_to_link.insert(j, (link_no + 1, ali_fm, ali_to, hmm_fm, hmm_to));
                }
            }
        }

        // ----- write best-hit + fasta. De-dup labels in fasta side.
        let mut seen_full: HashMap<String, ()> = HashMap::new();
        let mut seen_slice: HashMap<String, ()> = HashMap::new();
        for (idx, info) in binfos.iter().enumerate() {
            let (link_lab, ali_fm, ali_to, hmm_fm, hmm_to) =
                if let Some(&(n, afm, ato, hfm, hto)) = idx_to_link.get(&idx) {
                    (format!("link{n}_{}", info.hmm), afm, ato, hfm, hto)
                } else {
                    (
                        "".to_string(),
                        info.ali_fm,
                        info.ali_to,
                        info.hmm_fm,
                        info.hmm_to,
                    )
                };
            let hmm_len = ph.hmm_meta.get(&info.hmm).map(|m| m.len).unwrap_or(0);
            let hmmlen = hmm_to - hmm_fm + 1;
            let alilen = ali_to - ali_fm + 1;
            if hmmlen < args.min_hmm_len {
                continue;
            }
            if alilen < args.min_ali_len {
                continue;
            }
            if hmm_len == 0 || (hmmlen as f64) / (hmm_len as f64) < args.min_hmm_cov {
                continue;
            }
            if (alilen as f64) / (glen as f64) < args.min_ali_cov {
                continue;
            }

            // best-hit row keeps the *original* per-domain hmm.fm/hmm.to/
            // ali.fm/ali.to of `info`. Ruby uses the merged ranges only for
            // (a) the linked-length filter above and (b) the region_name.
            let region_name = format!("{gid}_fm{ali_fm}_to{ali_to}");
            write_best_row(
                &mut fwb,
                gid,
                glen,
                &ginfo,
                info,
                ph,
                &link_lab,
                &region_name,
            )?;

            // Fasta outputs
            let seq = fa.seq.get(gid).map(String::as_str).unwrap_or("");
            let lab_w = if ginfo.is_empty() {
                gid.clone()
            } else {
                format!("{gid} {ginfo}")
            };
            if seen_full.insert(lab_w.clone(), ()).is_none() {
                writeln!(fwbw, ">{lab_w}\n{seq}")?;
            }
            let lab_f = if ginfo.is_empty() {
                region_name.clone()
            } else {
                format!("{region_name} {ginfo}")
            };
            if seen_slice.insert(lab_f.clone(), ()).is_none() {
                let lo = ali_fm.saturating_sub(1);
                let hi = ali_to.min(seq.len());
                if lo < hi {
                    writeln!(fwbf, ">{lab_f}\n{slice}", slice = &seq[lo..hi])?;
                }
            }
        }
    }

    if args.create_evalue_table {
        write_evalue_table(args, fa, ph)?;
    }

    Ok(())
}

fn write_all_row<W: Write>(
    w: &mut W,
    gid: &str,
    glen: usize,
    ginfo: &str,
    info: &DomainHit,
    ph: &ParsedHmmsearch,
) -> Result<()> {
    let m = ph.hmm_meta.get(&info.hmm);
    let acc = m.map(|m| m.acc.as_str()).unwrap_or("");
    let desc = m.map(|m| m.desc.as_str()).unwrap_or("");
    let hmm_len = m.map(|m| m.len).unwrap_or(0);
    writeln!(
        w,
        "{gid}\t{glen}\t{ginfo}\t{hmm}\t{acc}\t{desc}\t{hlen}\t{score}\t{bias}\t{c}\t{i}\t{hf}\t{ht}\t{af}\t{at}\t{ef}\t{et}\t{ac}\t{fe}\t{fs}",
        gid = gid, glen = glen, ginfo = ginfo,
        hmm = info.hmm, acc = acc, desc = desc, hlen = hmm_len,
        score = info.score, bias = info.bias, c = info.c_evalue, i = info.i_evalue_str,
        hf = info.hmm_fm, ht = info.hmm_to, af = info.ali_fm, at = info.ali_to,
        ef = info.env_fm, et = info.env_to, ac = info.acc,
        fe = info.full_evalue_str, fs = info.full_score,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_best_row<W: Write>(
    w: &mut W,
    gid: &str,
    glen: usize,
    ginfo: &str,
    info: &DomainHit,
    ph: &ParsedHmmsearch,
    link: &str,
    region_name: &str,
) -> Result<()> {
    let m = ph.hmm_meta.get(&info.hmm);
    let acc = m.map(|m| m.acc.as_str()).unwrap_or("");
    let desc = m.map(|m| m.desc.as_str()).unwrap_or("");
    let hmm_len = m.map(|m| m.len).unwrap_or(0);
    writeln!(
        w,
        "{gid}\t{glen}\t{ginfo}\t{hmm}\t{acc}\t{desc}\t{hlen}\t{score}\t{bias}\t{c}\t{i}\t{hf}\t{ht}\t{af}\t{at}\t{ef}\t{et}\t{ac}\t{fe}\t{fs}\t{link}\t{region}",
        gid = gid, glen = glen, ginfo = ginfo,
        hmm = info.hmm, acc = acc, desc = desc, hlen = hmm_len,
        score = info.score, bias = info.bias, c = info.c_evalue, i = info.i_evalue_str,
        hf = info.hmm_fm, ht = info.hmm_to, af = info.ali_fm, at = info.ali_to,
        ef = info.env_fm, et = info.env_to, ac = info.acc,
        fe = info.full_evalue_str, fs = info.full_score,
        link = link, region = region_name,
    )?;
    Ok(())
}

fn write_evalue_table(args: &Args, fa: &FastaStore, ph: &ParsedHmmsearch) -> Result<()> {
    let mut w = BufWriter::new(File::create(args.outdir.join("evalues.tsv"))?);
    let mut hdr: Vec<String> = vec!["seq".to_string()];
    for hmm in &ph.hmm_order {
        let acc = ph.hmm_meta.get(hmm).map(|m| m.acc.as_str()).unwrap_or("");
        hdr.push(if acc.is_empty() {
            hmm.clone()
        } else {
            acc.to_string()
        });
    }
    writeln!(w, "{}", hdr.join("\t"))?;

    for gid in &fa.order {
        let mut row: Vec<String> = vec![gid.clone()];
        let per_hmm = ph.hits.get(gid);
        for hmm in &ph.hmm_order {
            // Mirror the Ruby reference exactly: argsort the i-Evalues to
            // get r = [original indices sorted ascending], then bidx =
            // position in r where r[idx] == 0 (= the rank of the
            // first-encountered domain hit). The emitted value is
            // i_evalue_str[bidx], NOT min(e). When the first-encountered
            // hit is also the most significant (the common case) these
            // coincide, but in other cases we deliberately reproduce the
            // Ruby behavior so verification stays byte-identical.
            let val = per_hmm
                .and_then(|h| h.get(hmm))
                .map(|v| {
                    let mut idx: Vec<usize> = (0..v.len()).collect();
                    idx.sort_by(|&a, &b| {
                        v[a].i_evalue
                            .partial_cmp(&v[b].i_evalue)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    let bidx = idx.iter().position(|&i| i == 0).unwrap_or(0);
                    // Ruby reformats evalues here via Float#to_s, so we
                    // run the same string-level normalization.
                    ruby_format_normalize(&v[bidx].i_evalue_str)
                })
                .unwrap_or_default();
            row.push(val);
        }
        writeln!(w, "{}", row.join("\t"))?;
    }
    Ok(())
}

// ---- entry point ---------------------------------------------------------

pub fn run(args: Args) -> Result<()> {
    if !args.hmmsearch_out.is_file() {
        return Err(anyhow!(
            "hmmsearch output not found: {}",
            args.hmmsearch_out.display()
        ));
    }
    if !args.fasta.is_file() {
        return Err(anyhow!("fasta not found: {}", args.fasta.display()));
    }

    eprintln!("pass 1: scanning summary for candidate gene IDs");
    let want = pass1_collect_gids(&args.hmmsearch_out, args.gene_evalue)?;
    eprintln!(
        "  {} candidate gene IDs above gene-evalue threshold",
        want.len()
    );

    eprintln!("pass 2: reading fasta");
    let fa = pass2_read_fasta(&args.fasta, &want)?;
    eprintln!("  {} sequences kept", fa.order.len());

    eprintln!("pass 3: parsing domain hits");
    let ph = pass3_parse_hmmsearch(&args.hmmsearch_out, &args, &fa)?;
    eprintln!(
        "  {} hmms, {} genes with hits",
        ph.hmm_order.len(),
        ph.hits.len()
    );

    eprintln!("writing outputs to {}", args.outdir.display());
    write_outputs(&args, &fa, &ph)?;
    eprintln!("done.");
    Ok(())
}

// ---- tests ---------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlap_len_basic() {
        assert_eq!(overlap_len((1, 10), (5, 15)), 6); // overlap 5..10
        assert_eq!(overlap_len((1, 10), (11, 20)), 0); // touching = no overlap
        assert_eq!(overlap_len((1, 10), (3, 5)), 3); // contained
        assert_eq!(overlap_len((10, 20), (1, 5)), 0); // disjoint
    }

    #[test]
    fn merge_pair_ranges_chains_consecutive() {
        // (0,1) and (1,2) share index 1 → should merge to (0,2)
        let m = merge_pair_ranges(vec![(0, 1), (1, 2), (4, 5)]);
        assert_eq!(m, vec![(0, 2), (4, 5)]);
    }

    #[test]
    fn ruby_format_normalize_matches_ruby_to_s() {
        // Adds .0 to integer mantissa in scientific notation.
        assert_eq!(ruby_format_normalize("4e-41"), "4.0e-41");
        assert_eq!(ruby_format_normalize("3e-34"), "3.0e-34");
        // Pads exponent to two digits.
        assert_eq!(ruby_format_normalize("1e-6"), "1.0e-06");
        assert_eq!(ruby_format_normalize("2.6e-7"), "2.6e-07");
        // Already canonical — pass through.
        assert_eq!(ruby_format_normalize("4.4e-41"), "4.4e-41");
        assert_eq!(ruby_format_normalize("1.0e-06"), "1.0e-06");
        // Positive exponent gets explicit +.
        assert_eq!(ruby_format_normalize("1e16"), "1.0e+16");
        // Decimal with `.` passes through.
        assert_eq!(ruby_format_normalize("0.0083"), "0.0083");
        // Integer-valued decimal gets `.0`.
        assert_eq!(ruby_format_normalize("100"), "100.0");
        // Empty stays empty.
        assert_eq!(ruby_format_normalize(""), "");
    }

    #[test]
    fn starts_with_indent_digit_handles_summary_rows() {
        assert!(starts_with_indent_digit("    4.6e-12   51.2"));
        assert!(starts_with_indent_digit("\t  3 0.5"));
        assert!(!starts_with_indent_digit("Query: name"));
        assert!(!starts_with_indent_digit(""));
        assert!(!starts_with_indent_digit("    abc"));
    }

    #[test]
    fn parse_domain_row_happy_path() {
        // 16 cols: idx, '?', score, bias, c-Evalue, i-Evalue,
        //          hmmfrom, hmmto, '..', alifrom, alito, '..', envfrom, envto, '..', acc
        let cols = [
            "1", "?", "51.2", "0.4", "1e-12", "6.3e-12", "1", "40", "..", "16", "57", "..", "5",
            "74", "..", "0.79",
        ];
        let h = parse_domain_row(&cols, "MyHmm").expect("should parse");
        assert_eq!(h.hmm, "MyHmm");
        assert_eq!(h.score, "51.2");
        assert_eq!(h.bias, "0.4");
        assert_eq!(h.c_evalue, "1e-12");
        assert_eq!(h.i_evalue_str, "6.3e-12");
        assert_eq!(h.i_evalue, 6.3e-12);
        assert_eq!(h.hmm_fm, 1);
        assert_eq!(h.hmm_to, 40);
        assert_eq!(h.ali_fm, 16);
        assert_eq!(h.ali_to, 57);
        assert_eq!(h.env_fm, "5");
        assert_eq!(h.env_to, "74");
        assert_eq!(h.acc, "0.79");
    }

    #[test]
    fn parse_domain_row_returns_none_on_malformed() {
        let cols = [
            "1", "?", "51.2", "0.4", "1e-12", "abc", // i-Evalue not parseable
            "1", "40", "..", "16", "57", "..", "5", "74", "..", "0.79",
        ];
        assert!(parse_domain_row(&cols, "MyHmm").is_none());
        let short = ["1", "?", "51.2"]; // not enough cols
        assert!(parse_domain_row(&short, "MyHmm").is_none());
    }
}
