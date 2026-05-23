//! End-to-end test for `pipp_util parse-hmmsearch`.
//!
//! Builds a synthetic hmmsearch output + fasta, runs the importer through
//! the library API, and checks the produced TSV / FASTA outputs.

use pipp_util::cmd_parse_hmmsearch::{run, Args};
use std::fs;
use std::path::PathBuf;

const HMMSEARCH_OUT: &str = "\
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.4 (Aug 2023); http://hmmer.org/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       MyHmm  [M=40]
Accession:   PF10417.9
Description: synthetic HMM for tests
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence  Description
    ------- ------ -----    ------- ------ -----   ---- --  --------  -----------
    4.6e-12   51.2   0.4    6.3e-12   50.7   0.4    1.2  1  gene1
    1e-05     22.6   0.5    2e-05     20.6   0.0    2.4  1  gene2

Domain annotation for each sequence (and alignments):
>> gene1
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   50.7   0.4    1e-13   6.3e-12       1      40 ..       5      45 ..       1      48 .. 0.99

>> gene2
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   20.6   0.0    1e-06     2e-05       3      38 ..      10      45 ..       7      48 .. 0.91

[ok]
";

const FASTA: &str = "\
>gene1 some description
ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM
>gene2
MKLAPCDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHI
>gene3_above_evalue_threshold
SSSSSSSSSSSSSSSSSSSS
";

fn default_args(hmm: PathBuf, fa: PathBuf, out: PathBuf) -> Args {
    Args {
        hmmsearch_out: hmm,
        fasta: fa,
        outdir: out,
        gene_evalue: 1e-2,
        evalue: 10.0,
        create_evalue_table: false,
        min_hmm_len: 0,
        min_hmm_cov: 0.0,
        min_ali_len: 0,
        min_ali_cov: 0.0,
        min_hmm_len_dom: 0,
        min_hmm_cov_dom: 0.0,
        min_ali_len_dom: 0,
        min_ali_cov_dom: 0.0,
        max_ali_ovp_frc: 0.2,
        max_ali_ovp_len: 100_000,
        max_hmm_ovp_frc: 0.2,
        max_hmm_ovp_len: 100_000,
    }
}

#[test]
fn end_to_end_produces_expected_outputs() {
    let tmp = tempfile::tempdir().unwrap();
    let hmm_path = tmp.path().join("hmmsearch.out");
    let fa_path = tmp.path().join("query.fa");
    let out_dir = tmp.path().join("out");
    fs::write(&hmm_path, HMMSEARCH_OUT).unwrap();
    fs::write(&fa_path, FASTA).unwrap();

    let args = default_args(hmm_path, fa_path, out_dir.clone());
    run(args).expect("parse-hmmsearch should succeed");

    // Files exist
    for f in [
        "all-hit.tsv",
        "best-hit.tsv",
        "best-hit.whole.fa",
        "best-hit.fa",
    ] {
        let p = out_dir.join(f);
        assert!(p.is_file(), "missing output: {}", p.display());
    }

    // all-hit.tsv: 1 header + 2 data rows (gene1, gene2 each have 1 domain hit)
    let all = fs::read_to_string(out_dir.join("all-hit.tsv")).unwrap();
    let mut lines = all.lines();
    let header = lines.next().expect("header");
    assert!(
        header.starts_with("protein\tlength(aa)"),
        "header: {header}"
    );
    let data: Vec<&str> = lines.collect();
    assert_eq!(data.len(), 2, "expected 2 data rows, got {}", data.len());

    // best-hit.tsv: 2 rows, with link + region_name columns
    let best = fs::read_to_string(out_dir.join("best-hit.tsv")).unwrap();
    let best_lines: Vec<&str> = best.lines().collect();
    assert_eq!(best_lines.len(), 3, "expected header + 2 best-hit rows");
    let header = best_lines[0];
    assert!(header.ends_with("\tlink\tregion_name"));

    // gene1 row: ali_fm=5, ali_to=45 → region "gene1_fm5_to45"
    let gene1_row = best_lines
        .iter()
        .find(|l| l.starts_with("gene1\t"))
        .expect("gene1 in best-hit");
    assert!(
        gene1_row.contains("\tgene1_fm5_to45"),
        "missing region_name: {gene1_row}"
    );

    // best-hit.whole.fa contains the full sequences
    let whole = fs::read_to_string(out_dir.join("best-hit.whole.fa")).unwrap();
    assert!(whole.contains(">gene1 some description"));
    assert!(whole.contains("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM"));
    assert!(whole.contains(">gene2"));

    // best-hit.fa is sliced: gene1 region is residues 5..=45 of the full seq
    let slice = fs::read_to_string(out_dir.join("best-hit.fa")).unwrap();
    let gene1_full = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM";
    let expected_slice = &gene1_full[4..45]; // 1-based [5..45] -> 0-based [4..45)
    assert!(
        slice.contains(expected_slice),
        "expected slice {expected_slice:?} not in best-hit.fa:\n{slice}"
    );

    // gene3 was below gene-evalue threshold → must not appear anywhere.
    assert!(!all.contains("gene3"));
    assert!(!best.contains("gene3"));
}

#[test]
fn evalue_table_is_optional() {
    let tmp = tempfile::tempdir().unwrap();
    let hmm_path = tmp.path().join("hmmsearch.out");
    let fa_path = tmp.path().join("query.fa");
    let out_dir = tmp.path().join("out");
    fs::write(&hmm_path, HMMSEARCH_OUT).unwrap();
    fs::write(&fa_path, FASTA).unwrap();

    // First run without the flag — evalues.tsv must not appear.
    let mut args = default_args(hmm_path.clone(), fa_path.clone(), out_dir.clone());
    run(args).unwrap();
    assert!(!out_dir.join("evalues.tsv").is_file());

    // Re-run with the flag (overwrite dir contents).
    args = default_args(hmm_path, fa_path, out_dir.clone());
    args.create_evalue_table = true;
    run(args).unwrap();
    let table = fs::read_to_string(out_dir.join("evalues.tsv")).unwrap();
    assert!(table.starts_with("seq\t"));
    assert!(table.contains("\ngene1\t"));
    assert!(table.contains("\ngene2\t"));
}
