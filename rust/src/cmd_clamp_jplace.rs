//! `pipp_util clamp-jplace` — make an APPLES-2 jplace consumable by gappa.
//!
//! Two problems can appear in jplace files coming out of the APPLES-2 path:
//!
//!  1. Non-finite numbers (`-nan`, `inf`, ...) in the `placements` arrays.
//!     These are produced by `gappa prepare unchunkify` (and sometimes by
//!     APPLES-2 itself on trees with negative branch lengths) and are not
//!     valid JSON, so neither serde_json nor genesis/gappa can even parse
//!     the file.
//!  2. Negative branch lengths in the `tree` string. genesis refuses these
//!     ("Invalid float number ...").
//!
//! Fix order matters: first replace non-finite tokens *in raw text* (only
//! outside JSON strings, so taxon names containing "inf"/"nan" are safe and
//! negative likelihoods stay intact), which makes the file valid JSON; then
//! parse it and clamp negative branch lengths in the `tree` field to 0.

use anyhow::{Context, Result};
use std::fs;
use std::path::PathBuf;

#[derive(clap::Args)]
pub struct Args {
    /// jplace file to rewrite in place
    pub jplace: PathBuf,
}

pub fn run(args: Args) -> Result<()> {
    let text = fs::read_to_string(&args.jplace)
        .with_context(|| format!("reading {}", args.jplace.display()))?;

    // Phase 1: non-finite -> 0 (raw text, outside strings) so the file parses.
    let (fixed, n_nonfinite) = replace_nonfinite_outside_strings(&text);

    let mut v: serde_json::Value = serde_json::from_str(&fixed)
        .with_context(|| format!("parsing {} (after non-finite fix)", args.jplace.display()))?;

    // Phase 2: clamp negative branch lengths in the tree field.
    let mut n_branch = 0;
    if let Some(tree) = v.get("tree").and_then(|t| t.as_str()) {
        let (clamped, n) = clamp_negative_branch_lengths(tree);
        n_branch = n;
        if n > 0 {
            v["tree"] = serde_json::Value::String(clamped);
        }
    }

    if n_nonfinite == 0 && n_branch == 0 {
        return Ok(()); // already clean; leave the file untouched
    }

    let out = serde_json::to_string(&v).context("serializing sanitized jplace")?;
    fs::write(&args.jplace, out).with_context(|| format!("writing {}", args.jplace.display()))?;
    eprintln!(
        "clamp-jplace: {n_nonfinite} non-finite value(s) -> 0, {n_branch} negative branch length(s) clamped in {}",
        args.jplace.display()
    );
    Ok(())
}

/// Replace `nan`/`inf`/`infinity` tokens (optionally signed, case-insensitive)
/// with `0`, but only when they occur *outside* a JSON string. Returns
/// (fixed_text, count). Tracks string state (with backslash escapes) so taxon
/// names inside the `tree`/`n` strings are never touched, and negative
/// numbers like `-1234.5` are left alone (the sign is only consumed when a
/// non-finite word follows it).
fn replace_nonfinite_outside_strings(text: &str) -> (String, usize) {
    let b = text.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(b.len());
    let mut i = 0;
    let mut in_str = false;
    let mut n = 0usize;
    while i < b.len() {
        let c = b[i];
        if in_str {
            out.push(c);
            if c == b'\\' && i + 1 < b.len() {
                out.push(b[i + 1]); // copy escaped char verbatim
                i += 2;
                continue;
            }
            if c == b'"' {
                in_str = false;
            }
            i += 1;
            continue;
        }
        if c == b'"' {
            in_str = true;
            out.push(c);
            i += 1;
            continue;
        }
        // optional sign, then a non-finite keyword
        let mut j = i;
        if b[j] == b'-' || b[j] == b'+' {
            j += 1;
        }
        if let Some(len) = match_nonfinite_word(&b[j..]) {
            out.push(b'0');
            n += 1;
            i = j + len;
            continue;
        }
        out.push(c);
        i += 1;
    }
    (
        String::from_utf8(out).expect("ascii-only edits keep utf-8 valid"),
        n,
    )
}

/// If `b` starts with "infinity"/"inf"/"nan" (case-insensitive), return the
/// matched length; otherwise None. Checks "infinity" before "inf".
fn match_nonfinite_word(b: &[u8]) -> Option<usize> {
    let n = b.len().min(8);
    let lower: Vec<u8> = b[..n].iter().map(|c| c.to_ascii_lowercase()).collect();
    if lower.starts_with(b"infinity") {
        Some(8)
    } else if lower.starts_with(b"inf") {
        Some(3)
    } else if lower.starts_with(b"nan") {
        Some(3)
    } else {
        None
    }
}

/// Replace every `:-<number>` branch length in a newick/jplace tree string
/// with `:0`. Returns (clamped_tree, count). Edge-number annotations `{N}`
/// are non-negative integers, so they are never matched; positive scientific
/// notation like `:1.5e-05` is not matched because the `-` does not
/// immediately follow `:`.
fn clamp_negative_branch_lengths(tree: &str) -> (String, usize) {
    let b = tree.as_bytes();
    let mut out: Vec<u8> = Vec::with_capacity(b.len());
    let mut i = 0;
    let mut n = 0usize;
    while i < b.len() {
        if b[i] == b':' && i + 1 < b.len() && b[i + 1] == b'-' {
            out.extend_from_slice(b":0");
            n += 1;
            i += 2;
            while i < b.len()
                && (b[i].is_ascii_digit() || matches!(b[i], b'.' | b'e' | b'E' | b'+' | b'-'))
            {
                i += 1;
            }
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    (String::from_utf8(out).expect("tree stays valid utf-8"), n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clamps_only_negative_branch_lengths() {
        let tree = "((A:-0.05772{1},B:0.17333{2}):0.11938{3},C:-1.5e-05{4}):0{5};";
        let (out, n) = clamp_negative_branch_lengths(tree);
        assert_eq!(out, "((A:0{1},B:0.17333{2}):0.11938{3},C:0{4}):0{5};");
        assert_eq!(n, 2);
    }

    #[test]
    fn positive_scientific_not_matched() {
        let tree = "(A:1.5e-05{0},B:2.0e-3{1});";
        let (out, n) = clamp_negative_branch_lengths(tree);
        assert_eq!(out, tree);
        assert_eq!(n, 0);
    }

    #[test]
    fn nonfinite_replaced_outside_strings_only() {
        // -nan in a placement, "inf" inside a name string must be preserved,
        // and a negative likelihood (-1234.5) must survive.
        let input = r#"{"tree":"(infA:-0.5{0},nanB:0.2{1}):0{2};","placements":[{"p":[[0,-1234.5,1,-nan,0]],"n":["infseq"]}]}"#;
        let (out, n) = replace_nonfinite_outside_strings(input);
        assert_eq!(n, 1); // only the -nan
        assert!(out.contains("infA:-0.5")); // name preserved
        assert!(out.contains("nanB:0.2")); // name preserved
        assert!(out.contains("\"infseq\"")); // name preserved
        assert!(out.contains("-1234.5")); // negative number preserved
        assert!(out.contains("[0,-1234.5,1,0,0]")); // -nan -> 0
        assert!(!out.contains("nan,"));
    }

    #[test]
    fn run_fixes_nan_and_negative_branches() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.jplace");
        let jplace = r#"{"version":3,"tree":"(A:-0.5{0},B:0.2{1}):0{2};","placements":[{"p":[[0,-1234.5,1,-nan,0]],"n":["q"]}],"fields":["edge_num","likelihood","like_weight_ratio","distal_length","pendant_length"]}"#;
        fs::write(&p, jplace).unwrap();
        run(Args { jplace: p.clone() }).unwrap();
        let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(&p).unwrap()).unwrap();
        assert_eq!(v["tree"], "(A:0{0},B:0.2{1}):0{2};");
        assert_eq!(v["placements"][0]["p"][0][1], -1234.5); // negative likelihood kept
        assert_eq!(v["placements"][0]["p"][0][3], 0); // -nan -> 0
    }

    #[test]
    fn clean_file_untouched() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("clean.jplace");
        let jplace = r#"{"tree":"(A:0.1{0},B:0.2{1}):0{2};","placements":[{"p":[[0,-100.0,1,0.5,0.1]],"n":["q"]}]}"#;
        fs::write(&p, jplace).unwrap();
        let before = fs::read_to_string(&p).unwrap();
        run(Args { jplace: p.clone() }).unwrap();
        let after = fs::read_to_string(&p).unwrap();
        assert_eq!(before, after); // byte-identical, not rewritten
    }
}
