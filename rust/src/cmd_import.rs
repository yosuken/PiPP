use anyhow::{anyhow, Context, Result};
use duckdb::{params, Connection};
use std::path::{Path, PathBuf};

#[derive(clap::Args)]
pub struct Args {
    /// Path to result/<refpkg> directory (containing assign/, alignment/, feature/aa/)
    pub result_dir: PathBuf,

    /// Output DuckDB file path. Defaults to <result_dir>/pipp.duckdb
    #[arg(long)]
    pub db: Option<PathBuf>,

    /// Name of the refpkg recorded in each row. Defaults to basename of <result_dir>
    #[arg(long)]
    pub refpkg: Option<String>,

    /// Path to the refpkg's backbone.json (validate_refpkg metadata). When
    /// given and present, one row is loaded into the refpkgs table.
    #[arg(long = "refpkg-json")]
    pub refpkg_json: Option<PathBuf>,

    /// Path to run_params.json ({params:{...}, software:[...]}). When given and
    /// present, loads run-time options into run_params and tool path/version
    /// into the software table.
    #[arg(long = "run-json")]
    pub run_json: Option<PathBuf>,

    /// Overwrite an existing DB file
    #[arg(long)]
    pub overwrite: bool,
}

/// run_params.json: every run-time option as key/value plus per-tool metadata.
#[derive(serde::Deserialize)]
struct RunJson {
    #[serde(default)]
    params: std::collections::BTreeMap<String, String>,
    #[serde(default)]
    software: Vec<SoftwareEntry>,
}

#[derive(serde::Deserialize)]
struct SoftwareEntry {
    name: String,
    #[serde(default)]
    path: Option<String>,
    #[serde(default)]
    version: Option<String>,
}

/// Subset of validate_refpkg's backbone.json that is worth persisting:
/// identity + provenance. The run-local derived paths (faln/ftreME/ppdir...)
/// are intentionally omitted — they are cache locations, not provenance.
#[derive(serde::Deserialize)]
struct BackboneMeta {
    refpkg: Option<String>,
    hmmname: Option<String>,
    hmmlen: Option<String>, // stored as a string in backbone.json
    #[serde(rename = "falnO")]
    faln_o: Option<String>,
    #[serde(rename = "ftreO")]
    ftre_o: Option<String>,
    #[serde(rename = "fhmmO")]
    fhmm_o: Option<String>,
}

pub fn run(args: Args) -> Result<()> {
    let result_dir = args
        .result_dir
        .canonicalize()
        .with_context(|| format!("result_dir not accessible: {}", args.result_dir.display()))?;

    let refpkg = args.refpkg.unwrap_or_else(|| {
        result_dir
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "unknown".to_string())
    });

    let db_path = args.db.unwrap_or_else(|| result_dir.join("pipp.duckdb"));
    if db_path.exists() {
        if args.overwrite {
            std::fs::remove_file(&db_path)
                .with_context(|| format!("removing existing DB: {}", db_path.display()))?;
        } else {
            return Err(anyhow!(
                "DB file already exists: {} (use --overwrite to replace)",
                db_path.display()
            ));
        }
    }

    let conn =
        Connection::open(&db_path).with_context(|| format!("opening DB: {}", db_path.display()))?;
    crate::schema::create_all(&conn)?;

    let assign_tsv = result_dir.join("assign/per_query.tsv");
    if assign_tsv.is_file() {
        let n = import_assignments(&conn, &refpkg, &assign_tsv)?;
        eprintln!("imported {} rows into assignments", n);
    } else {
        eprintln!("skip assignments: {} not found", assign_tsv.display());
    }

    let aa_tsv = result_dir.join("feature/aa/feature.tsv");
    if aa_tsv.is_file() {
        let n = import_aa_features(&conn, &refpkg, &aa_tsv)?;
        eprintln!("imported {} rows into aa_features", n);
    } else {
        eprintln!("skip aa_features: {} not found", aa_tsv.display());
    }

    let pos_tsv = result_dir.join("alignment/aligned_position.tsv");
    if pos_tsv.is_file() {
        let n = import_aligned_positions(&conn, &refpkg, &pos_tsv)?;
        eprintln!("imported {} rows into aligned_positions", n);
    } else {
        eprintln!("skip aligned_positions: {} not found", pos_tsv.display());
    }

    // jplace clamp records (one *.clamp.tsv per sanitized jplace; absent when
    // nothing was clamped).
    let n = import_jplace_clamps(&conn, &refpkg, &result_dir.join("placement"))?;
    if n > 0 {
        eprintln!("imported {} rows into jplace_clamps", n);
    }

    // refpkg identity/provenance from backbone.json (when supplied).
    if let Some(json) = &args.refpkg_json {
        if json.is_file() {
            import_refpkg_meta(&conn, &refpkg, json)?;
            eprintln!("imported 1 row into refpkgs");
        } else {
            eprintln!("skip refpkgs: {} not found", json.display());
        }
    }

    // run-time options + software path/version from run_params.json.
    if let Some(json) = &args.run_json {
        if json.is_file() {
            let (np, ns) = import_run_json(&conn, &refpkg, json)?;
            eprintln!("imported {np} rows into run_params, {ns} into software");
        } else {
            eprintln!("skip run_params/software: {} not found", json.display());
        }
    }

    eprintln!("wrote {}", db_path.display());
    Ok(())
}

fn import_jplace_clamps(conn: &Connection, refpkg: &str, placement_dir: &Path) -> Result<usize> {
    if !placement_dir.is_dir() {
        return Ok(0);
    }
    let mut app = conn.appender("jplace_clamps")?;
    let mut n = 0usize;
    let mut entries: Vec<PathBuf> = std::fs::read_dir(placement_dir)
        .with_context(|| format!("reading {}", placement_dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.to_string_lossy().ends_with(".clamp.tsv"))
        .collect();
    entries.sort();
    for path in entries {
        let mut rdr = tsv_reader(&path)?;
        for rec in rdr.records() {
            let rec = rec.with_context(|| format!("reading {}", path.display()))?;
            if rec.len() < 3 {
                continue;
            }
            let jplace = &rec[0];
            let nf: i64 = rec[1].trim().parse().unwrap_or(0);
            let nb: i64 = rec[2].trim().parse().unwrap_or(0);
            app.append_row(params![refpkg, jplace, nf, nb])?;
            n += 1;
        }
    }
    app.flush()?;
    Ok(n)
}

fn import_refpkg_meta(conn: &Connection, refpkg: &str, json_path: &Path) -> Result<()> {
    let txt = std::fs::read_to_string(json_path)
        .with_context(|| format!("reading {}", json_path.display()))?;
    let m: BackboneMeta = serde_json::from_str(&txt)
        .with_context(|| format!("parsing backbone.json: {}", json_path.display()))?;
    let hmmlen: Option<i64> = m.hmmlen.as_deref().and_then(|s| s.trim().parse().ok());
    let mut app = conn.appender("refpkgs")?;
    app.append_row(params![
        refpkg, m.refpkg, m.hmmname, hmmlen, m.faln_o, m.ftre_o, m.fhmm_o
    ])?;
    app.flush()?;
    Ok(())
}

fn import_run_json(conn: &Connection, refpkg: &str, json_path: &Path) -> Result<(usize, usize)> {
    let txt = std::fs::read_to_string(json_path)
        .with_context(|| format!("reading {}", json_path.display()))?;
    let r: RunJson = serde_json::from_str(&txt)
        .with_context(|| format!("parsing run_params.json: {}", json_path.display()))?;

    let mut np = 0usize;
    {
        let mut app = conn.appender("run_params")?;
        for (k, v) in &r.params {
            app.append_row(params![refpkg, k, v])?;
            np += 1;
        }
        app.flush()?;
    }

    let mut ns = 0usize;
    {
        let mut app = conn.appender("software")?;
        for s in &r.software {
            app.append_row(params![refpkg, s.name, s.path, s.version])?;
            ns += 1;
        }
        app.flush()?;
    }
    Ok((np, ns))
}

fn tsv_reader(path: &Path) -> Result<csv::Reader<std::fs::File>> {
    let rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_path(path)
        .with_context(|| format!("opening TSV: {}", path.display()))?;
    Ok(rdr)
}

fn parse_opt_f64(s: &str) -> Result<Option<f64>> {
    let t = s.trim();
    if t.is_empty() || t == "NA" || t == "nan" || t == "NaN" {
        return Ok(None);
    }
    Ok(Some(
        t.parse::<f64>()
            .with_context(|| format!("not a float: {s:?}"))?,
    ))
}

fn parse_opt_i64(s: &str) -> Result<Option<i64>> {
    let t = s.trim();
    if t.is_empty() || t == "NA" {
        return Ok(None);
    }
    Ok(Some(
        t.parse::<i64>()
            .with_context(|| format!("not an int: {s:?}"))?,
    ))
}

fn import_assignments(conn: &Connection, refpkg: &str, path: &Path) -> Result<usize> {
    let mut rdr = tsv_reader(path)?;
    let mut app = conn.appender("assignments")?;
    let mut n = 0usize;
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec.with_context(|| format!("reading row {} of {}", i + 2, path.display()))?;
        if rec.len() < 6 {
            return Err(anyhow!(
                "{}: row {} has {} columns, expected 6",
                path.display(),
                i + 2,
                rec.len()
            ));
        }
        let name = &rec[0];
        let lwr: f64 = rec[1]
            .parse()
            .with_context(|| format!("LWR at row {}", i + 2))?;
        let fract: f64 = rec[2]
            .parse()
            .with_context(|| format!("fract at row {}", i + 2))?;
        let a_lwr: f64 = rec[3]
            .parse()
            .with_context(|| format!("aLWR at row {}", i + 2))?;
        let a_fract: f64 = rec[4]
            .parse()
            .with_context(|| format!("afract at row {}", i + 2))?;
        let taxopath = &rec[5];
        app.append_row(params![refpkg, name, lwr, fract, a_lwr, a_fract, taxopath])?;
        n += 1;
    }
    app.flush()?;
    Ok(n)
}

fn import_aa_features(conn: &Connection, refpkg: &str, path: &Path) -> Result<usize> {
    // Header (28 cols):
    // gene len len_of_std_aa avg_MW N-ARSC C-ARSC S-ARSC
    //   K R H D E N Q S T Y A V L I P F M W G C others
    let mut rdr = tsv_reader(path)?;
    let mut app = conn.appender("aa_features")?;
    let mut n = 0usize;
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec.with_context(|| format!("reading row {} of {}", i + 2, path.display()))?;
        if rec.len() < 28 {
            return Err(anyhow!(
                "{}: row {} has {} columns, expected 28",
                path.display(),
                i + 2,
                rec.len()
            ));
        }
        let gene = &rec[0];
        let len = parse_opt_i64(&rec[1])?;
        let len_std_aa = parse_opt_i64(&rec[2])?;
        let avg_mw = parse_opt_f64(&rec[3])?;
        let n_arsc = parse_opt_f64(&rec[4])?;
        let c_arsc = parse_opt_f64(&rec[5])?;
        let s_arsc = parse_opt_f64(&rec[6])?;
        // 20 AAs + others
        let mut aa: Vec<Option<i64>> = Vec::with_capacity(21);
        for j in 7..28 {
            aa.push(parse_opt_i64(&rec[j])?);
        }
        app.append_row(params![
            refpkg, gene, len, len_std_aa, avg_mw, n_arsc, c_arsc, s_arsc, aa[0], aa[1], aa[2],
            aa[3], aa[4], aa[5], aa[6], aa[7], aa[8], aa[9], aa[10], aa[11], aa[12], aa[13],
            aa[14], aa[15], aa[16], aa[17], aa[18], aa[19], aa[20],
        ])?;
        n += 1;
    }
    app.flush()?;
    Ok(n)
}

fn import_aligned_positions(conn: &Connection, refpkg: &str, path: &Path) -> Result<usize> {
    // Header: query <pos_label_1> ... <pos_label_k> fract taxpath
    // Each pos column cell may be a comma-separated set of residues.
    let mut rdr = tsv_reader(path)?;
    let header = rdr.headers()?.clone();
    if header.len() < 3 {
        return Err(anyhow!(
            "{}: header has {} columns, expected >= 3 (query, <pos labels>, fract, taxpath)",
            path.display(),
            header.len()
        ));
    }
    let n_pos = header.len() - 3;
    let pos_labels: Vec<String> = (1..=n_pos).map(|i| header[i].to_string()).collect();
    let fract_idx = n_pos + 1;
    let taxpath_idx = n_pos + 2;

    let mut app = conn.appender("aligned_positions")?;
    let mut n = 0usize;
    for (i, rec) in rdr.records().enumerate() {
        let rec = rec.with_context(|| format!("reading row {} of {}", i + 2, path.display()))?;
        if rec.len() != header.len() {
            return Err(anyhow!(
                "{}: row {} has {} columns, header has {}",
                path.display(),
                i + 2,
                rec.len(),
                header.len()
            ));
        }
        let query = &rec[0];
        let fract = parse_opt_f64(&rec[fract_idx])?;
        let taxpath_raw = rec[taxpath_idx].trim();
        let taxpath: Option<&str> = if taxpath_raw.is_empty() || taxpath_raw == "NA" {
            None
        } else {
            Some(taxpath_raw)
        };

        for (j, label) in pos_labels.iter().enumerate() {
            let residues = &rec[j + 1];
            let pos_index = j as i64; // 0-based column order in the source TSV
            app.append_row(params![
                refpkg, query, pos_index, label, residues, fract, taxpath
            ])?;
            n += 1;
        }
    }
    app.flush()?;
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn fresh_conn() -> Connection {
        let conn = Connection::open_in_memory().unwrap();
        crate::schema::create_all(&conn).unwrap();
        conn
    }

    fn write_tsv(contents: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".tsv").tempfile().unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        f.flush().unwrap();
        f
    }

    fn count(conn: &Connection, table: &str) -> i64 {
        conn.query_row(&format!("SELECT COUNT(*) FROM {table}"), [], |r| r.get(0))
            .unwrap()
    }

    // --- parse_opt_* helpers ---

    #[test]
    fn parse_opt_f64_handles_values_and_missing() {
        assert_eq!(parse_opt_f64("1.5").unwrap(), Some(1.5));
        assert_eq!(parse_opt_f64("0").unwrap(), Some(0.0));
        assert_eq!(parse_opt_f64("  2.0 ").unwrap(), Some(2.0));
        assert_eq!(parse_opt_f64("").unwrap(), None);
        assert_eq!(parse_opt_f64("NA").unwrap(), None);
        assert_eq!(parse_opt_f64("nan").unwrap(), None);
        assert_eq!(parse_opt_f64("NaN").unwrap(), None);
        assert!(parse_opt_f64("abc").is_err());
    }

    #[test]
    fn parse_opt_i64_handles_values_and_missing() {
        assert_eq!(parse_opt_i64("42").unwrap(), Some(42));
        assert_eq!(parse_opt_i64("-3").unwrap(), Some(-3));
        assert_eq!(parse_opt_i64("").unwrap(), None);
        assert_eq!(parse_opt_i64("NA").unwrap(), None);
        assert!(parse_opt_i64("1.5").is_err());
    }

    // --- import_assignments ---

    #[test]
    fn import_assignments_loads_rows_and_records_refpkg() {
        let conn = fresh_conn();
        let f = write_tsv(
            "name\tLWR\tfract\taLWR\tafract\ttaxopath\n\
             seq1\t0.0\t0.0\t1.0\t1.0\tBacteria\n\
             seq1\t0.5\t0.5\t0.5\t0.5\tBacteria;E.coli\n\
             seq2\t0.4\t0.4\t1.0\t1.0\tEukaryota\n",
        );
        let n = import_assignments(&conn, "rp_test", f.path()).unwrap();
        assert_eq!(n, 3);
        assert_eq!(count(&conn, "assignments"), 3);

        let refpkg: String = conn
            .query_row("SELECT DISTINCT refpkg FROM assignments", [], |r| r.get(0))
            .unwrap();
        assert_eq!(refpkg, "rp_test");

        let (lwr, taxopath): (f64, String) = conn
            .query_row(
                "SELECT LWR, taxopath FROM assignments WHERE query_name='seq1' AND aLWR=0.5",
                [],
                |r| Ok((r.get(0)?, r.get(1)?)),
            )
            .unwrap();
        assert_eq!(lwr, 0.5);
        assert_eq!(taxopath, "Bacteria;E.coli");
    }

    #[test]
    fn import_assignments_rejects_short_rows() {
        let conn = fresh_conn();
        let f = write_tsv("name\tLWR\tfract\taLWR\tafract\ttaxopath\nseq1\t0.5\t0.5\n");
        // csv reader is strict (flexible=false), so the short row triggers
        // an error during record read; we just verify the import fails.
        assert!(import_assignments(&conn, "rp", f.path()).is_err());
    }

    // --- import_aa_features ---

    #[test]
    fn import_aa_features_loads_28_cols() {
        let conn = fresh_conn();
        let header = "gene\tlen\tlen_of_std_aa\tavg_MW\tN-ARSC\tC-ARSC\tS-ARSC\t\
                      K\tR\tH\tD\tE\tN\tQ\tS\tT\tY\tA\tV\tL\tI\tP\tF\tM\tW\tG\tC\tothers";
        let row1 = "seq1\t100\t100\t120.5\t0.5\t0.3\t0.1\t\
                    5\t5\t2\t8\t8\t5\t3\t7\t6\t3\t7\t7\t8\t6\t5\t4\t2\t1\t5\t3\t0";
        let f = write_tsv(&format!("{header}\n{row1}\n"));

        let n = import_aa_features(&conn, "rp", f.path()).unwrap();
        assert_eq!(n, 1);

        let (gene, len, avg_mw, aa_k, aa_others): (String, i32, f64, i32, i32) = conn
            .query_row(
                "SELECT query_name, len, avg_mw, aa_K, aa_others FROM aa_features",
                [],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?, r.get(4)?)),
            )
            .unwrap();
        assert_eq!(gene, "seq1");
        assert_eq!(len, 100);
        assert_eq!(avg_mw, 120.5);
        assert_eq!(aa_k, 5);
        assert_eq!(aa_others, 0);
    }

    #[test]
    fn import_aa_features_accepts_missing_values() {
        let conn = fresh_conn();
        let header = "gene\tlen\tlen_of_std_aa\tavg_MW\tN-ARSC\tC-ARSC\tS-ARSC\t\
                      K\tR\tH\tD\tE\tN\tQ\tS\tT\tY\tA\tV\tL\tI\tP\tF\tM\tW\tG\tC\tothers";
        let row = "seq_empty\tNA\tNA\tNA\tNA\tNA\tNA\t\
                   NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
        let f = write_tsv(&format!("{header}\n{row}\n"));
        let n = import_aa_features(&conn, "rp", f.path()).unwrap();
        assert_eq!(n, 1);

        let len: Option<i32> = conn
            .query_row("SELECT len FROM aa_features", [], |r| r.get(0))
            .unwrap();
        assert_eq!(len, None);
    }

    // --- import_aligned_positions ---

    #[test]
    fn import_aligned_positions_transposes_to_long() {
        let conn = fresh_conn();
        let f = write_tsv(
            "query\tTM2_D51\tTM3_C78\tTM3_E102\tfract\ttaxpath\n\
             seq1\tN\tR\tE\t1.0\tEukaryota;Animalia\n\
             seq2\tK\tC\tD\t0.8\tBacteria\n",
        );
        let n = import_aligned_positions(&conn, "rp", f.path()).unwrap();
        // 2 queries x 3 positions
        assert_eq!(n, 6);
        assert_eq!(count(&conn, "aligned_positions"), 6);

        let (residues, idx): (String, i32) = conn
            .query_row(
                "SELECT residues, pos_index FROM aligned_positions WHERE query_name='seq1' AND pos_label='TM3_C78'",
                [],
                |r| Ok((r.get(0)?, r.get(1)?)),
            )
            .unwrap();
        assert_eq!(residues, "R");
        assert_eq!(idx, 1); // TM3_C78 is the 2nd position column (0-based index 1)

        // pos_index recovers the original wide-table column order.
        let order: Vec<String> = {
            let mut stmt = conn
                .prepare("SELECT DISTINCT pos_label FROM aligned_positions ORDER BY pos_index")
                .unwrap();
            let rows = stmt
                .query_map([], |r| r.get::<_, String>(0))
                .unwrap()
                .map(|r| r.unwrap())
                .collect();
            rows
        };
        assert_eq!(order, vec!["TM2_D51", "TM3_C78", "TM3_E102"]);
    }

    // --- import_refpkg_meta ---

    #[test]
    fn import_refpkg_meta_loads_identity_and_provenance() {
        let conn = fresh_conn();
        let json = r#"{
            "name": "subACV",
            "refpkg": "/db/refpkg/rhodopsin/subACV",
            "hmmlen": "229",
            "hmmname": "rhodopsin.subACV",
            "fhmmO": "/db/refpkg/rhodopsin/subACV/x.hmm",
            "falnO": "/db/refpkg/rhodopsin/subACV/x.mfa",
            "ftreO": "/db/refpkg/rhodopsin/subACV/x.tree",
            "faln": "/db/refpkg/rhodopsin/subACV/derived/backbone.mfa",
            "ppdir": "/db/refpkg/rhodopsin/subACV/derived/for_pplacer"
        }"#;
        let f = tempfile::Builder::new().suffix(".json").tempfile().unwrap();
        std::fs::write(f.path(), json).unwrap();

        import_refpkg_meta(&conn, "subACV", f.path()).unwrap();
        assert_eq!(count(&conn, "refpkgs"), 1);

        let (rp, dir, name, len, aln, tree, hmm): (
            String,
            String,
            String,
            i64,
            String,
            String,
            String,
        ) = conn
            .query_row(
                "SELECT refpkg, refpkg_dir, hmmname, hmmlen, aln_source, tree_source, hmm_source FROM refpkgs",
                [],
                |r| {
                    Ok((
                        r.get(0)?,
                        r.get(1)?,
                        r.get(2)?,
                        r.get(3)?,
                        r.get(4)?,
                        r.get(5)?,
                        r.get(6)?,
                    ))
                },
            )
            .unwrap();
        assert_eq!(rp, "subACV");
        assert_eq!(dir, "/db/refpkg/rhodopsin/subACV");
        assert_eq!(name, "rhodopsin.subACV");
        assert_eq!(len, 229);
        assert_eq!(aln, "/db/refpkg/rhodopsin/subACV/x.mfa");
        assert_eq!(tree, "/db/refpkg/rhodopsin/subACV/x.tree");
        assert_eq!(hmm, "/db/refpkg/rhodopsin/subACV/x.hmm");
    }

    #[test]
    fn import_run_json_loads_params_and_software() {
        let conn = fresh_conn();
        let json = r#"{
            "params": {"minhmmlen":"174","evalue":"1e-5","maxseqlen":"100000","placer":"pplacer"},
            "software": [
                {"name":"hmmsearch","path":"/env/bin/hmmsearch","version":"HMMER 3.4"},
                {"name":"pplacer","path":"/env/bin/pplacer","version":"ff7556d-dirty"}
            ]
        }"#;
        let f = tempfile::Builder::new().suffix(".json").tempfile().unwrap();
        std::fs::write(f.path(), json).unwrap();

        let (np, ns) = import_run_json(&conn, "subACV", f.path()).unwrap();
        assert_eq!(np, 4);
        assert_eq!(ns, 2);
        assert_eq!(count(&conn, "run_params"), 4);
        assert_eq!(count(&conn, "software"), 2);

        let v: String = conn
            .query_row(
                "SELECT value FROM run_params WHERE param='minhmmlen' AND refpkg='subACV'",
                [],
                |r| r.get(0),
            )
            .unwrap();
        assert_eq!(v, "174");

        let (path, ver): (String, String) = conn
            .query_row(
                "SELECT path, version FROM software WHERE name='pplacer'",
                [],
                |r| Ok((r.get(0)?, r.get(1)?)),
            )
            .unwrap();
        assert_eq!(path, "/env/bin/pplacer");
        assert_eq!(ver, "ff7556d-dirty");
    }

    #[test]
    fn import_aligned_positions_handles_na_taxpath() {
        let conn = fresh_conn();
        let f = write_tsv("query\tpos1\tfract\ttaxpath\nseq1\tA\tNA\tNA\n");
        let n = import_aligned_positions(&conn, "rp", f.path()).unwrap();
        assert_eq!(n, 1);

        let (fract, taxpath): (Option<f64>, Option<String>) = conn
            .query_row("SELECT fract, taxpath FROM aligned_positions", [], |r| {
                Ok((r.get(0)?, r.get(1)?))
            })
            .unwrap();
        assert_eq!(fract, None);
        assert_eq!(taxpath, None);
    }
}
