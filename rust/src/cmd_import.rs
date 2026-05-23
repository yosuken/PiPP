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

    /// Overwrite an existing DB file
    #[arg(long)]
    pub overwrite: bool,
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

    eprintln!("wrote {}", db_path.display());
    Ok(())
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
            app.append_row(params![refpkg, query, label, residues, fract, taxpath])?;
            n += 1;
        }
    }
    app.flush()?;
    Ok(n)
}
