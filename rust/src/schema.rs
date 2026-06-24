use anyhow::Result;
use duckdb::Connection;

pub fn create_all(conn: &Connection) -> Result<()> {
    conn.execute_batch(
        r#"
        CREATE TABLE IF NOT EXISTS assignments (
            refpkg     TEXT   NOT NULL,
            query_name TEXT   NOT NULL,
            LWR        DOUBLE NOT NULL,
            fract      DOUBLE NOT NULL,
            aLWR       DOUBLE NOT NULL,
            afract     DOUBLE NOT NULL,
            taxopath   TEXT   NOT NULL
        );

        CREATE TABLE IF NOT EXISTS aa_features (
            refpkg         TEXT NOT NULL,
            query_name     TEXT NOT NULL,
            len            INTEGER,
            len_std_aa     INTEGER,
            avg_mw         DOUBLE,
            n_arsc         DOUBLE,
            c_arsc         DOUBLE,
            s_arsc         DOUBLE,
            aa_K           INTEGER,
            aa_R           INTEGER,
            aa_H           INTEGER,
            aa_D           INTEGER,
            aa_E           INTEGER,
            aa_N           INTEGER,
            aa_Q           INTEGER,
            aa_S           INTEGER,
            aa_T           INTEGER,
            aa_Y           INTEGER,
            aa_A           INTEGER,
            aa_V           INTEGER,
            aa_L           INTEGER,
            aa_I           INTEGER,
            aa_P           INTEGER,
            aa_F           INTEGER,
            aa_M           INTEGER,
            aa_W           INTEGER,
            aa_G           INTEGER,
            aa_C           INTEGER,
            aa_others      INTEGER
        );

        CREATE TABLE IF NOT EXISTS aligned_positions (
            refpkg     TEXT    NOT NULL,
            query_name TEXT    NOT NULL,
            pos_index  INTEGER NOT NULL,  -- 0-based column order in the source TSV (recovers wide-table layout)
            pos_label  TEXT    NOT NULL,
            residues   TEXT    NOT NULL,
            fract      DOUBLE,
            taxpath    TEXT
        );

        -- records when `pipp_util clamp-jplace` actually sanitized a placement
        -- jplace (non-finite values and/or negative branch lengths). A row
        -- here means the avoidance fix was used for that placement.
        CREATE TABLE IF NOT EXISTS jplace_clamps (
            refpkg       TEXT   NOT NULL,
            jplace       TEXT   NOT NULL,
            n_nonfinite  BIGINT NOT NULL,
            n_neg_branch BIGINT NOT NULL
        );

        -- full jplace JSON produced by the selected placer, one row per
        -- result/<refpkg>/placement/*.jplace file.
        CREATE TABLE IF NOT EXISTS jplaces (
            refpkg  TEXT NOT NULL,
            jplace  TEXT NOT NULL,
            content TEXT NOT NULL
        );

        -- task-level completion state captured from log/tasks/<task>/exit and
        -- GNU parallel joblogs while they still exist before cleanup.
        CREATE TABLE IF NOT EXISTS task_runs (
            refpkg        TEXT NOT NULL,
            task          TEXT NOT NULL,
            status        TEXT NOT NULL,
            exit_code     INTEGER,
            n_jobs        BIGINT,
            started_epoch DOUBLE,
            runtime_sec   DOUBLE,
            finished_at   TEXT,
            exit_marker   TEXT,
            parallel_log  TEXT
        );

        -- one row per refpkg: identity + provenance taken from the
        -- validate_refpkg metadata (backbone.json). The source_* columns point
        -- at the original refpkg inputs (not the run-local derived cache).
        CREATE TABLE IF NOT EXISTS refpkgs (
            refpkg        TEXT NOT NULL,
            refpkg_dir    TEXT,
            hmmname       TEXT,
            hmmlen        BIGINT,
            aln_source    TEXT,
            tree_source   TEXT,
            hmm_source    TEXT
        );

        -- all run-time options as key/value (so new CLI flags need no schema
        -- change). One row per (refpkg, param); params also include metadata
        -- like pipp_version, run_datetime, command_line, query.
        CREATE TABLE IF NOT EXISTS run_params (
            refpkg TEXT NOT NULL,
            param  TEXT NOT NULL,
            value  TEXT
        );

        -- resolved path + version of each external tool used by the run.
        CREATE TABLE IF NOT EXISTS software (
            refpkg  TEXT NOT NULL,
            name    TEXT NOT NULL,
            path    TEXT,
            version TEXT
        );

        -- prefilter best hits (hmmsearch parsed: non-overlapping, linked-merged,
        -- filter-passing). One row per detected region for THIS refpkg's HMM
        -- (best-hit.tsv rows where hmm_name == the refpkg's HMM name).
        CREATE TABLE IF NOT EXISTS prefilter_hits (
            refpkg       TEXT NOT NULL,
            region_name  TEXT,
            protein      TEXT,
            protein_len  BIGINT,
            protein_info TEXT,
            hmm_name     TEXT,
            hmm_acc      TEXT,
            hmm_desc     TEXT,
            hmm_len      BIGINT,
            score        DOUBLE,
            bias         DOUBLE,
            c_evalue     DOUBLE,
            i_evalue     DOUBLE,
            hmm_fm       BIGINT,
            hmm_to       BIGINT,
            ali_fm       BIGINT,
            ali_to       BIGINT,
            env_fm       BIGINT,
            env_to       BIGINT,
            acc          DOUBLE,
            full_evalue  DOUBLE,
            full_score   DOUBLE,
            link         TEXT
        );

        -- whole (ungapped) query protein sequence, one row per detected protein
        -- (result/<refpkg>/seq/whole.fa). region sequences are NOT stored; they
        -- are a slice of the whole sequence (see the query_region view).
        CREATE TABLE IF NOT EXISTS query_whole (
            refpkg     TEXT NOT NULL,
            query_name TEXT NOT NULL,
            sequence   TEXT NOT NULL
        );

        -- aligned (gapped) query sequence in backbone column space, one row per
        -- placed region (result/<refpkg>/alignment/aligned_wo_ref.fa). The full
        -- alignment (aligned.fa) is this plus the refpkg backbone alignment.
        CREATE TABLE IF NOT EXISTS query_aligned (
            refpkg      TEXT NOT NULL,
            region_name TEXT NOT NULL,
            sequence    TEXT NOT NULL
        );

        -- region (ungapped) sequence = the whole protein sliced at the best-hit
        -- alignment coordinates. Derived, so seq/region.fa needs no storage.
        CREATE VIEW IF NOT EXISTS query_region AS
            SELECT h.refpkg,
                   h.region_name,
                   h.protein AS query_name,
                   substr(w.sequence, h.ali_fm, h.ali_to - h.ali_fm + 1) AS sequence
            FROM prefilter_hits h
            JOIN query_whole w
              ON w.refpkg = h.refpkg AND w.query_name = h.protein;

        -- per-sequence best i-Evalue against this refpkg's HMM (evalues.tsv),
        -- only cells that actually had a value (sparse matrix -> long format).
        CREATE TABLE IF NOT EXISTS prefilter_evalues (
            refpkg   TEXT NOT NULL,
            seq      TEXT NOT NULL,
            i_evalue DOUBLE
        );
        "#,
    )?;
    Ok(())
}
