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
            refpkg     TEXT   NOT NULL,
            query_name TEXT   NOT NULL,
            pos_label  TEXT   NOT NULL,
            residues   TEXT   NOT NULL,
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
        "#,
    )?;
    Ok(())
}
