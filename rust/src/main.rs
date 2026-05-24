use anyhow::Result;
use clap::{Parser, Subcommand};
use pipp_util::{cmd_clamp_jplace, cmd_import, cmd_parse_hmmsearch, cmd_validate_query};

#[derive(Parser)]
#[command(name = "pipp_util", version, about = "Utilities for PiPP")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Import TSV outputs from one refpkg result directory into a DuckDB file
    Import(cmd_import::Args),
    /// Parse an hmmsearch output and emit best-hit / all-hit tables and fastas
    #[command(name = "parse-hmmsearch")]
    ParseHmmsearch(cmd_parse_hmmsearch::Args),
    /// Validate + clean a query FASTA (single pass, streaming MD5)
    #[command(name = "validate-query")]
    ValidateQuery(cmd_validate_query::Args),
    /// Clamp negative branch lengths in a jplace tree to 0 (apples-2 fix)
    #[command(name = "clamp-jplace")]
    ClampJplace(cmd_clamp_jplace::Args),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Import(args) => cmd_import::run(args),
        Cmd::ParseHmmsearch(args) => cmd_parse_hmmsearch::run(args),
        Cmd::ValidateQuery(args) => cmd_validate_query::run(args),
        Cmd::ClampJplace(args) => cmd_clamp_jplace::run(args),
    }
}
