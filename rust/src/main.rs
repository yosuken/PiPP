use anyhow::Result;
use clap::{Parser, Subcommand};

mod cmd_import;
mod cmd_parse_hmmsearch;
mod schema;

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
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Import(args) => cmd_import::run(args),
        Cmd::ParseHmmsearch(args) => cmd_parse_hmmsearch::run(args),
    }
}
