use anyhow::Result;
use clap::{Parser, Subcommand};

mod cmd_import;
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
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Import(args) => cmd_import::run(args),
    }
}
