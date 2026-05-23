//! Library crate for `pipp_util`.
//!
//! The binary entry point is in `src/main.rs`. Subcommand implementations
//! and helpers live here so they can be unit/integration tested without
//! shelling out to the binary.

pub mod cmd_import;
pub mod cmd_parse_hmmsearch;
pub mod schema;
