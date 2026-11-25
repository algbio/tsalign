#![deny(clippy::mod_module_files)]

use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};

mod align;
mod alphabet;
mod preprocess;
mod show;
mod util;

/// jemalloc makes the classic aligner 15-23 % faster than glibc malloc,
/// because it allocates and frees millions of small search nodes (#278, #288).
/// It does not support MSVC targets, so Windows builds keep the system allocator.
#[cfg(not(target_os = "windows"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[derive(Parser)]

struct Cli {
    #[command(subcommand)]
    subcommand: Subcommand,
}

#[derive(clap::Subcommand)]
enum Subcommand {
    Align(Box<align::Cli>),
    Preprocess(preprocess::Cli),
    Show(show::Cli),
}

fn main() -> Result<()> {
    // The argument matches are kept to tell explicitly given arguments from default values.
    let matches = Cli::command().get_matches();
    let cli = Cli::from_arg_matches(&matches).unwrap_or_else(|error| error.exit());
    let subcommand_matches = matches
        .subcommand()
        .map(|(_, matches)| matches)
        .expect("a subcommand is required");

    match cli.subcommand {
        Subcommand::Align(cli) => align::cli(*cli, subcommand_matches),
        Subcommand::Preprocess(cli) => preprocess::cli(cli, subcommand_matches),
        Subcommand::Show(cli) => show::cli(cli),
    }
}
