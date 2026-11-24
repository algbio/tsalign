#![deny(clippy::mod_module_files)]

use anyhow::Result;
use clap::Parser;

mod align;
mod alphabet;
mod preprocess;
mod show;
mod util;

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
    let cli = Cli::parse();

    match cli.subcommand {
        Subcommand::Align(cli) => align::cli(*cli),
        Subcommand::Preprocess(cli) => preprocess::cli(cli),
        Subcommand::Show(cli) => show::cli(cli),
    }
}
