use clap::{Parser, Subcommand};

pub mod run;

pub use run::{Mode, RunArgs, RunMode};

#[derive(Parser, Debug)]
#[command(name = "kira-riboqc")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    Run(RunArgs),
}
