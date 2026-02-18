use clap::{Args, ValueEnum};
use std::path::PathBuf;

#[derive(Args, Debug, Clone)]
pub struct RunArgs {
    #[arg(long)]
    pub input: PathBuf,
    #[arg(long)]
    pub out: PathBuf,
    #[arg(long, value_enum, default_value_t = Mode::Cell)]
    pub mode: Mode,
    #[arg(long)]
    pub metadata: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = RunMode::Standalone)]
    pub run_mode: RunMode,
    #[arg(long, default_value_t = false)]
    pub disable_translation_extension: bool,
}

#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    Cell,
    Sample,
}

#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Standalone,
    Pipeline,
}
