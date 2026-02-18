use clap::Parser;
use std::time::Instant;
use tracing::info;

use kira_riboqc::cli::{Cli, Commands};
use kira_riboqc::pipeline::stage_translation_regime::run_stage_translation_regime;
use kira_riboqc::pipeline::stage1_load::run_stage1;
use kira_riboqc::pipeline::stage2_axes::run_stage2;
use kira_riboqc::pipeline::stage3_scores::run_stage3;
use kira_riboqc::pipeline::stage4_classify::run_stage4;
use kira_riboqc::pipeline::stage5_report::run_stage5;

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();
    tracing_subscriber::fmt().with_target(false).try_init().ok();
    info!("SIMD backend: {}", kira_riboqc::simd::SIMD_KIND);

    match cli.command {
        Commands::Run(args) => {
            let start = Instant::now();
            info!("Stage 1 start");
            let (input, stage1_stats) = run_stage1(&args)?;
            let stage2 = run_stage2(&input)?;
            let stage3 = run_stage3(&stage2)?;
            let stage4 = run_stage4(&stage2, &stage3)?;
            let stage_translation = if args.disable_translation_extension {
                None
            } else {
                Some(run_stage_translation_regime(&input, &stage2)?)
            };
            run_stage5(
                &args.out,
                args.run_mode,
                &input,
                &stage1_stats,
                &stage2,
                &stage3,
                &stage4,
                stage_translation.as_ref(),
            )?;
            let elapsed = start.elapsed();
            info!(?elapsed, "Run end");
        }
    }

    Ok(())
}
