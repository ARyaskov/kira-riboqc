use tracing::info;

use crate::model::scores::{fragility_from_tpc, lfti, ras, tss};
use crate::pipeline::stage2_axes::Stage2Output;

pub struct CellScores {
    pub tss: f64,
    pub lfti: f64,
    pub ras: f64,
    pub ras_red_flag: bool,
    pub tss_explain: String,
    pub lfti_explain: String,
    pub ras_explain: String,
}

pub struct Stage3Output {
    pub scores: Vec<CellScores>,
}

pub fn run_stage3(stage2: &Stage2Output) -> anyhow::Result<Stage3Output> {
    let _span = tracing::info_span!("stage3_scores").entered();
    info!("Stage 3 start");

    let mut scores = Vec::with_capacity(stage2.axes.len());
    let mut red_flags = 0u32;

    for ax in &stage2.axes {
        let tl = ax.tl;
        let st = ax.st;
        let rqc = ax.rqc;
        let tpc = ax.tpc;
        let fragility = fragility_from_tpc(tpc);

        let tss_val = tss(tl, rqc, st);
        let lfti_val = lfti(rqc, st, fragility);
        let ras_val = ras(tl, rqc, fragility);
        let ras_red_flag = tl > 0.70 && rqc > 0.70 && tpc < -0.20;

        if ras_red_flag {
            red_flags += 1;
        }

        let tss_explain = explain(vec![
            ("TL", 0.45 * tl),
            ("RQC", 0.35 * rqc),
            ("ST", 0.20 * st),
        ]);
        let lfti_explain = explain(vec![
            ("RQC", 0.40 * rqc),
            ("ST", 0.35 * st),
            ("fragility", 0.25 * fragility),
        ]);
        let ras_explain = explain(vec![
            ("TL", 0.50 * tl),
            ("RQC", 0.30 * rqc),
            ("fragility", 0.20 * fragility),
        ]);

        scores.push(CellScores {
            tss: tss_val,
            lfti: lfti_val,
            ras: ras_val,
            ras_red_flag,
            tss_explain,
            lfti_explain,
            ras_explain,
        });
    }

    let tss_median = median_non_nan(scores.iter().map(|s| s.tss));
    let lfti_median = median_non_nan(scores.iter().map(|s| s.lfti));
    let ras_median = median_non_nan(scores.iter().map(|s| s.ras));
    let frac_red_flag = if scores.is_empty() {
        0.0
    } else {
        red_flags as f64 / scores.len() as f64
    };

    info!(
        tss_median = tss_median,
        lfti_median = lfti_median,
        ras_median = ras_median,
        ras_red_flag_frac = frac_red_flag,
        "Stage 3 summary"
    );
    info!("Stage 3 end");

    Ok(Stage3Output { scores })
}

pub fn explain(components: Vec<(&str, f64)>) -> String {
    let mut parts: Vec<(&str, f64, f64)> = components
        .into_iter()
        .map(|(name, value)| (name, value, value.abs()))
        .collect();

    parts.sort_by(|a, b| {
        b.2.partial_cmp(&a.2)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(b.0))
    });

    let mut out = String::new();
    for (idx, (name, value, _abs)) in parts.iter().enumerate() {
        if idx > 0 {
            out.push(' ');
        }
        let sign = if *value < 0.0 { '-' } else { '+' };
        out.push(sign);
        out.push_str(name);
        out.push('(');
        out.push_str(&format!("{:.2}", value.abs()));
        out.push(')');
    }
    out
}

fn median_non_nan<I>(iter: I) -> f64
where
    I: Iterator<Item = f64>,
{
    let mut values: Vec<f64> = iter.filter(|v| !v.is_nan()).collect();
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = values.len() / 2;
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        (values[mid - 1] + values[mid]) / 2.0
    }
}
