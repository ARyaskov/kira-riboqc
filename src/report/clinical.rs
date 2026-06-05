use std::fmt::Write;

use crate::report::summary::Summary;

pub fn build_report(summary: &Summary) -> String {
    let mut out = String::with_capacity(512);
    out.push_str("kira-riboqc - Translation State Quality Control Report\n");
    let _ = writeln!(out, "Cells analyzed: {}", summary.input.n_cells);
    let _ = writeln!(out, "Species: {}\n", summary.input.species);

    out.push_str("Dominant translation regimes:\n");
    let mut pairs: Vec<(&str, f64)> = summary
        .regimes
        .fractions
        .iter()
        .map(|(k, v)| (k.as_str(), *v))
        .collect();
    pairs.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(b.0))
    });
    for (regime, frac) in pairs {
        let _ = writeln!(out, "- {}: {:.1}%", regime, frac * 100.0);
    }

    out.push_str("\nDistributions:\n");
    let _ = writeln!(
        out,
        "- translation_load median {:.3}, p90 {:.3}, p99 {:.3}",
        summary.distributions.translation_load.median,
        summary.distributions.translation_load.p90,
        summary.distributions.translation_load.p99
    );
    let _ = writeln!(
        out,
        "- ribosome_density median {:.3}, p90 {:.3}, p99 {:.3}",
        summary.distributions.ribosome_density.median,
        summary.distributions.ribosome_density.p90,
        summary.distributions.ribosome_density.p99
    );
    let _ = writeln!(
        out,
        "- stress_translation_index median {:.3}, p90 {:.3}, p99 {:.3}",
        summary.distributions.stress_translation_index.median,
        summary.distributions.stress_translation_index.p90,
        summary.distributions.stress_translation_index.p99
    );

    out
}
