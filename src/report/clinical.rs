use crate::report::summary::Summary;

pub fn build_report(summary: &Summary) -> String {
    let mut out = String::new();
    out.push_str("kira-riboqc - Translation State Quality Control Report\n");
    out.push_str(&format!("Cells analyzed: {}\n", summary.input.n_cells));
    out.push_str(&format!("Species: {}\n\n", summary.input.species));

    out.push_str("Dominant translation regimes:\n");
    let mut pairs = summary
        .regimes
        .fractions
        .iter()
        .map(|(k, v)| (k.clone(), *v))
        .collect::<Vec<_>>();
    pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then_with(|| a.0.cmp(&b.0)));
    for (regime, frac) in pairs {
        out.push_str(&format!("- {}: {:.1}%\n", regime, frac * 100.0));
    }

    out.push_str("\nDistributions:\n");
    out.push_str(&format!(
        "- translation_load median {:.3}, p90 {:.3}, p99 {:.3}\n",
        summary.distributions.translation_load.median,
        summary.distributions.translation_load.p90,
        summary.distributions.translation_load.p99
    ));
    out.push_str(&format!(
        "- ribosome_density median {:.3}, p90 {:.3}, p99 {:.3}\n",
        summary.distributions.ribosome_density.median,
        summary.distributions.ribosome_density.p90,
        summary.distributions.ribosome_density.p99
    ));
    out.push_str(&format!(
        "- stress_translation_index median {:.3}, p90 {:.3}, p99 {:.3}\n",
        summary.distributions.stress_translation_index.median,
        summary.distributions.stress_translation_index.p90,
        summary.distributions.stress_translation_index.p99
    ));

    out
}
