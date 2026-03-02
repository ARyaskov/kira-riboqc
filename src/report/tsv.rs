use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::pipeline::stage_translation_regime::StageTranslationRegimeOutput;
use crate::report::pipeline_contract::PipelineCellRow;

pub const PIPELINE_TSV_HEADER: &str = "barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\ttranslation_load\tribosome_density\telongation_pressure\tinitiation_bias\tribosomal_specialization\tstress_translation_index\tribosome_core\tinitiation_core\tbio_core\tmtor_core\tisr_core\tTPI\tRBL\tmTOR_P\tISR_A\tTPIB\tTSM\ttranslation_high\tbiogenesis_high\tisr_active\tproteotoxic_risk\ttranslational_stress_mode\tregime\tflags\tconfidence";

pub fn write_riboqc_tsv(path: &Path, rows: &[PipelineCellRow]) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(8 * 1024 * 1024, file);

    writeln!(w, "{PIPELINE_TSV_HEADER}")?;

    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
            r.barcode,
            r.sample,
            r.condition,
            r.species,
            r.libsize,
            r.nnz,
            r.expressed_genes,
            r.translation_load,
            r.ribosome_density,
            r.elongation_pressure,
            r.initiation_bias,
            r.ribosomal_specialization,
            r.stress_translation_index,
            format_nan6(r.ribosome_core),
            format_nan6(r.initiation_core),
            format_nan6(r.bio_core),
            format_nan6(r.mtor_core),
            format_nan6(r.isr_core),
            format_nan6(r.tpi),
            format_nan6(r.rbl),
            format_nan6(r.mtor_p),
            format_nan6(r.isr_a),
            format_nan6(r.tpib),
            format_nan6(r.tsm),
            bool_num(r.translation_high),
            bool_num(r.biogenesis_high),
            bool_num(r.isr_active),
            bool_num(r.proteotoxic_risk),
            bool_num(r.translational_stress_mode),
            r.regime,
            r.flags,
            r.confidence,
        )?;
    }

    Ok(())
}

pub const TRANSLATION_REGIME_HEADER: &str = "cell_id\tribosome_loading_heterogeneity\ttranslation_selectivity_index\tisr_like_signature_score\tcodon_bias_proxy\ttranslation_commitment_score\tribosome_core\tinitiation_core\tbio_core\tmtor_core\tisr_core\tTPI\tRBL\tmTOR_P\tISR_A\tTPIB\tTSM\ttranslation_high\tbiogenesis_high\tisr_active\tproteotoxic_risk\ttranslational_stress_mode\tmissing_panel_gene_count\ttranslation_regime";

pub fn write_translation_regime_tsv(
    path: &Path,
    stage: &StageTranslationRegimeOutput,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, file);

    writeln!(w, "{TRANSLATION_REGIME_HEADER}")?;
    for cell in &stage.cells {
        writeln!(
            w,
            "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            cell.cell_id,
            cell.ribosome_loading_heterogeneity,
            cell.translation_selectivity_index,
            cell.isr_like_signature_score,
            cell.codon_bias_proxy,
            cell.translation_commitment_score,
            format_nan6(cell.ribosome_core),
            format_nan6(cell.initiation_core),
            format_nan6(cell.bio_core),
            format_nan6(cell.mtor_core),
            format_nan6(cell.isr_core),
            format_nan6(cell.tpi),
            format_nan6(cell.rbl),
            format_nan6(cell.mtor_p),
            format_nan6(cell.isr_a),
            format_nan6(cell.tpib),
            format_nan6(cell.tsm),
            bool_num(cell.translation_high),
            bool_num(cell.biogenesis_high),
            bool_num(cell.isr_active),
            bool_num(cell.proteotoxic_risk),
            bool_num(cell.translational_stress_mode),
            cell.missing_panel_gene_count,
            cell.translation_regime
        )?;
    }

    Ok(())
}

fn bool_num(flag: bool) -> u8 {
    if flag { 1 } else { 0 }
}

fn format_nan6(value: f64) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else {
        format!("{value:.6}")
    }
}
