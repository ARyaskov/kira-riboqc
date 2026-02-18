use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::pipeline::stage_translation_regime::StageTranslationRegimeOutput;
use crate::report::pipeline_contract::PipelineCellRow;

pub const PIPELINE_TSV_HEADER: &str = "barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\ttranslation_load\tribosome_density\telongation_pressure\tinitiation_bias\tribosomal_specialization\tstress_translation_index\tregime\tflags\tconfidence";

pub fn write_riboqc_tsv(path: &Path, rows: &[PipelineCellRow]) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(8 * 1024 * 1024, file);

    writeln!(w, "{PIPELINE_TSV_HEADER}")?;

    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}",
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
            r.regime,
            r.flags,
            r.confidence,
        )?;
    }

    Ok(())
}

pub const TRANSLATION_REGIME_HEADER: &str = "cell_id\tribosome_loading_heterogeneity\ttranslation_selectivity_index\tisr_like_signature_score\tcodon_bias_proxy\ttranslation_commitment_score\ttranslation_regime";

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
            "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}",
            cell.cell_id,
            cell.ribosome_loading_heterogeneity,
            cell.translation_selectivity_index,
            cell.isr_like_signature_score,
            cell.codon_bias_proxy,
            cell.translation_commitment_score,
            cell.translation_regime
        )?;
    }

    Ok(())
}
