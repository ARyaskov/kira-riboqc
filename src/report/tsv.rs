use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::pipeline::stage_translation_regime::StageTranslationRegimeOutput;
use crate::report::pipeline_contract::PipelineCellRow;

pub const PIPELINE_TSV_HEADER: &str = "barcode\tsample\tcondition\tspecies\tlibsize\tnnz\texpressed_genes\ttranslation_load\tribosome_density\telongation_pressure\tinitiation_bias\tribosomal_specialization\tstress_translation_index\tribosome_core\tinitiation_core\tbio_core\tmtor_core\tisr_core\tTPI\tRBL\tmTOR_P\tISR_A\tTPIB\tTSM\ttranslation_high\tbiogenesis_high\tisr_active\tproteotoxic_risk\ttranslational_stress_mode\tregime\tflags\tconfidence";

pub fn write_riboqc_tsv(path: &Path, rows: &[PipelineCellRow]) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    w.write_all(PIPELINE_TSV_HEADER.as_bytes())?;
    w.write_all(b"\n")?;

    let mut ibuf = itoa::Buffer::new();

    for r in rows {
        w.write_all(r.barcode.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(r.sample.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(r.condition.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(r.species.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(ibuf.format(r.libsize).as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(ibuf.format(r.nnz).as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(ibuf.format(r.expressed_genes).as_bytes())?;
        w.write_all(b"\t")?;
        write_f6(&mut w, r.translation_load)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, r.ribosome_density)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, r.elongation_pressure)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, r.initiation_bias)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, r.ribosomal_specialization)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, r.stress_translation_index)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.ribosome_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.initiation_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.bio_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.mtor_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.isr_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.tpi)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.rbl)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.mtor_p)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.isr_a)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.tpib)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, r.tsm)?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(r.translation_high))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(r.biogenesis_high))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(r.isr_active))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(r.proteotoxic_risk))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(r.translational_stress_mode))?;
        w.write_all(b"\t")?;
        w.write_all(r.regime.as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(r.flags.as_bytes())?;
        w.write_all(b"\t")?;
        write_f6(&mut w, r.confidence)?;
        w.write_all(b"\n")?;
    }

    Ok(())
}

pub const TRANSLATION_REGIME_HEADER: &str = "cell_id\tribosome_loading_heterogeneity\ttranslation_selectivity_index\tisr_like_signature_score\tcodon_bias_proxy\ttranslation_commitment_score\tribosome_core\tinitiation_core\tbio_core\tmtor_core\tisr_core\tTPI\tRBL\tmTOR_P\tISR_A\tTPIB\tTSM\ttranslation_high\tbiogenesis_high\tisr_active\tproteotoxic_risk\ttranslational_stress_mode\tmissing_panel_gene_count\ttranslation_regime";

pub fn write_translation_regime_tsv(
    path: &Path,
    stage: &StageTranslationRegimeOutput,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    w.write_all(TRANSLATION_REGIME_HEADER.as_bytes())?;
    w.write_all(b"\n")?;

    let mut ibuf = itoa::Buffer::new();

    for cell in &stage.cells {
        w.write_all(cell.cell_id.as_bytes())?;
        w.write_all(b"\t")?;
        write_f6(&mut w, cell.ribosome_loading_heterogeneity)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, cell.translation_selectivity_index)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, cell.isr_like_signature_score)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, cell.codon_bias_proxy)?;
        w.write_all(b"\t")?;
        write_f6(&mut w, cell.translation_commitment_score)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.ribosome_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.initiation_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.bio_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.mtor_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.isr_core)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.tpi)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.rbl)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.mtor_p)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.isr_a)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.tpib)?;
        w.write_all(b"\t")?;
        write_f6_nan(&mut w, cell.tsm)?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(cell.translation_high))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(cell.biogenesis_high))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(cell.isr_active))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(cell.proteotoxic_risk))?;
        w.write_all(b"\t")?;
        w.write_all(bool_byte(cell.translational_stress_mode))?;
        w.write_all(b"\t")?;
        w.write_all(ibuf.format(cell.missing_panel_gene_count).as_bytes())?;
        w.write_all(b"\t")?;
        w.write_all(cell.translation_regime.as_bytes())?;
        w.write_all(b"\n")?;
    }

    Ok(())
}

#[inline]
fn bool_byte(flag: bool) -> &'static [u8] {
    if flag { b"1" } else { b"0" }
}

fn write_f6<W: Write>(w: &mut W, v: f64) -> std::io::Result<()> {
    if v.is_nan() {
        return w.write_all(b"NaN");
    }
    let mut buf = ryu::Buffer::new();
    let rounded = (v * 1_000_000.0).round() / 1_000_000.0;
    let s = buf.format(rounded);
    write_fixed6(w, s)
}

fn write_f6_nan<W: Write>(w: &mut W, v: f64) -> std::io::Result<()> {
    write_f6(w, v)
}

fn write_fixed6<W: Write>(w: &mut W, s: &str) -> std::io::Result<()> {
    if let Some(dot) = s.find('.') {
        let frac_len = s.len() - dot - 1;
        if frac_len >= 6 {
            return w.write_all(&s.as_bytes()[..dot + 7]);
        }
        w.write_all(s.as_bytes())?;
        for _ in frac_len..6 {
            w.write_all(b"0")?;
        }
        Ok(())
    } else {
        w.write_all(s.as_bytes())?;
        w.write_all(b".000000")
    }
}
