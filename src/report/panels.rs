use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::input::{InputBundle, SharedCacheData};
use crate::model::axes::{
    CHAPERONE, LUXURY, MACHINERY, METAB, RQC_MACH, RRNA_PROXY, STRESS, UBIQUITIN, is_ribosomal,
};

struct PanelDef {
    panel_id: &'static str,
    panel_name: &'static str,
    panel_group: &'static str,
    symbols_defined: &'static [&'static str],
    ribosomal_pattern: bool,
}

const PANEL_DEFS: &[PanelDef] = &[
    PanelDef {
        panel_id: "RP",
        panel_name: "Ribosomal proteins",
        panel_group: "translation_core",
        symbols_defined: &[],
        ribosomal_pattern: true,
    },
    PanelDef {
        panel_id: "MACH",
        panel_name: "Translation machinery",
        panel_group: "translation_core",
        symbols_defined: MACHINERY,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "RRNA",
        panel_name: "rRNA processing",
        panel_group: "translation_core",
        symbols_defined: RRNA_PROXY,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "STRESS",
        panel_name: "Stress translation",
        panel_group: "stress",
        symbols_defined: STRESS,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "METAB",
        panel_name: "Metabolic survival",
        panel_group: "stress",
        symbols_defined: METAB,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "LUXURY",
        panel_name: "Luxury anabolic",
        panel_group: "growth",
        symbols_defined: LUXURY,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "RQC",
        panel_name: "RQC machinery",
        panel_group: "rqc",
        symbols_defined: RQC_MACH,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "CHAPERONE",
        panel_name: "Chaperone",
        panel_group: "proteostasis",
        symbols_defined: CHAPERONE,
        ribosomal_pattern: false,
    },
    PanelDef {
        panel_id: "UBIQUITIN",
        panel_name: "Ubiquitin stress",
        panel_group: "proteostasis",
        symbols_defined: UBIQUITIN,
        ribosomal_pattern: false,
    },
];

struct PanelStats {
    panel_size_defined: usize,
    panel_size_mappable: usize,
    missing_genes: String,
    coverage_values: Vec<f64>,
    sum_values: Vec<f64>,
}

enum MatrixSource<'a> {
    Matrix {
        col_ptr: &'a [u32],
        row_idx: &'a [u32],
        values: &'a [u32],
        n_cols: usize,
    },
    Cache(&'a SharedCacheData),
}

impl<'a> MatrixSource<'a> {
    fn from_input(input: &'a InputBundle) -> Self {
        if let Some(cache) = input.shared_cache.as_ref() {
            Self::Cache(cache)
        } else {
            Self::Matrix {
                col_ptr: &input.matrix.col_ptr,
                row_idx: &input.matrix.row_idx,
                values: &input.matrix.values,
                n_cols: input.matrix.n_cols as usize,
            }
        }
    }

    fn n_cols(&self) -> usize {
        match self {
            MatrixSource::Matrix { n_cols, .. } => *n_cols,
            MatrixSource::Cache(cache) => cache.n_cells as usize,
        }
    }

    fn column_bounds(&self, col: usize) -> (usize, usize) {
        match self {
            MatrixSource::Matrix { col_ptr, .. } => {
                (col_ptr[col] as usize, col_ptr[col + 1] as usize)
            }
            MatrixSource::Cache(cache) => {
                (cache.col_ptr[col] as usize, cache.col_ptr[col + 1] as usize)
            }
        }
    }

    fn row_idx(&self) -> &[u32] {
        match self {
            MatrixSource::Matrix { row_idx, .. } => row_idx,
            MatrixSource::Cache(cache) => cache.row_idx.as_slice(),
        }
    }

    fn values(&self) -> &[u32] {
        match self {
            MatrixSource::Matrix { values, .. } => values,
            MatrixSource::Cache(cache) => cache.values_u32.as_slice(),
        }
    }
}

pub fn write_panels_report(path: &Path, input: &InputBundle) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(8 * 1024 * 1024, file);

    writeln!(
        w,
        "panel_id\tpanel_name\tpanel_group\tpanel_size_defined\tpanel_size_mappable\tmissing_genes\tcoverage_median\tcoverage_p10\tsum_median\tsum_p90\tsum_p99"
    )?;

    let source = MatrixSource::from_input(input);
    let n_cells = source.n_cols();
    let n_genes = input.gene_index.genes.len();

    let mut panel_stats = Vec::with_capacity(PANEL_DEFS.len());
    let mut gene_panel_mask = vec![0u16; n_genes];

    for (panel_idx, def) in PANEL_DEFS.iter().enumerate() {
        let (gene_ids, panel_size_defined, missing_genes) = resolve_panel(input, def);
        for gid in &gene_ids {
            let idx = *gid as usize;
            if idx < gene_panel_mask.len() {
                gene_panel_mask[idx] |= 1u16 << panel_idx;
            }
        }

        panel_stats.push(PanelStats {
            panel_size_defined,
            panel_size_mappable: gene_ids.len(),
            missing_genes,
            coverage_values: vec![0.0; n_cells],
            sum_values: vec![0.0; n_cells],
        });
    }

    let row_idx = source.row_idx();
    let values = source.values();

    for col in 0..n_cells {
        let (start, end) = source.column_bounds(col);
        let mut hit_counts = vec![0u32; PANEL_DEFS.len()];

        for k in start..end {
            let row = row_idx[k] as usize;
            let gid = input.gene_index.row_to_gene[row] as usize;
            let mask = *gene_panel_mask.get(gid).unwrap_or(&0);
            if mask == 0 {
                continue;
            }
            for panel_idx in 0..PANEL_DEFS.len() {
                if (mask & (1u16 << panel_idx)) != 0 {
                    hit_counts[panel_idx] += 1;
                    panel_stats[panel_idx].sum_values[col] += values[k] as f64;
                }
            }
        }

        for panel_idx in 0..PANEL_DEFS.len() {
            let mappable = panel_stats[panel_idx].panel_size_mappable;
            panel_stats[panel_idx].coverage_values[col] = if mappable == 0 {
                0.0
            } else {
                hit_counts[panel_idx] as f64 / mappable as f64
            };
        }
    }

    for (panel_idx, def) in PANEL_DEFS.iter().enumerate() {
        let stats = &panel_stats[panel_idx];
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            def.panel_id,
            def.panel_name,
            def.panel_group,
            stats.panel_size_defined,
            stats.panel_size_mappable,
            stats.missing_genes,
            percentile(&stats.coverage_values, 0.50),
            percentile(&stats.coverage_values, 0.10),
            percentile(&stats.sum_values, 0.50),
            percentile(&stats.sum_values, 0.90),
            percentile(&stats.sum_values, 0.99)
        )?;
    }

    Ok(())
}

fn resolve_panel(input: &InputBundle, def: &PanelDef) -> (Vec<u32>, usize, String) {
    if def.ribosomal_pattern {
        let mut ids = Vec::new();
        for g in &input.gene_index.genes {
            if is_ribosomal(&g.symbol) {
                ids.push(g.gene_id);
            }
        }
        let len = ids.len();
        return (ids, len, String::new());
    }

    let mut ids = Vec::new();
    let mut missing = Vec::new();

    for symbol in def.symbols_defined {
        if let Some(gid) = input.gene_index.map.get(*symbol) {
            ids.push(*gid);
        } else {
            missing.push((*symbol).to_string());
        }
    }

    (
        ids,
        def.symbols_defined.len(),
        if missing.is_empty() {
            String::new()
        } else {
            missing.join(",")
        },
    )
}

fn percentile(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let rank = (p * sorted.len() as f64).ceil() as usize;
    let idx = if rank == 0 { 0 } else { rank - 1 };
    sorted[idx.min(sorted.len() - 1)]
}
