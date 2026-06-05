use std::path::Path;

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;
use kira_scio::model::CanonicalData;

use super::InputError;
use super::gene_index::FeatureRow;

#[derive(Debug, Clone)]
pub struct CscMatrix {
    pub n_rows: u32,
    pub n_cols: u32,
    pub col_ptr: Vec<u32>,
    pub row_idx: Vec<u32>,
    pub values: Vec<u32>,
}

pub struct MtxDataset {
    pub matrix: CscMatrix,
    pub barcodes: Vec<String>,
    pub features: Vec<FeatureRow>,
}

pub fn load_mtx_dataset(input_path: &Path) -> Result<MtxDataset, InputError> {
    let canonical = Reader::with_options(
        input_path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_all()
    .map_err(|e| InputError::Parse(e.message))?;

    canonical_into_dataset(canonical)
}

fn canonical_into_dataset(canonical: CanonicalData) -> Result<MtxDataset, InputError> {
    let CanonicalData { metadata, matrix } = canonical;

    let n_rows = u32::try_from(matrix.n_genes)
        .map_err(|_| InputError::Parse("n_rows exceeds u32".to_string()))?;
    let n_cols = u32::try_from(matrix.n_cells)
        .map_err(|_| InputError::Parse("n_cols exceeds u32".to_string()))?;

    let mut col_ptr = Vec::with_capacity(matrix.col_ptr.len());
    for v in &matrix.col_ptr {
        col_ptr.push(
            u32::try_from(*v).map_err(|_| InputError::Parse("col_ptr exceeds u32".to_string()))?,
        );
    }

    let row_idx = matrix.row_idx;

    let mut values = Vec::with_capacity(matrix.values.len());
    for v in matrix.values {
        if v < 0.0 {
            return Err(InputError::Parse(
                "matrix values must be non-negative".to_string(),
            ));
        }
        if (v.fract()).abs() > 1e-6 {
            return Err(InputError::Parse(
                "matrix values must be integer-like for riboqc".to_string(),
            ));
        }
        values.push(v as u32);
    }

    let csc = CscMatrix {
        n_rows,
        n_cols,
        col_ptr,
        row_idx,
        values,
    };

    let mut features = Vec::with_capacity(metadata.gene_symbols.len());
    for (idx, symbol) in metadata.gene_symbols.iter().enumerate() {
        let raw_id = metadata
            .gene_ids
            .get(idx)
            .cloned()
            .unwrap_or_else(|| format!("GENE_{}", idx + 1));
        features.push(FeatureRow::from_raw(raw_id, symbol.clone()));
    }

    Ok(MtxDataset {
        matrix: csc,
        barcodes: metadata.barcodes,
        features,
    })
}
