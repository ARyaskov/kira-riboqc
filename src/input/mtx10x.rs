use std::fs::File;
use std::io::Read;
use std::path::Path;

use flate2::read::GzDecoder;
use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;

use super::InputError;

#[derive(Debug, Clone)]
pub struct CscMatrix {
    pub n_rows: u32,
    pub n_cols: u32,
    pub col_ptr: Vec<u32>,
    pub row_idx: Vec<u32>,
    pub values: Vec<u32>,
}

pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>, InputError> {
    let file = File::open(path)?;
    if path.extension().map(|ext| ext == "gz").unwrap_or(false) {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

pub fn load_barcodes(path: &Path) -> Result<Vec<String>, InputError> {
    let md = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_metadata()
    .map_err(|e| InputError::Parse(e.message))?;
    Ok(md.barcodes)
}

pub fn load_mtx(path: &Path) -> Result<CscMatrix, InputError> {
    let matrix = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::Mtx10x),
            strict: true,
        },
    )
    .read_matrix()
    .map_err(|e| InputError::Parse(e.message))?;

    let n_rows = u32::try_from(matrix.n_genes)
        .map_err(|_| InputError::Parse("n_rows exceeds u32".to_string()))?;
    let n_cols = u32::try_from(matrix.n_cells)
        .map_err(|_| InputError::Parse("n_cols exceeds u32".to_string()))?;

    let col_ptr = matrix
        .col_ptr
        .into_iter()
        .map(|v| u32::try_from(v).map_err(|_| InputError::Parse("col_ptr exceeds u32".to_string())))
        .collect::<Result<Vec<_>, _>>()?;
    let row_idx = matrix
        .row_idx
        .into_iter()
        .map(|v| u32::try_from(v).map_err(|_| InputError::Parse("row_idx exceeds u32".to_string())))
        .collect::<Result<Vec<_>, _>>()?;

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

    Ok(CscMatrix {
        n_rows,
        n_cols,
        col_ptr,
        row_idx,
        values,
    })
}
