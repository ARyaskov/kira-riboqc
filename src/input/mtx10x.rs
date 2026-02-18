use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::GzDecoder;

use super::InputError;

#[derive(Debug, Clone)]
pub struct CscMatrix {
    pub n_rows: u32,
    pub n_cols: u32,
    pub col_ptr: Vec<u32>,
    pub row_idx: Vec<u32>,
    pub values: Vec<u32>,
}

#[derive(Debug, Clone)]
struct Triplet {
    col: u32,
    row: u32,
    val: u32,
    idx: u32,
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
    let reader = BufReader::new(open_maybe_gz(path)?);
    let mut barcodes = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }
        barcodes.push(trimmed.to_string());
    }
    Ok(barcodes)
}

pub fn load_mtx(path: &Path) -> Result<CscMatrix, InputError> {
    let reader = BufReader::new(open_maybe_gz(path)?);
    let mut lines = reader.lines();

    let header_line = loop {
        match lines.next() {
            Some(Ok(line)) => {
                let trimmed = line.trim();
                if trimmed.is_empty() || trimmed.starts_with('%') {
                    continue;
                }
                break trimmed.to_string();
            }
            Some(Err(e)) => return Err(InputError::Io(e)),
            None => return Err(InputError::Parse("matrix file is empty".to_string())),
        }
    };

    let dims: Vec<&str> = header_line.split_whitespace().collect();
    if dims.len() != 3 {
        return Err(InputError::Parse(
            "matrix header line must have 3 integers".to_string(),
        ));
    }
    let n_rows: u32 = dims[0]
        .parse()
        .map_err(|_| InputError::Parse("invalid n_rows".to_string()))?;
    let n_cols: u32 = dims[1]
        .parse()
        .map_err(|_| InputError::Parse("invalid n_cols".to_string()))?;
    let nnz_expected: u32 = dims[2]
        .parse()
        .map_err(|_| InputError::Parse("invalid nnz".to_string()))?;

    let mut triplets = Vec::with_capacity(nnz_expected as usize);
    let mut idx: u32 = 0;

    for line in lines {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('%') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 3 {
            return Err(InputError::Parse("invalid matrix entry line".to_string()));
        }
        let row_raw: i64 = parts[0]
            .parse()
            .map_err(|_| InputError::Parse("invalid row index".to_string()))?;
        let col_raw: i64 = parts[1]
            .parse()
            .map_err(|_| InputError::Parse("invalid col index".to_string()))?;
        let val_raw: i64 = parts[2]
            .parse()
            .map_err(|_| InputError::Parse("invalid value".to_string()))?;
        if row_raw < 1 || col_raw < 1 {
            return Err(InputError::Parse(
                "matrix indices must be 1-based".to_string(),
            ));
        }
        if val_raw < 0 {
            return Err(InputError::Parse(
                "matrix values must be non-negative".to_string(),
            ));
        }

        let row = (row_raw - 1) as u32;
        let col = (col_raw - 1) as u32;
        if row >= n_rows || col >= n_cols {
            return Err(InputError::Parse(
                "matrix indices out of bounds".to_string(),
            ));
        }

        triplets.push(Triplet {
            col,
            row,
            val: val_raw as u32,
            idx,
        });
        idx = idx.wrapping_add(1);
    }

    if triplets.len() != nnz_expected as usize {
        return Err(InputError::Parse(format!(
            "nnz mismatch: expected {nnz_expected}, got {}",
            triplets.len()
        )));
    }

    triplets.sort_by(|a, b| {
        a.col
            .cmp(&b.col)
            .then(a.row.cmp(&b.row))
            .then(a.idx.cmp(&b.idx))
    });

    let mut counts = vec![0u32; n_cols as usize];
    for t in &triplets {
        counts[t.col as usize] += 1;
    }
    let mut col_ptr = vec![0u32; n_cols as usize + 1];
    for i in 0..n_cols as usize {
        col_ptr[i + 1] = col_ptr[i] + counts[i];
    }

    let mut row_idx = Vec::with_capacity(triplets.len());
    let mut values = Vec::with_capacity(triplets.len());
    for t in triplets {
        row_idx.push(t.row);
        values.push(t.val);
    }

    Ok(CscMatrix {
        n_rows,
        n_cols,
        col_ptr,
        row_idx,
        values,
    })
}
