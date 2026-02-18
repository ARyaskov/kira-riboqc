use std::collections::BTreeMap;
use std::path::Path;

use csv::ReaderBuilder;

use super::InputError;
use super::mtx10x::open_maybe_gz;

#[derive(Debug, Clone)]
pub struct CellMeta {
    pub cell_id: String,
    pub fields: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
pub struct MetadataTable {
    pub rows: BTreeMap<String, CellMeta>,
}

pub fn load_metadata(path: &Path) -> Result<MetadataTable, InputError> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let (delim, fallback) = match ext.as_str() {
        "tsv" => (b'\t', None),
        "csv" => (b',', None),
        _ => (b'\t', Some(b',')),
    };

    match load_metadata_with_delim(path, delim) {
        Ok(table) => Ok(table),
        Err(err) => {
            if let Some(fallback_delim) = fallback {
                let table = load_metadata_with_delim(path, fallback_delim)?;
                Ok(table)
            } else {
                Err(err)
            }
        }
    }
}

fn load_metadata_with_delim(path: &Path, delim: u8) -> Result<MetadataTable, InputError> {
    let reader = open_maybe_gz(path)?;
    let mut rdr = ReaderBuilder::new().delimiter(delim).from_reader(reader);

    let headers = rdr
        .headers()
        .map_err(|e| InputError::Metadata(format!("failed to read headers: {e}")))?
        .clone();
    if headers.is_empty() {
        return Err(InputError::Metadata(
            "metadata headers are empty".to_string(),
        ));
    }

    let mut cell_id_idx: Option<usize> = None;
    let mut barcode_idx: Option<usize> = None;
    for (i, h) in headers.iter().enumerate() {
        let key = h.trim().to_ascii_lowercase();
        if key == "cell_id" {
            cell_id_idx = Some(i);
        }
        if key == "barcode" {
            barcode_idx = Some(i);
        }
    }

    let id_idx = cell_id_idx.or(barcode_idx).ok_or_else(|| {
        InputError::Metadata("metadata must have a cell_id or barcode column".to_string())
    })?;

    let mut rows: BTreeMap<String, CellMeta> = BTreeMap::new();
    for result in rdr.records() {
        let record =
            result.map_err(|e| InputError::Metadata(format!("failed to read record: {e}")))?;
        let cell_id = record.get(id_idx).unwrap_or("").to_string();
        if cell_id.is_empty() {
            continue;
        }

        if rows.contains_key(&cell_id) {
            continue;
        }

        let mut fields = BTreeMap::new();
        for (i, header) in headers.iter().enumerate() {
            if i == id_idx {
                continue;
            }
            let value = record.get(i).unwrap_or("");
            fields.insert(header.to_string(), value.to_string());
        }

        rows.insert(cell_id.clone(), CellMeta { cell_id, fields });
    }

    Ok(MetadataTable { rows })
}
