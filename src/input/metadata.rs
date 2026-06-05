use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use rustc_hash::FxHashMap;

use super::InputError;

#[derive(Debug, Clone)]
pub struct CellMeta {
    pub cell_id: String,
    pub fields: FxHashMap<String, String>,
}

impl CellMeta {
    pub fn field(&self, key: &str) -> Option<&str> {
        self.fields
            .get(&key.to_ascii_lowercase())
            .map(String::as_str)
    }
}

#[derive(Debug, Clone, Default)]
pub struct MetadataTable {
    pub rows: FxHashMap<String, CellMeta>,
    pub header_keys_lower: Vec<String>,
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
                load_metadata_with_delim(path, fallback_delim)
            } else {
                Err(err)
            }
        }
    }
}

fn open_buffered(path: &Path) -> Result<BufReader<Box<dyn Read>>, InputError> {
    let file = File::open(path)?;
    let raw: Box<dyn Read> = if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(BufReader::with_capacity(64 * 1024, raw))
}

fn load_metadata_with_delim(path: &Path, delim: u8) -> Result<MetadataTable, InputError> {
    let reader = open_buffered(path)?;
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

    let header_keys_lower: Vec<String> = headers
        .iter()
        .map(|h| h.trim().to_ascii_lowercase())
        .collect();

    let id_idx = header_keys_lower
        .iter()
        .position(|k| k == "cell_id")
        .or_else(|| header_keys_lower.iter().position(|k| k == "barcode"))
        .ok_or_else(|| {
            InputError::Metadata("metadata must have a cell_id or barcode column".to_string())
        })?;

    let mut rows: FxHashMap<String, CellMeta> = FxHashMap::default();
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

        let mut fields = FxHashMap::default();
        for (i, key_lower) in header_keys_lower.iter().enumerate() {
            if i == id_idx {
                continue;
            }
            let value = record.get(i).unwrap_or("");
            fields.insert(key_lower.clone(), value.to_string());
        }

        rows.insert(cell_id.clone(), CellMeta { cell_id, fields });
    }

    Ok(MetadataTable {
        rows,
        header_keys_lower,
    })
}
