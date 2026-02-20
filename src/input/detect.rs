use std::path::{Path, PathBuf};

use super::InputError;
use super::types::InputFormat;

pub fn detect_input(input_path: &Path) -> Result<InputFormat, InputError> {
    if !input_path.exists() {
        return Err(InputError::InvalidFormat(format!(
            "input path does not exist: {}",
            input_path.display()
        )));
    }
    if !input_path.is_dir() {
        return Err(InputError::InvalidFormat(format!(
            "input path is not a directory: {}",
            input_path.display()
        )));
    }

    detect_input_for_prefix(input_path, None)
}

pub fn detect_input_for_prefix(
    input_path: &Path,
    _prefix: Option<&str>,
) -> Result<InputFormat, InputError> {
    let ds = kira_scio::discover(input_path).map_err(|e| InputError::InvalidFormat(e.message))?;

    let barcodes_path = ds
        .barcodes
        .ok_or_else(|| InputError::MissingFile("barcodes.tsv(.gz)".to_string()))?;
    let features_path = ds.features.or(ds.genes).ok_or_else(|| {
        InputError::MissingFile("features.tsv(.gz) or genes.tsv(.gz)".to_string())
    })?;

    Ok(InputFormat::TenXDir {
        matrix_path: ds.matrix,
        barcodes_path,
        features_path,
    })
}

pub fn detect_prefix(input_path: &Path) -> Result<Option<String>, InputError> {
    kira_scio::detect_prefix(input_path).map_err(|e| InputError::InvalidFormat(e.to_string()))
}

pub fn resolve_cache_path(input_path: &Path, prefix: Option<&str>) -> PathBuf {
    input_path.join(kira_scio::resolve_shared_cache_filename(prefix))
}
