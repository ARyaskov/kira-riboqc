use std::collections::BTreeSet;
use std::fs;
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
    prefix: Option<&str>,
) -> Result<InputFormat, InputError> {
    let matrix_name = with_prefix(prefix, "matrix.mtx");
    let barcodes_name = with_prefix(prefix, "barcodes.tsv");
    let features_name = with_prefix(prefix, "features.tsv");
    let genes_name = with_prefix(prefix, "genes.tsv");

    let matrix_path = detect_file(input_path, &matrix_name)
        .ok_or_else(|| InputError::MissingFile(format!("{matrix_name}(.gz)")))?;
    let barcodes_path = detect_file(input_path, &barcodes_name)
        .ok_or_else(|| InputError::MissingFile(format!("{barcodes_name}(.gz)")))?;
    let features_path = detect_file(input_path, &features_name)
        .or_else(|| detect_file(input_path, &genes_name))
        .ok_or_else(|| {
            InputError::MissingFile(format!("{features_name}(.gz) or {genes_name}(.gz)"))
        })?;

    Ok(InputFormat::TenXDir {
        matrix_path,
        barcodes_path,
        features_path,
    })
}

pub fn detect_prefix(input_path: &Path) -> Result<Option<String>, InputError> {
    let mut prefixes: BTreeSet<String> = BTreeSet::new();
    let suffixes = [
        "_matrix.mtx",
        "_matrix.mtx.gz",
        "_features.tsv",
        "_features.tsv.gz",
        "_genes.tsv",
        "_genes.tsv.gz",
        "_barcodes.tsv",
        "_barcodes.tsv.gz",
    ];

    for entry in fs::read_dir(input_path)? {
        let entry = entry?;
        let file_name = entry.file_name();
        let Some(file_name) = file_name.to_str() else {
            continue;
        };

        for suffix in suffixes {
            if let Some(prefix) = file_name.strip_suffix(suffix) {
                if !prefix.is_empty() {
                    prefixes.insert(prefix.to_string());
                }
            }
        }
    }

    if prefixes.is_empty() {
        return Ok(None);
    }
    if prefixes.len() == 1 {
        return Ok(prefixes.into_iter().next());
    }
    Err(InputError::InvalidFormat(
        "multiple prefixed datasets found in input directory".to_string(),
    ))
}

pub fn resolve_cache_path(input_path: &Path, prefix: Option<&str>) -> PathBuf {
    match prefix {
        Some(p) => input_path.join(format!("{p}.kira-organelle.bin")),
        None => input_path.join("kira-organelle.bin"),
    }
}

fn with_prefix(prefix: Option<&str>, base: &str) -> String {
    match prefix {
        Some(p) => format!("{p}_{base}"),
        None => base.to_string(),
    }
}

fn detect_file(dir: &Path, filename: &str) -> Option<PathBuf> {
    let direct = dir.join(filename);
    if direct.exists() {
        return Some(direct);
    }
    let gz = dir.join(format!("{filename}.gz"));
    if gz.exists() {
        return Some(gz);
    }
    None
}
