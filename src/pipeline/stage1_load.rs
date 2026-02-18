use std::collections::BTreeMap;
use std::fs;

use tracing::{info, warn};

use crate::cli::{RunArgs, RunMode};
use crate::input::{
    InputBundle, InputError, InputFormat, MetadataTable, Stage1Stats, build_gene_index,
    detect_input, detect_input_for_prefix, detect_prefix, load_barcodes, load_features,
    load_metadata, load_mtx, read_shared_cache, resolve_cache_path,
};

pub fn run_stage1(args: &RunArgs) -> anyhow::Result<(InputBundle, Stage1Stats)> {
    let _span = tracing::info_span!("stage1_load").entered();

    fs::create_dir_all(&args.out)?;

    let (format, matrix, barcodes, features, gene_index, shared_cache) = match args.run_mode {
        RunMode::Standalone => {
            let format = detect_input(&args.input).map_err(anyhow::Error::from)?;
            let (matrix_path, barcodes_path, features_path) = match &format {
                InputFormat::TenXDir {
                    matrix_path,
                    barcodes_path,
                    features_path,
                } => (matrix_path, barcodes_path, features_path),
                InputFormat::SharedCache { .. } => {
                    return Err(anyhow::anyhow!(
                        "unexpected cache format in standalone mode"
                    ));
                }
            };
            let barcodes = load_barcodes(barcodes_path).map_err(anyhow::Error::from)?;
            let features = load_features(features_path).map_err(anyhow::Error::from)?;
            let gene_index = build_gene_index(&features);
            let matrix = load_mtx(matrix_path).map_err(anyhow::Error::from)?;
            (format, matrix, barcodes, features, gene_index, None)
        }
        RunMode::Pipeline => {
            let prefix = detect_prefix(&args.input).map_err(anyhow::Error::from)?;
            let cache_path = resolve_cache_path(&args.input, prefix.as_deref());
            if cache_path.exists() {
                let cache = read_shared_cache(&cache_path).map_err(anyhow::Error::from)?;
                let barcodes = cache.barcodes.clone();
                let features = cache
                    .genes
                    .iter()
                    .enumerate()
                    .map(|(idx, symbol)| crate::input::FeatureRow {
                        raw_id: format!("CACHE_GENE_{}", idx + 1),
                        raw_name: symbol.clone(),
                        raw_type: String::new(),
                        norm_symbol: crate::input::normalize_symbol(symbol),
                    })
                    .collect::<Vec<_>>();
                let gene_index = build_gene_index(&features);
                let matrix = csc_from_cache(&cache)?;
                let format = InputFormat::SharedCache {
                    cache_path: cache_path.clone(),
                };
                (format, matrix, barcodes, features, gene_index, Some(cache))
            } else {
                warn!(
                    cache_path = %cache_path.display(),
                    "shared cache not found; falling back to MTX input"
                );
                let format = detect_input_for_prefix(&args.input, prefix.as_deref())
                    .map_err(anyhow::Error::from)?;
                let (matrix_path, barcodes_path, features_path) = match &format {
                    InputFormat::TenXDir {
                        matrix_path,
                        barcodes_path,
                        features_path,
                    } => (matrix_path, barcodes_path, features_path),
                    InputFormat::SharedCache { .. } => {
                        return Err(anyhow::anyhow!(
                            "unexpected cache format in pipeline MTX fallback"
                        ));
                    }
                };
                let barcodes = load_barcodes(barcodes_path).map_err(anyhow::Error::from)?;
                let features = load_features(features_path).map_err(anyhow::Error::from)?;
                let gene_index = build_gene_index(&features);
                let matrix = load_mtx(matrix_path).map_err(anyhow::Error::from)?;
                (format, matrix, barcodes, features, gene_index, None)
            }
        }
    };

    validate_dimensions(&matrix, barcodes.len(), features.len())?;

    let metadata = if let Some(path) = &args.metadata {
        let table = load_metadata(path).map_err(anyhow::Error::from)?;
        Some(filter_metadata_to_barcodes(table, &barcodes))
    } else {
        None
    };

    let stats = Stage1Stats {
        n_genes_raw: features.len() as u32,
        n_genes_unique: gene_index.genes.len() as u32,
        n_cells: barcodes.len() as u32,
        nnz: shared_cache
            .as_ref()
            .map(|c| c.nnz as u32)
            .unwrap_or(matrix.values.len() as u32),
    };

    let format_label = match &format {
        InputFormat::TenXDir { .. } => "TenXDir",
        InputFormat::SharedCache { .. } => "SharedCache",
    };

    info!(
        format = format_label,
        n_genes_raw = stats.n_genes_raw,
        n_genes_unique = stats.n_genes_unique,
        n_cells = stats.n_cells,
        nnz = stats.nnz,
        duplicate_genes = gene_index.duplicates.len(),
        "Stage 1 summary"
    );

    let bundle = InputBundle {
        format,
        matrix,
        barcodes,
        features,
        gene_index,
        metadata,
        shared_cache,
    };

    Ok((bundle, stats))
}

fn csc_from_cache(
    cache: &crate::input::SharedCacheData,
) -> anyhow::Result<crate::input::CscMatrix> {
    let n_rows =
        u32::try_from(cache.n_genes).map_err(|_| anyhow::anyhow!("cache n_genes exceeds u32"))?;
    let n_cols =
        u32::try_from(cache.n_cells).map_err(|_| anyhow::anyhow!("cache n_cells exceeds u32"))?;
    Ok(crate::input::CscMatrix {
        n_rows,
        n_cols,
        // Cache-backed runs read sparse arrays from shared_cache directly in Stage 2+.
        // Keep only lightweight shape placeholders here to avoid duplicating large buffers.
        col_ptr: vec![0; n_cols as usize + 1],
        row_idx: Vec::new(),
        values: Vec::new(),
    })
}

fn validate_dimensions(
    matrix: &crate::input::CscMatrix,
    n_barcodes: usize,
    n_features: usize,
) -> Result<(), InputError> {
    if matrix.n_cols as usize != n_barcodes {
        return Err(InputError::Dimension(format!(
            "n_cols mismatch: matrix has {}, barcodes has {}",
            matrix.n_cols, n_barcodes
        )));
    }
    if matrix.n_rows as usize != n_features {
        return Err(InputError::Dimension(format!(
            "n_rows mismatch: matrix has {}, features has {}",
            matrix.n_rows, n_features
        )));
    }
    Ok(())
}

fn filter_metadata_to_barcodes(table: MetadataTable, barcodes: &[String]) -> MetadataTable {
    let mut rows: BTreeMap<String, crate::input::CellMeta> = BTreeMap::new();
    for barcode in barcodes {
        if let Some(row) = table.rows.get(barcode) {
            rows.insert(barcode.clone(), row.clone());
        }
    }
    MetadataTable { rows }
}
