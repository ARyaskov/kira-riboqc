use std::fs;

use rustc_hash::FxHashMap;
use tracing::{info, warn};

use crate::cli::{RunArgs, RunMode};
use crate::input::{
    CellMeta, CscMatrix, FeatureRow, InputBundle, InputError, InputFormat, MetadataTable,
    MtxDataset, SharedCacheData, Stage1Stats, build_gene_index, detect_input,
    detect_input_for_prefix, detect_prefix, load_metadata, load_mtx_dataset, normalize_symbol,
    read_shared_cache, resolve_cache_path,
};

pub fn run_stage1(args: &RunArgs) -> anyhow::Result<(InputBundle, Stage1Stats)> {
    let _span = tracing::info_span!("stage1_load").entered();

    fs::create_dir_all(&args.out)?;

    let (format, matrix, barcodes, features, gene_index, shared_cache) = match args.run_mode {
        RunMode::Standalone => load_standalone(args)?,
        RunMode::Pipeline => load_pipeline(args)?,
    };

    validate_dimensions(
        &matrix,
        barcodes.len(),
        features.len(),
        shared_cache.as_ref(),
    )?;

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
            .map(|c| c.nnz() as u32)
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

type LoadOutput = (
    InputFormat,
    CscMatrix,
    Vec<String>,
    Vec<FeatureRow>,
    crate::input::GeneIndex,
    Option<SharedCacheData>,
);

fn load_standalone(args: &RunArgs) -> anyhow::Result<LoadOutput> {
    let format = detect_input(&args.input).map_err(anyhow::Error::from)?;
    let InputFormat::TenXDir { .. } = &format else {
        return Err(anyhow::anyhow!(
            "unexpected cache format in standalone mode"
        ));
    };
    let MtxDataset {
        matrix,
        barcodes,
        features,
    } = load_mtx_dataset(&args.input).map_err(anyhow::Error::from)?;
    let gene_index = build_gene_index(&features);
    Ok((format, matrix, barcodes, features, gene_index, None))
}

fn load_pipeline(args: &RunArgs) -> anyhow::Result<LoadOutput> {
    let prefix = detect_prefix(&args.input).map_err(anyhow::Error::from)?;
    let cache_path = resolve_cache_path(&args.input, prefix.as_deref());

    if cache_path.exists() {
        let cache = read_shared_cache(&cache_path).map_err(anyhow::Error::from)?;
        let barcodes = cache.barcodes().to_vec();
        let features = cache
            .genes()
            .iter()
            .enumerate()
            .map(|(idx, symbol)| FeatureRow {
                raw_id: format!("CACHE_GENE_{}", idx + 1),
                raw_name: symbol.clone(),
                raw_type: String::new(),
                norm_symbol: normalize_symbol(symbol),
            })
            .collect::<Vec<_>>();
        let gene_index = build_gene_index(&features);
        let matrix = csc_placeholder(&cache)?;
        let format = InputFormat::SharedCache {
            cache_path: cache_path.clone(),
        };
        Ok((format, matrix, barcodes, features, gene_index, Some(cache)))
    } else {
        warn!(
            cache_path = %cache_path.display(),
            "shared cache not found; falling back to MTX input"
        );
        let format =
            detect_input_for_prefix(&args.input, prefix.as_deref()).map_err(anyhow::Error::from)?;
        let InputFormat::TenXDir { .. } = &format else {
            return Err(anyhow::anyhow!(
                "unexpected cache format in pipeline MTX fallback"
            ));
        };
        let MtxDataset {
            matrix,
            barcodes,
            features,
        } = load_mtx_dataset(&args.input).map_err(anyhow::Error::from)?;
        let gene_index = build_gene_index(&features);
        Ok((format, matrix, barcodes, features, gene_index, None))
    }
}

fn csc_placeholder(cache: &SharedCacheData) -> anyhow::Result<CscMatrix> {
    let n_rows =
        u32::try_from(cache.n_genes()).map_err(|_| anyhow::anyhow!("cache n_genes exceeds u32"))?;
    let n_cols =
        u32::try_from(cache.n_cells()).map_err(|_| anyhow::anyhow!("cache n_cells exceeds u32"))?;
    Ok(CscMatrix {
        n_rows,
        n_cols,
        col_ptr: vec![0; n_cols as usize + 1],
        row_idx: Vec::new(),
        values: Vec::new(),
    })
}

fn validate_dimensions(
    matrix: &CscMatrix,
    n_barcodes: usize,
    n_features: usize,
    cache: Option<&SharedCacheData>,
) -> Result<(), InputError> {
    let (n_cols, n_rows) = match cache {
        Some(c) => (c.n_cells(), c.n_genes()),
        None => (matrix.n_cols as usize, matrix.n_rows as usize),
    };
    if n_cols != n_barcodes {
        return Err(InputError::Dimension(format!(
            "n_cols mismatch: matrix has {n_cols}, barcodes has {n_barcodes}"
        )));
    }
    if n_rows != n_features {
        return Err(InputError::Dimension(format!(
            "n_rows mismatch: matrix has {n_rows}, features has {n_features}"
        )));
    }
    Ok(())
}

fn filter_metadata_to_barcodes(table: MetadataTable, barcodes: &[String]) -> MetadataTable {
    let mut rows: FxHashMap<String, CellMeta> = FxHashMap::default();
    for barcode in barcodes {
        if let Some(row) = table.rows.get(barcode) {
            rows.insert(barcode.clone(), row.clone());
        }
    }
    MetadataTable {
        rows,
        header_keys_lower: table.header_keys_lower,
    }
}
