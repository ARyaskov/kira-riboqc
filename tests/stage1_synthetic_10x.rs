use std::fs::{create_dir_all, write};

use kira_riboqc::cli::{Mode, RunArgs, RunMode};
use kira_riboqc::pipeline::stage1_load::run_stage1;
use tempfile::tempdir;

#[test]
fn synthetic_10x_loads_and_indexes() {
    let tmp = tempdir().unwrap();
    let input_dir = tmp.path().join("input");
    let out_dir = tmp.path().join("out");
    create_dir_all(&input_dir).unwrap();

    let matrix = "%%MatrixMarket matrix coordinate integer general\n% comment\n3 2 4\n2 1 3\n1 2 2\n1 1 5\n3 2 1\n";
    let barcodes = "CELL1\nCELL2\n";
    let features =
        "G1\tRPLP0\tGene Expression\nG2\tRPLP0\tGene Expression\nG3\tEEF2\tGene Expression\n";

    write(input_dir.join("matrix.mtx"), matrix).unwrap();
    write(input_dir.join("barcodes.tsv"), barcodes).unwrap();
    write(input_dir.join("features.tsv"), features).unwrap();

    let args = RunArgs {
        input: input_dir.clone(),
        out: out_dir,
        mode: Mode::Cell,
        metadata: None,
        run_mode: RunMode::Standalone,
        disable_translation_extension: false,
    };

    let (bundle, _stats) = run_stage1(&args).unwrap();

    assert_eq!(bundle.matrix.n_rows, 3);
    assert_eq!(bundle.matrix.n_cols, 2);
    assert_eq!(bundle.gene_index.genes.len(), 2);
    assert_eq!(bundle.gene_index.row_to_gene.len(), 3);
    assert_eq!(
        bundle.gene_index.row_to_gene[0],
        bundle.gene_index.row_to_gene[1]
    );
    assert_eq!(bundle.gene_index.duplicates.len(), 1);

    let dup = &bundle.gene_index.duplicates[0];
    assert_eq!(dup.symbol, "RPLP0");
    assert_eq!(dup.first_row, 0);
    assert_eq!(dup.dup_row, 1);

    assert_eq!(bundle.matrix.col_ptr, vec![0, 2, 4]);
    assert_eq!(bundle.matrix.row_idx, vec![0, 1, 0, 2]);
    assert_eq!(bundle.matrix.values, vec![5, 3, 2, 1]);
}

#[test]
fn metadata_join_filters_to_barcodes() {
    let tmp = tempdir().unwrap();
    let input_dir = tmp.path().join("input");
    let out_dir = tmp.path().join("out");
    create_dir_all(&input_dir).unwrap();

    let matrix = "%%MatrixMarket matrix coordinate integer general\n3 2 1\n1 1 7\n";
    let barcodes = "CELL1\nCELL2\n";
    let features = "G1\tRPLP0\nG2\tEEF2\nG3\tGAPDH\n";
    let metadata = "cell_id\tgroup\nCELL1\tA\nEXTRA\tB\n";

    write(input_dir.join("matrix.mtx"), matrix).unwrap();
    write(input_dir.join("barcodes.tsv"), barcodes).unwrap();
    write(input_dir.join("features.tsv"), features).unwrap();

    let meta_path = input_dir.join("meta.tsv");
    write(&meta_path, metadata).unwrap();

    let args = RunArgs {
        input: input_dir,
        out: out_dir,
        mode: Mode::Cell,
        metadata: Some(meta_path),
        run_mode: RunMode::Standalone,
        disable_translation_extension: false,
    };

    let (bundle, _stats) = run_stage1(&args).unwrap();
    let table = bundle.metadata.expect("metadata should be loaded");
    assert_eq!(table.rows.len(), 1);
    let row = table.rows.get("CELL1").expect("CELL1 metadata");
    assert_eq!(row.cell_id, "CELL1");
    assert_eq!(row.fields.get("group").map(String::as_str), Some("A"));
}
