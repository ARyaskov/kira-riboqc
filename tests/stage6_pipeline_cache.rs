use std::fs::{create_dir_all, write};
use std::path::Path;

use clap::Parser;
use kira_riboqc::cli::{Cli, Commands, Mode, RunArgs, RunMode};
use kira_riboqc::input::{InputFormat, detect_prefix, read_shared_cache, resolve_cache_path};
use kira_riboqc::pipeline::stage1_load::run_stage1;
use tempfile::tempdir;

#[test]
fn cli_run_mode_parsing() {
    let cli = Cli::try_parse_from(["kira-riboqc", "run", "--input", "in", "--out", "out"])
        .expect("cli parse");
    let Commands::Run(args) = cli.command;
    assert_eq!(args.run_mode, RunMode::Standalone);

    let cli = Cli::try_parse_from([
        "kira-riboqc",
        "run",
        "--input",
        "in",
        "--out",
        "out",
        "--run-mode",
        "pipeline",
    ])
    .expect("cli parse pipeline");
    let Commands::Run(args) = cli.command;
    assert_eq!(args.run_mode, RunMode::Pipeline);
}

#[test]
fn prefix_detection_prefixed_and_nonprefixed() {
    let tmp = tempdir().unwrap();
    let nonpref = tmp.path().join("nonpref");
    let pref = tmp.path().join("pref");
    create_dir_all(&nonpref).unwrap();
    create_dir_all(&pref).unwrap();

    write(nonpref.join("matrix.mtx"), "").unwrap();
    write(nonpref.join("barcodes.tsv"), "").unwrap();
    write(nonpref.join("features.tsv"), "").unwrap();

    write(pref.join("XYZ_matrix.mtx"), "").unwrap();
    write(pref.join("XYZ_barcodes.tsv"), "").unwrap();
    write(pref.join("XYZ_features.tsv"), "").unwrap();

    assert_eq!(detect_prefix(&nonpref).unwrap(), None);
    assert_eq!(detect_prefix(&pref).unwrap(), Some("XYZ".to_string()));
}

#[test]
fn shared_cache_filename_resolution() {
    let base = Path::new("/tmp/input");
    assert_eq!(
        resolve_cache_path(base, None),
        base.join("kira-organelle.bin")
    );
    assert_eq!(
        resolve_cache_path(base, Some("XYZ")),
        base.join("XYZ.kira-organelle.bin")
    );
}

#[test]
fn shared_cache_read_valid_and_traversal() {
    let tmp = tempdir().unwrap();
    let cache_path = tmp.path().join("kira-organelle.bin");

    write_valid_cache(
        &cache_path,
        &["RPLP0", "EEF2", "ATF4"],
        &["C1", "C2"],
        &[0, 2, 3],
        &[0, 2, 1],
        &[5, 1, 7],
        false,
    );

    let cache = read_shared_cache(&cache_path).unwrap();
    assert_eq!(cache.n_genes, 3);
    assert_eq!(cache.n_cells, 2);
    assert_eq!(cache.nnz, 3);
    assert_eq!(cache.genes, vec!["RPLP0", "EEF2", "ATF4"]);
    assert_eq!(cache.barcodes, vec!["C1", "C2"]);
    assert_eq!(cache.col_ptr, vec![0, 2, 3]);
    assert_eq!(cache.row_idx, vec![0, 2, 1]);
    assert_eq!(cache.values_u32, vec![5, 1, 7]);
}

#[test]
fn shared_cache_header_crc_validation() {
    let tmp = tempdir().unwrap();
    let cache_ok = tmp.path().join("ok.kira-organelle.bin");
    let cache_bad = tmp.path().join("bad.kira-organelle.bin");

    write_valid_cache(&cache_ok, &["G1"], &["C1"], &[0, 1], &[0], &[9], false);
    write_valid_cache(&cache_bad, &["G1"], &["C1"], &[0, 1], &[0], &[9], true);

    assert!(read_shared_cache(&cache_ok).is_ok());
    let err = read_shared_cache(&cache_bad).unwrap_err();
    assert!(err.to_string().contains("CRC"));
}

#[test]
fn pipeline_cache_present_uses_cache_path() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("out");
    create_dir_all(&input).unwrap();

    let cache_path = input.join("kira-organelle.bin");
    write_valid_cache(
        &cache_path,
        &["RPLP0", "EEF2", "ATF4"],
        &["C1", "C2"],
        &[0, 2, 3],
        &[0, 2, 1],
        &[5, 1, 7],
        false,
    );

    let args = RunArgs {
        input,
        out,
        mode: Mode::Cell,
        metadata: None,
        run_mode: RunMode::Pipeline,
        disable_translation_extension: false,
    };

    let (bundle, stats) = run_stage1(&args).unwrap();
    match bundle.format {
        InputFormat::SharedCache { cache_path: p } => assert_eq!(p, cache_path),
        _ => panic!("expected shared cache input format"),
    }
    assert!(bundle.shared_cache.is_some());
    assert_eq!(stats.n_cells, 2);
}

#[test]
fn pipeline_cache_missing_falls_back_to_mtx() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("out");
    create_dir_all(&input).unwrap();

    write(
        input.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n3 2 2\n1 1 5\n2 2 3\n",
    )
    .unwrap();
    write(input.join("barcodes.tsv"), "C1\nC2\n").unwrap();
    write(
        input.join("features.tsv"),
        "G1\tRPLP0\nG2\tEEF2\nG3\tATF4\n",
    )
    .unwrap();

    let args = RunArgs {
        input,
        out,
        mode: Mode::Cell,
        metadata: None,
        run_mode: RunMode::Pipeline,
        disable_translation_extension: false,
    };

    let (bundle, _) = run_stage1(&args).unwrap();
    assert!(matches!(bundle.format, InputFormat::TenXDir { .. }));
    assert!(bundle.shared_cache.is_none());
}

#[test]
fn pipeline_cache_invalid_is_hard_error() {
    let tmp = tempdir().unwrap();
    let input = tmp.path().join("input");
    let out = tmp.path().join("out");
    create_dir_all(&input).unwrap();

    let cache_path = input.join("kira-organelle.bin");
    write_valid_cache(&cache_path, &["G1"], &["C1"], &[0, 1], &[0], &[9], true);

    let args = RunArgs {
        input,
        out,
        mode: Mode::Cell,
        metadata: None,
        run_mode: RunMode::Pipeline,
        disable_translation_extension: false,
    };

    let err = run_stage1(&args).unwrap_err();
    assert!(err.to_string().contains("CRC"));
}

fn write_valid_cache(
    path: &Path,
    genes: &[&str],
    barcodes: &[&str],
    col_ptr: &[u64],
    row_idx: &[u32],
    values: &[u32],
    tamper_header_crc: bool,
) {
    let genes_table = encode_string_table(genes);
    let barcodes_table = encode_string_table(barcodes);

    let mut offset = 256usize;
    let genes_offset = align64(offset);
    offset = genes_offset + genes_table.len();

    let barcodes_offset = align64(offset);
    offset = barcodes_offset + barcodes_table.len();

    let col_ptr_offset = align64(offset);
    let col_ptr_bytes = col_ptr.len() * 8;
    offset = col_ptr_offset + col_ptr_bytes;

    let row_idx_offset = align64(offset);
    let row_idx_bytes = row_idx.len() * 4;
    offset = row_idx_offset + row_idx_bytes;

    let values_offset = align64(offset);
    let values_bytes = values.len() * 4;
    offset = values_offset + values_bytes;

    let file_bytes = offset;
    let mut file = vec![0u8; file_bytes];

    file[genes_offset..genes_offset + genes_table.len()].copy_from_slice(&genes_table);
    file[barcodes_offset..barcodes_offset + barcodes_table.len()].copy_from_slice(&barcodes_table);

    for (i, v) in col_ptr.iter().enumerate() {
        let start = col_ptr_offset + i * 8;
        file[start..start + 8].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in row_idx.iter().enumerate() {
        let start = row_idx_offset + i * 4;
        file[start..start + 4].copy_from_slice(&v.to_le_bytes());
    }
    for (i, v) in values.iter().enumerate() {
        let start = values_offset + i * 4;
        file[start..start + 4].copy_from_slice(&v.to_le_bytes());
    }

    file[0..4].copy_from_slice(b"KORG");
    file[4..6].copy_from_slice(&1u16.to_le_bytes());
    file[6..8].copy_from_slice(&0u16.to_le_bytes());
    file[8..12].copy_from_slice(&0x1234_5678u32.to_le_bytes());
    file[12..16].copy_from_slice(&256u32.to_le_bytes());

    put_u64(&mut file, 16, genes.len() as u64);
    put_u64(&mut file, 24, barcodes.len() as u64);
    put_u64(&mut file, 32, row_idx.len() as u64);

    put_u64(&mut file, 40, genes_offset as u64);
    put_u64(&mut file, 48, genes_table.len() as u64);
    put_u64(&mut file, 56, barcodes_offset as u64);
    put_u64(&mut file, 64, barcodes_table.len() as u64);

    put_u64(&mut file, 72, col_ptr_offset as u64);
    put_u64(&mut file, 80, row_idx_offset as u64);
    put_u64(&mut file, 88, values_offset as u64);

    put_u64(&mut file, 96, 0);
    put_u64(&mut file, 104, 0);
    put_u64(&mut file, 112, file_bytes as u64);
    put_u64(&mut file, 128, 0);

    let mut header = file[0..256].to_vec();
    header[120..128].fill(0);
    let mut crc = crc64_ecma(&header);
    if tamper_header_crc {
        crc = crc.wrapping_add(1);
    }
    put_u64(&mut file, 120, crc);

    write(path, file).unwrap();
}

fn encode_string_table(strings: &[&str]) -> Vec<u8> {
    let mut blob = Vec::new();
    let mut offsets = Vec::with_capacity(strings.len() + 1);
    offsets.push(0u32);
    for s in strings {
        blob.extend_from_slice(s.as_bytes());
        offsets.push(blob.len() as u32);
    }

    let mut out = Vec::new();
    out.extend_from_slice(&(strings.len() as u32).to_le_bytes());
    for off in offsets {
        out.extend_from_slice(&off.to_le_bytes());
    }
    out.extend_from_slice(&blob);
    out
}

fn align64(v: usize) -> usize {
    (v + 63) & !63
}

fn put_u64(buf: &mut [u8], offset: usize, value: u64) {
    buf[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn crc64_ecma(bytes: &[u8]) -> u64 {
    let mut crc = 0u64;
    for &b in bytes {
        crc ^= (b as u64) << 56;
        for _ in 0..8 {
            if (crc & 0x8000_0000_0000_0000) != 0 {
                crc = (crc << 1) ^ 0x42F0_E1EB_A9EA_3693;
            } else {
                crc <<= 1;
            }
        }
    }
    crc
}
