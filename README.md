# kira-riboqc

Deterministic ribosome and translation-state quality control for single-cell RNA-seq.

## Build requirements

- Rust >= 1.95

## Install

Install from crates.io:

```bash
cargo install kira-riboqc
```

## Usage examples

Standalone run (cell mode, default):

```bash
kira-riboqc run \
  --input ./data/pbmc3k \
  --out ./out/pbmc3k \
  --mode cell
```

Standalone run (sample mode):

```bash
kira-riboqc run \
  --input ./data/inf \
  --out ./out/inf \
  --mode sample
```

Standalone run with metadata:

```bash
kira-riboqc run \
  --input ./data/inf \
  --metadata ./data/inf/metadata.tsv \
  --out ./out/inf
```

Pipeline run (shared cache lookup + pipeline artifacts):

```bash
kira-riboqc run \
  --input ./data/inf \
  --out ./out/inf \
  --mode sample \
  --run-mode pipeline
```

Pipeline run without translation extension:

```bash
kira-riboqc run \
  --input ./data/inf \
  --out ./out/inf \
  --run-mode pipeline \
  --disable-translation-extension
```

## Modes

- `--mode cell` (default): per-cell QC and per-cell artifacts.
- `--mode sample`: aggregate-centric sample run.
- `--run-mode standalone` (default): existing standalone behavior and outputs.
- `--run-mode pipeline`: pipeline contract mode for `kira-organelle`.

## Pipeline cache lookup

In pipeline mode, `kira-riboqc` first searches the input directory for shared cache as specified in [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md):

- no prefix: `kira-organelle.bin`
- prefixed dataset: `<PREFIX>.kira-organelle.bin`

Prefix is detected non-recursively from names like `<PREFIX>_matrix.mtx(.gz)`, `<PREFIX>_features.tsv(.gz)`, `<PREFIX>_barcodes.tsv(.gz)`.

Behavior:

- cache exists and valid: use shared cache path.
- cache missing: warn once and fall back to MTX input.
- cache exists but invalid: hard error (no silent fallback).

## Pipeline output contract

In pipeline mode, outputs are written to:

- `--out <DIR>` -> `<DIR>/kira-riboqc/`

Required artifacts:

- `riboqc.tsv` (per-cell contract table)
- `summary.json` (run-level aggregates)
- `panels_report.tsv` (panel audit)
- `pipeline_step.json` (ingestion manifest for `kira-organelle`)
- `translation_regime_metrics.tsv` (additive translation extension artifact, when enabled)

Standalone mode also writes `report.txt`.

## Shared cache specification

- Cache format specification: [kira-shared-sc-cache/CACHE_FILE.md](https://github.com/ARyaskov/kira-shared-sc-cache/blob/main/CACHE_FILE.md)
- Reader validates header/magic/version/endian/header-size/file-bytes, header CRC64-ECMA, section bounds, string tables, and CSC invariants.

## SIMD note

- SIMD backend is selected at compile time.
- The chosen backend is logged at startup.
- Scalar fallback is always available.
