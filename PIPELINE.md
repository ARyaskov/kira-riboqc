Stage 1 - Input loading and gene index
Purpose: detect input source and build deterministic matrix/barcode/feature structures.
Inputs: input directory, run mode.
Outputs: InputBundle + Stage1Stats.
Determinism guarantees: stable file resolution, stable gene index ordering, strict cache validation in pipeline mode.

Stage 2 - Normalization and translation axes
Purpose: compute per-cell normalization and translation axis components.
Inputs: InputBundle.
Outputs: Stage2Output.
Determinism guarantees: fixed constants, stable sorting, no parallel randomness.

Stage 3 - Composite translation stress scores
Purpose: compute TSS/LFTI/RAS and explanation strings.
Inputs: Stage2Output.
Outputs: Stage3Output.
Determinism guarantees: fixed formulae and stable explanation ordering.

Stage 4 - Rule-based regime classification
Purpose: assign one regime per cell using threshold rules.
Inputs: Stage2Output + Stage3Output.
Outputs: Stage4Output.
Determinism guarantees: explicit rule order, direct comparisons, fixed thresholds.

Stage 5 - Translation-regime extension (expression-only proxies)
Purpose: compute derived translational stress/adaptation proxy metrics and classify translation-mode regimes.
Inputs: InputBundle + Stage2Output.
Outputs: StageTranslationRegimeOutput (optional, controlled by CLI flag).
Determinism guarantees: fixed gene panels from `resources/ribosome/translation`, fixed thresholds, no stochastic steps.
Notes: this stage is a proxy-only extension; it does not measure ribosome occupancy or direct translation rates.

Stage 6 - Reporting and pipeline contract artifacts
Purpose: emit machine-readable outputs for `kira-organelle` ingestion.
Inputs: InputBundle + Stage1Stats + Stage2Output + Stage3Output + Stage4Output + optional StageTranslationRegimeOutput.
Outputs:
- pipeline mode: `<out>/kira-riboqc/riboqc.tsv`, `summary.json`, `panels_report.tsv`, `pipeline_step.json`
- extension enabled: additional `translation_regime_metrics.tsv` and `summary.json.translation`
- standalone mode: same base artifacts in `<out>/` plus `report.txt`
Determinism guarantees: fixed TSV schema/order, sorted aggregate calculations, stable JSON keys.

Stage 7 - Performance optimizations
Purpose: membership-vector lookup and SIMD scaffolding.
Inputs: stage internals.
Outputs: same analytical values with lower overhead.
Determinism guarantees: scalar-equivalent math and compile-time backend selection.
