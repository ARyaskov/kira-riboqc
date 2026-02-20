# kira-riboqc Metrics Specification

This document defines canonical metrics produced by `kira-riboqc`, including formulas, constants, thresholds, directionality, and output contracts.

Scope:
- core axis/score/classification pipeline (`riboqc.tsv`, `summary.json`, `pipeline_step.json`)
- optional translation-extension (`translation_regime_metrics.tsv`)

## Canonical Conventions

1. Determinism
- No stochastic steps.
- Fixed constants, fixed rule order, fixed sort order (`barcode` / `cell_id` lexical).

2. Sample semantics
- All metrics are computed per cell/barcode.

3. Naming
- Stable snake_case metric ids and stable string regime labels.

4. Value domains
- Enrichment-like panel signals: `[0,1]` via `clamp01`.
- TL/ST/RQC/TSS/LFTI/RAS: `[0,1]` via `clamp01`.
- `tpc`: `[-1,1]` via `clamp_to_minus1_plus1`.
- `confidence`: `[0,1]`.

5. Directionality
- Higher means stronger intensity (or risk) for:
  - `translation_load`, `stress_translation_index`, `ribosomal_specialization`, `elongation_pressure`
  - `tss`, `lfti`, `ras`
  - translation-extension commitment/selectivity metrics
- Higher is better only for:
  - `confidence`

## Notation

Let:
- `count(g,c)` be raw UMI count for gene `g` in cell `c`.
- `libsize(c) = sum_g count(g,c)`.
- `cpm(g,c) = 1e6 * count(g,c) / max(libsize(c), 1)`.
- `x(g,c) = ln(1 + cpm(g,c))`.
- `mean(S,c)` arithmetic mean of `x(g,c)` over genes in set `S` observed in gene index.
- `clamp01(v) = min(max(v,0),1)`.
- `clamp11(v) = min(max(v,-1),1)`.

Panel enrichment transform used in Stage 2 and translation extension:
- `OFFSET = 0.5`
- `SCALE = 1.5`
- `enrich(set_mean, bg_mean) = clamp01((set_mean - bg_mean + OFFSET) / SCALE)`

Background gene set:
- top `BG_TOP_K = 256` by detected-cell count
- exclude genes with symbol prefix `RPL*`, `RPS*`, `MT-*`, `MT.*`

## Core Components and Axes (Stage 2)

Gene sets:
- `RRNA_PROXY = [NOP56, NOP58, FBL, DKC1]`
- `MACHINERY = [EEF1A1, EEF2, EEF1B2, EIF4A1, EIF4G1]`
- `STRESS = [ATF4, DDIT3, XBP1, HSPA5, DNAJB9]`
- `METAB = [SLC2A1, LDHA, PFKP, PDK1]`
- `LUXURY = [PCNA, HMGB2, TYMS]`
- `RQC_MACH = [ZNF598, ASCC3, PELO, HBS1L, NEMF, LTN1, TCF25]`
- `CHAPERONE = [HSPA1A, HSPA8, HSP90AA1, DNAJB1]`
- `UBIQUITIN = [UBB, UBC, SQSTM1]`

Per-cell component signals (all in `[0,1]` unless absent):
- `rp = enrich(mean(RIBOSOMAL), mean(BG))`
- `mach = enrich(mean(MACHINERY), mean(BG))`
- `rrna_proxy = enrich(mean(RRNA_PROXY), mean(BG))`
- `stress = enrich(mean(STRESS), mean(BG))`
- `metab = enrich(mean(METAB), mean(BG))`
- `luxury = enrich(mean(LUXURY), mean(BG))`, optional (`None` if unresolved)
- `rqc_mach = enrich(mean(RQC_MACH), mean(BG))`
- `chaperone = enrich(mean(CHAPERONE), mean(BG))`
- `ubiquitin = enrich(mean(UBIQUITIN), mean(BG))`

Axis definitions:
- `translation_load (tl) = clamp01(0.5*rp + 0.3*mach + 0.2*rrna_proxy)`
- `stress_translation_index (st)`:
  - if `luxury` exists:
    - `lux_supp = clamp01(1 - luxury)`
    - `st = clamp01(0.45*stress + 0.35*metab + 0.20*lux_supp)`
  - else:
    - `st = clamp01(0.55*stress + 0.45*metab)`
    - set `st_low_confidence = true`
- `ribosomal_specialization (rqc) = clamp01(0.7*rqc_mach + 0.3*stress)`
- `elongation_pressure (rqc_pressure) = clamp01(rqc_mach - rp)`
- `deg_cap = clamp01(0.6*chaperone + 0.4*ubiquitin)`
- `tpc = clamp11(deg_cap - tl)`

## Composite Scores (Stage 3)

Fragility proxy:
- `fragility = clamp01(0.5 - tpc/2)`

Score formulas:
- `tss = clamp01(0.45*tl + 0.35*rqc + 0.20*st)`
- `lfti = clamp01(0.40*rqc + 0.35*st + 0.25*fragility)`
- `ras = clamp01(0.50*tl + 0.30*rqc + 0.20*fragility)`

Red-flag predicate:
- `ras_red_flag = (tl > 0.70) && (rqc > 0.70) && (tpc < -0.20)`

## Regime Classification (Stage 4)

QC constants:
- `MIN_LIBSIZE = 500`
- `MIN_DETECTED_GENES = 200`

Threshold constants:
- `TL_LOW = 0.25`, `TL_MID_LOW = 0.35`, `TL_MID_HIGH = 0.70`, `TL_HIGH = 0.75`
- `ST_LOW = 0.35`, `ST_HIGH = 0.65`
- `RQC_LOW = 0.35`, `RQC_MID = 0.45`, `RQC_HIGH = 0.70`
- `TPC_NEG_SOFT = -0.05`, `TPC_NEG = -0.15`

Rule order is deterministic and priority-based:
1. `tl < TL_LOW && st < ST_LOW && rqc < RQC_LOW` -> `TranslationSuppressed`
2. `TL_MID_LOW <= tl <= TL_MID_HIGH && rqc < RQC_MID && tpc >= TPC_NEG_SOFT` -> `EfficientTranslation`
3. `st >= ST_HIGH && rqc >= RQC_MID` -> `SelectiveSurvivalTranslation`
4. `tl >= TL_HIGH && rqc >= RQC_MID && tpc < TPC_NEG` -> `OverloadedTranslation`
5. `rqc >= RQC_HIGH && tl >= TL_MID_LOW && tpc < 0` -> `RQC_Dependent`
6. else -> `Unclassified`

Run-level QC flags:
- `low_counts_cell = (libsize < MIN_LIBSIZE)`
- `few_detected_genes = (detected_genes < MIN_DETECTED_GENES)`

## Public Regime Mapping and Confidence

Internal regime -> output `regime` (in `riboqc.tsv`):
- `EfficientTranslation` -> `GrowthDrivenTranslation` if `translation_load >= 0.60`, else `HomeostaticTranslation`
- `SelectiveSurvivalTranslation` -> `StressAdaptiveTranslation`
- `OverloadedTranslation` -> `TranslationalOverdrive`
- `RQC_Dependent` -> `StressAdaptiveTranslation`
- `TranslationSuppressed` -> `TranslationalCollapse`
- `Unclassified` -> `Unclassified`

Flags per row (`flags`, comma-separated, deduplicated):
- `LOW_COUNTS_CELL`
- `FEW_DETECTED_GENES`
- `LOW_CONFIDENCE` (from `st_low_confidence` and/or translation extension low coverage)
- `LOW_RIBO_SIGNAL` (`ribosome_density < 0.2`)
- `RAS_RED_FLAG`

Confidence score:
- initialize `conf = 1.0`
- `-0.30` if `LOW_COUNTS_CELL`
- `-0.25` if `FEW_DETECTED_GENES`
- `-0.20` if `LOW_CONFIDENCE`
- `-0.15` if `LOW_RIBO_SIGNAL`
- final `confidence = clamp01(conf)`

## Translation Extension Metrics (Optional)

Output file: `translation_regime_metrics.tsv`.

Coverage constants:
- `LOW_COVERAGE_MIN_RATIO = 0.35`
- panel low-confidence when `found == 0` or `found/total < 0.35` for any tracked set
- if any set is low coverage, every row gets `low_confidence = true` (propagates to core `LOW_CONFIDENCE` flag)

Per-cell panel enrichments (`[0,1]`):
- `ribosome = enrich(mean(RIBOSOME_CORE), mean(BG))`
- `housekeeping = enrich(mean(HOUSEKEEPING), mean(BG))`
- `global = enrich(mean(GLOBAL_TRANSLATION), mean(BG))`
- `stress = enrich(mean(STRESS_RESPONSIVE), mean(BG))`
- `isr = enrich(mean(ISR_LIKE), mean(BG))`
- `codon = enrich(mean(CODON_SUBOPTIMAL), mean(BG))`

Metrics:
- `ribosome_loading_heterogeneity = clamp01(0.6*abs(ribosome-housekeeping) + 0.4*abs(global-housekeeping))`
- `translation_selectivity_index = clamp01((0.6*stress + 0.4*isr) / (0.5*global + 0.5*ribosome + 0.2))`
- `isr_like_signature_score = clamp01(0.7*isr + 0.3*stress - 0.15*global + 0.15)`
- `codon_bias_proxy = clamp01(0.7*codon + 0.3*stress - 0.25*housekeeping + 0.1)`
- `translation_commitment_score = clamp01(0.40*heterogeneity + 0.35*selectivity + 0.25*isr_score)`

Threshold:
- `HIGH_SELECTIVE_THRESHOLD = 0.65` for run-level high selective fraction.

Classification (priority order):
1. `commitment >= 0.78 && selectivity >= 0.62 && isr >= 0.58 && codon_bias >= 0.50` -> `TranslationFixation`
2. `isr >= 0.68 && selectivity >= 0.45` -> `ISR_DominatedTranslation`
3. `selectivity >= 0.58 && heterogeneity >= 0.42` -> `SelectiveTranslationStress`
4. else -> `BalancedTranslation`

Run-level translation summary:
- `regime_fractions`: fractions for `BalancedTranslation`, `SelectiveTranslationStress`, `ISR_DominatedTranslation`, `TranslationFixation`
- `mean_translation_commitment = mean(translation_commitment_score)`
- `high_selective_translation_fraction = frac(translation_selectivity_index >= 0.65)`

## Output Contracts

## `riboqc.tsv`

Header:
- `barcode`, `sample`, `condition`, `species`, `libsize`, `nnz`, `expressed_genes`,
  `translation_load`, `ribosome_density`, `elongation_pressure`, `initiation_bias`,
  `ribosomal_specialization`, `stress_translation_index`, `regime`, `flags`, `confidence`

Formatting:
- numeric metric columns emitted with 6 decimals.
- row order sorted lexicographically by `barcode`.

## `translation_regime_metrics.tsv` (optional)

Header:
- `cell_id`, `ribosome_loading_heterogeneity`, `translation_selectivity_index`,
  `isr_like_signature_score`, `codon_bias_proxy`, `translation_commitment_score`,
  `translation_regime`

Formatting:
- numeric metric columns emitted with 6 decimals.
- row order sorted lexicographically by `cell_id`.

## `summary.json`

Object fields:
- `tool`:
  - `name`: `"kira-riboqc"`
  - `version`: cargo package version
  - `simd`: compiled SIMD backend id
- `input`:
  - `n_cells`
  - `species`
- `distributions`:
  - `translation_load`: `{median,p90,p99}`
  - `ribosome_density`: `{median,p90,p99}`
  - `stress_translation_index`: `{median,p90,p99}`
- `regimes`:
  - `counts`: fixed keys `HomeostaticTranslation`, `GrowthDrivenTranslation`, `StressAdaptiveTranslation`, `TranslationalOverdrive`, `TranslationalCollapse`, `Unclassified`
  - `fractions`: same keys, rounded to 6 decimals
- `qc`:
  - `low_confidence_fraction`
  - `low_ribo_signal_fraction`
- optional `translation`:
  - `regime_fractions`
  - `mean_translation_commitment`
  - `high_selective_translation_fraction`

Percentile rule:
- nearest-rank with `rank = ceil(p*N)`, index `rank-1`, finite-safe sorting (non-finite treated as `0`).

## `pipeline_step.json` (pipeline run mode)

Contract:
- `tool`: `{name:"kira-riboqc", stage:"translation", version}`
- `artifacts`:
  - `summary: "summary.json"`
  - `primary_metrics: "riboqc.tsv"`
  - `panels: "panels_report.tsv"`
  - optional `translation_metrics: "translation_regime_metrics.tsv"`
- `cell_metrics`:
  - `file: "riboqc.tsv"`
  - `id_column: "barcode"`
  - `regime_column: "regime"`
  - `confidence_column: "confidence"`
  - `flag_column: "flags"`
  - optional `translation_mode_column: "translation_regime"`
- `regimes`: canonical public regime list
