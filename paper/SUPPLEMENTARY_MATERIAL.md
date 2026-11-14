# Supplementary Material

All tables referenced here are machine-readable CSV files in `research/results/`; the commands that produce them are listed in `research_docs/EXPERIMENT_LOG.md`.

## S1. Data processing
* Build script: `research/preprocessing/build_dataset.py` (verifies raw SHA-256 before use). Counts: raw 1,402,829; range-only QC 512,384; strict QC 428,264 (drops: calculated or flagged DIC 35,819; DIC without secondary QC 67,839; salinity or oxygen flag/QC 49,351, overlapping). Statistics in `research/data/build_stats.json`.
* The repository's original processed dataset is reproduced exactly from the raw file by range-only QC.

## S2. Estimator validation
* `e01_synthetic_recovery.csv`: bias and coverage for 48 settings (Fig. S1).
* `e01b_coverage_followup.csv`: 200-replicate comparison of bootstrap and CR1 intervals in the worst settings.

## S3. Main results tables
* `e02_layer_trends_primary.csv`: per basin × layer — n, cells, cruises, median within-cell span, detection limit, β_DIC, β_Cant (with 95 % CIs and CR1 standard errors), R, δ, and direct TRACE dCant/dt.
* `summary_h1.csv`: H1 verdicts.
* `summary_robustness.csv`: δ, CI, sign stability, range over variants, detection limit.
* `summary_h4.csv`: strict vs relaxed QC deep-layer comparison.

## S4. Sensitivity analyses
`e02_layer_trends_<variant>.csv` for: S1 relaxed QC; S2 cell size 2°/6°/8°; S3 span 5/15 yr; S4 linear conditioning; S5 no AOU; S6 no θ, S; S7 post-1990; S9 seasonal harmonics; S10 phosphate instead of AOU; B1 no conditioning. Leave-one-cruise-out: `e04_loco.csv`.

## S5. eMLR(C*) comparison
`e05_emlr_comparison.csv`: 1989–2012 window trends, TRACE 1994–2007 rate, eMLR(C*) rate (standard case, σ, ensemble range), mechanical verdict.

## S6. Static empirical models (H3)
* `e03_temporal_{relaxed_none-year, strict_none-year-canth}.csv`: train/validation/external RMSE, bias and MAE per basin × model × covariate.
* `e03_blocked_summary_*.csv`, `e03_blocked_folds_*.csv`: spatial- and cruise-blocked CV.
* `e03_paired_*.csv`: paired cruise-cluster bootstrap comparisons (±4 µmol kg⁻¹ equivalence).
* `summary_h3_*.csv`, `summary_h3_rankings_*.csv`: bias-removal fractions and model rankings by validation design.

## S7. Audit of the original repository analysis
Documented in `research_docs/REPOSITORY_INTELLIGENCE.md` §8, including:
* the incorrect cluster bootstrap (`np.isin` drops duplicate cluster draws);
* the anticonservative per-observation TOST;
* leakage of the frozen symbolic-regression candidate into blocked CV;
* the string-based equation-family classifier;
* hard-coded results;
* unused GLODAP QC flags.

## S8. Software
Python 3.12.16; numpy, pandas, scipy, scikit-learn, xarray; TRACE-Python (package version 1.1.1, release tag v1.1.2; zip SHA-256 01cf923e…0b51) with PyCO2SYS 2.0.0b8.

---
# Round-2 supplementary material (2026-11-14)

## S-R2.1 Differential estimator: derivation and solver
For one basin, let T be the target-layer samples and R the reference-layer samples (γⁿ > 27.8), restricted to cruises with samples in both. The model is y = Dₐα + D_cγ + Xφ + Δβ·t·1_T + ε, where:
- Dₐ: dummies for cell × {T, R};
- D_c: cruise dummies;
- X: layer-specific quadratic conditioning terms.

**Solution.**
1. Sweep out Dₐ by within-(cell×layer) demeaning: Z̃ = Z − P_a Z.
2. With weights w (bootstrap multiplicities), the normal equations for (φ, Δβ, γ) partition into
 - S_ZZ = Z̃ᵀWZ̃;
 - S_DZ = per-cruise weighted sums of Z̃;
 - S_DD = diag(w_c) − Nᵀ diag(1/w_a) N, where N is the (cell × cruise) weighted count matrix.
3. Eliminate γ with the pseudo-inverse of S_DD (rank-deficient by one): S = S_ZZ − S_DZᵀ S_DD⁺ S_DZ. The coefficients follow from the X-block of S.

**Checks and cost.**
- Agrees with alternating projections to 1.8·10⁻¹⁴ (`tests/test_research_estimators.py`).
- Cost is O(n·k) per bootstrap replicate plus a G×G pseudo-inverse (G ≤ 300 cruises).

**Identification.** Δβ comes from how the within-cruise target-minus-reference contrast varies across cruises of different dates, net of cell effects. Any cruise-level additive offset, constant or trending in time, is absorbed by γ.

## S-R2.2 Synthetic drift test (E06s)
**Setup.** Three groups (Pacific 26.0–26.5, Atlantic 26.5–27.0, Indian 26.5–27.0). Each synthetic series is built from:
- the fitted static two-way-FE component of the observed data;
- true trends of 0.8 (target) and 0.1 (reference), so the true Δβ = 0.7;
- cruise offsets d·(t_cruise − 2000) + N(0, 3²);
- permuted residuals.

60 replicates; B = 150.

**Results.**
- Differential estimator: bias ≤ 0.002 µmol kg⁻¹ yr⁻¹; coverage 93–100 % for d ∈ {0, 0.05, 0.1}.
- Single-layer estimator: absorbs d one-for-one.

File: `research/results/e06_synthetic_drift.csv`.

## S-R2.3 TRACE structural variants (E07) and boundary condition (E12)
TRACE inputs (lon, lat, depth, S, T, decimal year) are unchanged across variants.

| Variant | Change | Upper-layer Δδ range (µmol kg⁻¹ yr⁻¹) |
|---|---|---|
| Δ/Γ 1.0 / 1.6 | TTD width (ages re-estimated by TRACE's network) | ≤ 0.015 |
| θ input | potential instead of in-situ T | ≤ 0.004 |
| mean age ×0.8 / ×1.2 | scale factors × constant | −0.04 / +0.03 |
| diseq ×0.9 / ×1.1 | scales atmospheric CO₂ anomaly | +0.03…+0.07 / −0.03…−0.07 |
| unadjusted atm. trajectory (E12) | CO2Trajectories.txt instead of CO2TrajectoriesAdjusted.txt | −0.05…−0.12 |

TRACE's adjusted trajectory has anomaly/atmospheric-anomaly ratios of 0.912 (1960), 0.906 (1990), 0.903 (2000), 0.898 (2010) and 0.893 (2020). Its 1972–2021 slope is 1.57 vs 1.78 ppm yr⁻¹ for the atmosphere (ratio 0.883).

## S-R2.4 C* and alkalinity control (E08)
C* = DIC − 117·PO₄ − ½(TAlk + 16·PO₄), on samples with PO₄ and TAlk WOCE flag 2 and secondary-QC 1. Full table: `research/results/e08_cstar_controls.csv`.

TAlk trends at fixed θ, S with 95 % intervals excluding 0 are all negative:
- Pacific 27.5–27.8, 27.8–28.0, 28.0–28.1, 28.1–28.2: −0.057, −0.057, −0.048, −0.074;
- Indian 27.5–27.8: −0.082 µmol kg⁻¹ yr⁻¹.

## S-R2.5 Pooled statistics (E09, E12) and eMLR products (E10)
Files:
- `e09_pooled.csv`, `e09_pooled_atm.csv`: pooled δ and R, global cruise bootstrap B = 500.
- `summary_round2_emlr_pooled.csv`: sample-weighted product/TRACE ratios, with sensitivity-case minima and maxima.
- `e09_piecewise_time.csv`: slopes before and after 2000.

## S-R2.6 Reproduction
Round-1 results were reproduced bit-for-bit on a second machine (≤ 2·10⁻¹⁴); see `EXPERIMENT_LOG.md` R0.

Full Round-2 pipeline after `build_dataset` and `e02a_run_trace` (about 45 min on a laptop):

```
python -m research.experiments.e06_differential_trends --boot 300
python -m research.experiments.e06_differential_trends --synthetic --reps 60 --boot 150
python -m research.experiments.e07_trace_sensitivity --boot 300
python -m research.experiments.e08_cstar_controls --boot 300
python -m research.experiments.e09_pooled_and_timeshape --boot 500
python -m research.experiments.e10_emlr_products
python -m research.experiments.e11_scale_factor --boot 500
python -m research.experiments.e12_trace_unadjusted --fit --boot 300
python -m research.experiments.e09_pooled_and_timeshape --boot 500 --cant atm
python -m research.experiments.e11_scale_factor --boot 500 --cant atm
python -m research.evaluation.summarize_round2
python -m research.visualization.figures_round2
```

# Round-3 supplementary material (2026-11-14)

## S-R3.1 Earth-system-model sampling (E13)

**Data.** MRI-ESM2-0, CMIP6 r1i2p1f1, regular 1° grid (`gr`), version v20200303. Historical for 1972–2014; ssp585 for 2015–2021 (scenario differences before 2021 are negligible). Variables:
- annual `dissic`, `dissicnat`, `o2`;
- monthly `thetao`, `so`, averaged to calendar years.

**Sampling.** Each sample takes model values for its calendar year:
- horizontally, the nearest model column (great-circle distance) wet at both bracketing levels;
- vertically, linear interpolation between levels; depths beyond the deepest wet level (6150 m) are clamped.

**Conversions** (all from model fields, via TEOS-10/GSW):
- concentrations from mol m⁻³ to µmol kg⁻¹ with in-situ density;
- AOU = O₂sol(S, θ) − O₂ (Garcia & Gordon 1992 fit);
- in-situ temperature for TRACE from θ.

The analysis uses the observed basin, layer, cell, cruise and date of every sample, so the sampling design is identical to the observational analysis.

**Fidelity.**
- Model vs observed correlations at the samples: θ 0.975, S 0.86, AOU 0.89, DIC 0.90.
- Mean biases: θ +0.8 °C, AOU −14 µmol kg⁻¹, DIC −26 µmol kg⁻¹.
- TRACE returns NaN for model T, S at 1.8 % of samples (warm, salty upper-thermocline waters); these are excluded.

## S-R3.2 Perfect-model results (E13, E13d)

Pooled upper-ocean ratios (global shared cruise bootstrap, 500 replicates):

| Comparison | All | Atlantic | Indian | Pacific | Southern Ocean |
|---|---|---|---|---|---|
| model DIC / C_ant^true | 0.97 | 1.00 | 1.00 | 0.95 | 0.99 |
| C_ant^true / TRACE (default) | 1.12 | 1.04 | 1.01 | 1.18 | 1.06 |
| C_ant^true / TRACE (unadjusted atm.) | 1.02 | 0.94 | 0.92 | 1.08 | 0.97 |

Differential scale factors m (all 20 layers): DIC vs truth −0.01; truth vs TRACE 0.125 (default) and 0.024 (unadjusted). Full tables are in `research/results/e13_*.csv`.

## S-R3.3 Alternative reference layers (E14)

| Reference | m (all basins) | Atlantic | Indian | Pacific | Southern Ocean |
|---|---|---|---|---|---|
| γⁿ > 27.8 (primary) | 0.18 [0.11, 0.23] | 0.20 | 0.19 | 0.16 | 0.18 |
| γⁿ 27.8–28.0 | 0.16 [0.11, 0.21] | 0.18 | 0.15 | 0.15 | 0.17 |
| γⁿ > 28.0 | 0.19 [0.12, 0.25] | 0.19 | 0.21 | 0.18 | 0.19 |

## S-R4.1 Within-cruise depth-dependent artefacts (E16, E16b)

**Near reference (E16).** Targets are the layers with γⁿ ≤ 27.0; the reference is 27.0 < γⁿ ≤ 27.8.

| Reference | m (all basins) | Atlantic | Indian | Pacific | Southern Ocean |
|---|---|---|---|---|---|
| γⁿ > 27.8 (primary, γⁿ ≤ 27.0 targets) | 0.18 [0.11, 0.23] | — | — | — | — |
| γⁿ 27.0–27.8 (near) | 0.13 [0.06, 0.20] | 0.09 [−0.03, 0.24] | 0.07 [−0.05, 0.25] | 0.11 [0.02, 0.19] | 0.48 [0.18, 0.78] |

**Scale-drift sign test (E16b).** For a drift DIC_obs = DIC_true·(1 + κ(t − 2000)):
- the single-layer excess implies κ = δ / DIC_target (positive: +2 to +8 × 10⁻⁵ yr⁻¹);
- the differential excess implies κ = Δδ / (DIC_target − DIC_ref) (negative: −0.4 to −4.6 × 10⁻³ yr⁻¹).

The signs conflict in 11 of 12 upper layers (the exception has δ ≈ 0). The full table is in `research/results/e16_scale_drift_sign.csv`.

## S-R4.2 Second model: NorESM2-LM (E15, pre-registered)

Samples ≤ 2014; 265,516 used (0.8 % without a TRACE value). Pooled upper-ocean ratios (global shared cruise bootstrap, 500 replicates):

| Comparison | All | Atlantic | Indian | Pacific | Southern Ocean |
|---|---|---|---|---|---|
| model DIC / C_ant^true | 1.01 [0.97, 1.04] | 0.98 | 0.92 | 1.04 | 1.00 |
| C_ant^true / TRACE (default) | 0.94 [0.91, 0.97] | 0.93 | 1.02 | 0.92 | 1.02 |
| C_ant^true / TRACE (unadjusted atm.) | 0.86 [0.84, 0.88] | 0.85 | 0.93 | 0.84 | 0.93 |
| MRI-ESM2-0 ≤ 2014, C_ant^true / TRACE (default) | 1.12 [1.10, 1.14] | — | — | — | — |

Differential scale factors m (all 20 layers): DIC vs truth −0.008 [−0.05, 0.02]; truth vs TRACE −0.02 [−0.05, 0.01] (default) and −0.11 [−0.13, −0.08] (unadjusted).

**Outcome against the pre-registered expectations.**
- H15a (estimator unbiased): met.
- H15b (shortfall > 5 % that the unadjusted boundary removes): not met.

**Descriptive comparison (post hoc).** At the same upper-ocean samples (≤ 2014), mean C_ant^true relative to TRACE is 1.01 in MRI-ESM2-0 and 0.85 in NorESM2-LM, and NorESM2-LM holds 0.90 × MRI-ESM2-0's C_ant. The difference in rate ratios between the models follows their ventilation. Full tables are in `research/results/e13_model_*_noresm2_lm*.csv`.
