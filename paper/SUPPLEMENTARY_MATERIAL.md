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
