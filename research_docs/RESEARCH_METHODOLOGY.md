# Research Methodology (Phase 6) — pre-registered before main experiments

**Status:** written 2026-11-04 after only the exploratory coverage probe E00 (no uncertainty, not used for inference). Any later deviation from this document is recorded in `EXPERIMENT_LOG.md` with a reason.

## 1. Title (working)
*Does interior-ocean carbon accumulate as tracer ages predict? A sample-matched test of transit-time-based anthropogenic CO2 rates against fifty years of GLODAP DIC.*

## 2. Motivation
Empirical DIC models built on hydrographic predictors (including this repository’s S–T–AOU regressions) are routinely validated on later years, yet DIC at fixed hydrography rises because of anthropogenic CO2 uptake. Independent tracer-based (TTD) methods such as TRACE predict *how fast* Cant should accumulate in each water mass under steady circulation. Observed DIC trends at fixed properties provide a direct, falsifiable check of those predictions — and deviations indicate either non-Cant (natural) DIC change, circulation change (TSS violation), or method bias.

## 3. Research question
Within fixed locations and at fixed (θ, S, AOU), is the observed rate of DIC change, by basin and neutral-density layer, consistent with the rate of Cant accumulation predicted by TRACE at the same samples?

## 4. Hypotheses
* **H1 (primary).** In layers with γⁿ < 27.8 kg m⁻³ (surface, thermocline, mode and intermediate waters — where the expected signal ≥ 0.3 µmol kg⁻¹ yr⁻¹ ≫ detection limit), the ratio R = β_DIC / β_Cant is consistent with 1. *Decision rule:* H1 is supported in a basin × layer if the 95 % CI of R overlaps [0.85, 1.15] (TRACE’s stated ±15 % structural uncertainty); “contradicted” if the CI lies entirely outside it.
* **H2.** Where H1 fails, deviations are concentrated in the Southern Ocean and the subpolar North Atlantic (pre-registered expectation from Steinfeldt et al., 2024; Gruber et al., 2019). Tested descriptively (sign, magnitude, location) — not a confirmatory test.
* **H3 (secondary, repository correction).** For static DIC regressions trained on < 2015 data, ≥ 70 % of the mean temporal-holdout bias (2015–17 and ≥ 2018) is removed by adding the TRACE Cant estimate (or a linear time term) as a covariate, and conclusions about functional form (rank-1 vs rank-2 vs full quadratic) change once the trend is modelled and dependence-correct uncertainty is used.
* **H4 (data quality).** Deep-layer (γⁿ ≥ 27.8) trend estimates change by more than their own CI half-width between strict and relaxed QC — i.e., deep trends are QC-sensitive.

## 5. Study regions and scales
* Basins: Atlantic, Indian, Pacific (north of 35 °S, GLODAP regions 1, 16, 8) and Southern Ocean (south of 35 °S), as in the repository. Arctic excluded (sparse DIC, TRACE validity poorly constrained).
* Vertical: neutral-density layers, edges γⁿ = [−, 26.0, 26.5, 27.0, 27.5, 27.8, 28.0, 28.1, 28.2, −] kg m⁻³ (density layers avoid aliasing isopycnal heave as a trend).
* Spatial unit: 4° × 4° cells (sensitivity 2°, 6°, 8°).
* Time: 1972–2021, decimal year.
* Inference target: a layer-mean within-cell trend over the full period (not interannual variability; not forecasting).

## 6. Data and QC
* Source: GLODAPv2.2023 merged master file, SHA-256 `002881fa…8766`, read-only.
* **Strict QC (primary):** DIC WOCE flag 2 and secondary-QC 1; salinity and oxygen flag 2 and qc 1; finite θ, pressure, γⁿ; repository range filters retained. → 428,264 samples (before layer/cell restrictions).
* **Relaxed QC (ablation):** repository range-only QC (512,384).
* Inclusion: cells whose observations span ≥ 10 years within a basin × layer (sensitivity 5 and 15 years). Basin × layer groups with < 500 samples or < 10 qualifying cells are not reported.
* Built by `research/preprocessing/build_dataset.py` (checksum-verified).

## 7. Estimator
For each basin b and layer ℓ:

  y_i = α_c(i) + f_bℓ(θ_i, S_i, AOU_i) + β_bℓ · (t_i − 2000) + ε_i

* α_c: cell fixed effects (within-cell identification — removes changes in where cruises sampled).
* f: full quadratic in standardised θ, S, AOU (10 terms incl. interactions; sensitivity: linear).
* Estimated by OLS after within-cell demeaning (Frisch–Waugh–Lovell).
* Run with **y = observed DIC** and, at the same samples, with **y = TRACE Cant** (`tracepy.trace`, historical atmospheric trajectory, default disequilibrium, inputs lon/lat/depth/S/T/decimal year). Because TRACE Cant depends on T and S, applying the identical estimator yields β_Cant: Cant accumulation at fixed properties, like-for-like.
* Estimands per basin × layer: β_DIC, β_Cant, Δβ = β_DIC − β_Cant, R = β_DIC / β_Cant.

## 8. Uncertainty quantification
* **Cruise-cluster bootstrap** (primary): resample cruises *with replacement, preserving multiplicity* (fixing the repository’s `np.isin` bug), B = 500, seed 20261104; DIC and Cant fits use the same resample so Δβ and R get joint CIs (percentile 95 %).
* Cluster-robust (CR1, cruise) standard errors as an analytic cross-check.
* TRACE structural uncertainty: ±15 % band on β_Cant (decision rule in H1).
* **Detection limit** per layer: δ_det = 4 µmol kg⁻¹ / (median within-cell time span) — trends below it are labelled “below offset-equivalent detection limit”.

## 9. Estimator validation (must pass before inference)
* **E01 synthetic recovery:** replace DIC by f̂(θ,S,AOU) + α̂_c + β_true·(t−2000) + cruise offsets ~ N(0, σ_c²) + residual noise resampled from the fitted model; β_true ∈ {0, 0.1, 0.5, 1.0}; σ_c ∈ {0, 2, 4} µmol kg⁻¹; 100 replicates per setting. Pass: |mean bias| < 0.02 µmol kg⁻¹ yr⁻¹ and bootstrap CI coverage within 90–99 %.
* **Like-for-like check:** β_Cant from the estimator vs the direct TRACE dCant/dt (Cant(t+½) − Cant(t−½) at fixed inputs) averaged over the same samples — they should agree if θ,S conditioning does not distort the Cant trend.

## 10. Baselines and comparisons
* Baseline 1: naive pooled regression (no cell FE) — shows sampling-composition bias.
* Baseline 2: depth-bin instead of density-layer analysis.
* Repository baseline (H3): mean, linear, rank-1, rank-2, full-quadratic, cubic regressions per basin, with and without (a) linear year, (b) TRACE Cant covariate.

## 11. Sensitivity analyses and ablations
| ID | Change | Purpose |
|---|---|---|
| S1 | relaxed QC | H4 |
| S2 | cell size 2°, 6°, 8° | spatial identification |
| S3 | min span 5 / 15 yr | inclusion rule |
| S4 | f linear | conditioning flexibility |
| S5 | drop AOU from f | role of respiration conditioning |
| S6 | drop θ,S from f | role of physical conditioning |
| S7 | exclude pre-1990 | early-data / calibration-era influence |
| S8 | leave-one-cruise-out influence (top-10 influential cruises) | single-cruise leverage |
| S9 | add seasonal harmonic (sin/cos of month) | seasonal sampling bias (upper layers) |
| S10 | nutrients (PO4) instead of AOU | AOU disequilibrium (Southern Ocean) |

## 12. H3 validation design (repository correction)
Splits: temporal (train < 2015; val 2015–17; ext ≥ 2018), spatial-block (5°×10°), cruise-block — GroupKFold(5); the repository’s models refit per fold (the frozen SR candidate is excluded because it cannot be refit). Metrics: RMSE, bias, MAE; paired per-cruise mean |error| differences with correct cruise-cluster bootstrap; equivalence judged against ±4 µmol kg⁻¹ (DIC measurement consistency), not ±1.

## 13. Leakage and confounding controls
* No random splits for skill claims. Trend estimation uses all years by design (it is inference, not prediction).
* Time–space confounding → cell fixed effects + density layers.
* Cruise offsets → cruise-cluster inference + detection limit + LOCO.
* Same-sample comparison (TRACE evaluated at the observed samples, not on a grid).
* AOU disequilibrium → S10.
* Interpretation discipline: β_DIC is an *association* at fixed properties; Δβ ≠ 0 is attributed to “non-Cant DIC change and/or TTD/TSS error and/or residual offsets” unless a sensitivity uniquely discriminates.

## 14. Expected contribution
1. A sample-matched, density-resolved observational test of tracer-based Cant accumulation rates over the full GLODAP record.
2. A map (basin × layer) of where observed carbon accumulation departs from steady-state TTD expectations.
3. A quantified correction of the repository’s published conclusions on temporal transferability and functional form.

## 15. Software and reproducibility
Python 3.12.16 (.venv), numpy/pandas/scipy; TRACE-Python v1.1.2 (zip SHA-256 `01cf923e…0b51`; install pending approval). All scripts under `research/`; every run appends to `EXPERIMENT_LOG.md` with command, seed, data hash and outputs.
