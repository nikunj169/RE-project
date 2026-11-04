# Experiment Log

Environment for all entries unless stated: macOS (Darwin 25.6.0), `.venv` Python 3.12.16, numpy/pandas/scipy as installed by `uv` from `environment/requirements.txt`. Raw data: GLODAPv2.2023 merged master `.mat`, SHA-256 `002881fa71923d15c1bd5aca3b99bc220ad5d1c77a90f130ba1ff0d3910d8766` (NOAA NCEI accession 0283442, downloaded 2026-11-04). Commands are run from the repository root.

Status labels: **exploratory** (not used for inference) · **validation** (tests the method) · **confirmatory** (pre-registered test) · **sensitivity** · **failed/aborted**.

---

## P0 — Provenance check (validation)
* Command: inline script (see conversation record); logic reproduced in `research/preprocessing/build_dataset.py`.
* Result: repository QC applied to raw file → 512,384 rows; values identical to `data/processed/glodap_qc_filtered.csv` (TCO2, year exact). Provenance of the repository dataset confirmed.
* Of these rows: 35,819 DIC flag 0 (calculated); 67,839 DIC secondary-QC 0.

## P1 — Baseline reproduction (validation)
* Script: `scratchpad/repro_probe.py` (independent code; not repository code).
* Result: rank-1 validation RMSE 20.17/14.48/16.00/10.58 and external 17.67/18.41/13.82/8.99 (Atl/Ind/Pac/SO) — identical to repository. Full quadratic 23.41/13.73/15.43/10.25 vs repository 23.41/13.73/15.50/10.46 (difference attributed to ill-conditioning from unscaled AOU in repository).
* Probe: full quadratic + linear year term → external bias 0.47/1.41/1.72/−0.91 (vs 8.26/7.65/7.29/4.23 without); external RMSE 17.45/12.83/12.25/7.78 (vs 19.62/14.92/14.25/8.86); year coefficient 0.43/0.29/0.35/0.30 µmol kg⁻¹ yr⁻¹. **Exploratory**; motivates H3.

## P2 — Residual drift by depth (exploratory)
* Script: `scratchpad/drift_probe.py`. Naive pooled OLS slope of full-quadratic residuals vs year, by basin × depth bin.
* Result: ≈0.6–1.0 µmol kg⁻¹ yr⁻¹ (0–100 m) decreasing to ≈−0.06–0.2 (>3000 m). Confounded by sampling geography; motivates within-cell design.

## D1 — Research dataset build
* Command: `python -I research/preprocessing/build_dataset.py`
* Output: `research/data/glodap_research.csv.gz`, `research/data/build_stats.json`.
* Result: relaxed (repo) QC 512,384; strict QC 428,264. Strict drops: DIC calc/flagged 35,819; DIC no secondary QC 67,839; S/O2 flag/qc 49,351 (overlapping).

## E00 — Coverage / first-look within-cell trend (exploratory)
* Command: `python -I research/experiments/e00_feasibility.py` → `research/data/e00_feasibility.csv`
* Settings: 4°×4° cells, span ≥ 10 yr, linear conditioning on θ, S, AOU, no uncertainty.
* Result: 72–281 qualifying cells per basin × layer for γⁿ < 28.2; point trends ≈1.0 µmol kg⁻¹ yr⁻¹ (γⁿ < 26.5) declining to ≈0.0–0.3 (γⁿ > 27.8); deep Atlantic > deep Indo-Pacific. Two layers (γⁿ > 28.2, Indian/Pacific) too sparse.
* Use: feasibility only. Methodology written after this run (and before E01–E03).

## Methodology registered
* `research_docs/RESEARCH_METHODOLOGY.md` written 2026-11-04 before E01–E03.

---
## E01 — Synthetic recovery of the trend estimator (validation)
* Command: `python -m research.experiments.e01_synthetic_recovery --reps 100 --boot 200`
* Groups: Pacific γ 26.0–26.5, Atlantic 27.8–28.0, Southern Ocean 28.0–28.1, Indian 27.0–27.5. β_true ∈ {0, 0.1, 0.5, 1.0}; cruise-offset SD ∈ {0, 2, 4} µmol kg⁻¹. Seed 20261104.
* Output: `research/results/e01_synthetic_recovery.csv`, `e01_verdict.json`.
* Result: max |bias| 0.012 µmol kg⁻¹ yr⁻¹ → **bias criterion passed**. Coverage 82–99 %; 46/48 settings in [90, 99] → **coverage criterion failed as pre-registered** (Atlantic deep β=0 σ=4: 82 %; Pacific β=0.5 σ=4: 89 %).

## E01b — Coverage follow-up (deviation from protocol, logged)
* Reason: E01 coverage failure; decide interval method before inference.
* Command: `python -m research.experiments.e01b_coverage_followup --reps 200 --boot 200` (seed 777).
* Result: worst-case settings (σ_cruise = 4) — percentile bootstrap 90.5 / 91.0 / 91.0 %; CR1 86.5 / 86.0 / 86.5 %; union = bootstrap.
* Decision: keep pre-registered percentile bootstrap. State in the manuscript that nominal 95 % CIs have ≈ 91 % coverage under 4 µmol kg⁻¹ cruise offsets; results whose CI barely excludes a threshold are labelled fragile.

## E02a — TRACE evaluation
* Software: TRACE-Python package `trace-python` 1.1.1 (setup.py) from GitHub release tag v1.1.2 (zip SHA-256 01cf923ef5ba9068ab94b2e83cacaf2a9bc9e265db959dcee38d0964e94a0b51), PyCO2SYS 2.0.0b8; install approved by user 2026-11-04.
* Command: `python -m research.experiments.e02a_run_trace --chunk 25000`; settings: historical atmospheric trajectory (1), canth_diseq 1, EOS-80.
* Result: 402,747 strict-QC samples (study basins, finite γⁿ); Cant finite for 400,755 (1,992 NaN, 0.5 %, outside TRACE input validity). Runtime 18 s. Sanity check: subtropical surface Atlantic 51.7 (1995) → 72.4 (2015) µmol kg⁻¹.

## E02 — Primary layer trends, observed DIC vs TRACE Cant (confirmatory, H1/H2)
* Command: `python -m research.experiments.e02_layer_trends --variant primary` (B = 500, seed 20261104).
* Outputs: `research/results/e02_layer_trends_primary.csv/json`; summaries by `python -m research.evaluation.summarize` → `summary_h1.csv`, `summary_robustness.csv`, `summary_h4.csv`.
* Note: an earlier run of the same command before TRACE existed produced DIC-only output; it was overwritten by the joint run (β_DIC identical by construction except for the 0.5 % samples lacking Cant).
* **H1**: 20 basin × layer units with γⁿ ≤ 27.8 — 18 supported, 2 contradicted (Pacific 27.0–27.5: R = 1.68 [1.40, 2.25]; Pacific 27.5–27.8: R = 5.2 [2.9, 13.0]).
* **Systematic pattern (beyond H1 verdicts)**: δ = β_DIC − β_Cant > 0 with CI excluding 0 in light layers of all basins except the Atlantic surface (γ < 26): δ ≈ +0.10 to +0.17 µmol kg⁻¹ yr⁻¹ (R ≈ 1.13–1.25). Deep Pacific δ ≈ +0.05 (CI excludes 0). Southern Ocean bottom (γ > 28.2) δ = −0.066 [−0.107, −0.019]; Atlantic 28.0–28.1 δ = −0.072 [−0.125, −0.022].
* No |δ| exceeds the offset-equivalent detection limit (0.16–0.19 µmol kg⁻¹ yr⁻¹).
* Like-for-like check: estimator β_Cant vs direct TRACE dCant/dt agree within 0–0.06 µmol kg⁻¹ yr⁻¹; estimator is lower in deep layers (part of the Cant trend projects on θ/S conditioning) — reason for applying the identical estimator to both.

## E02-S — Sensitivity variants (sensitivity)
* Command: loop over S1, S2 (2/6/8°), S3 (5/15 yr), S4, S5, S6, S7, S9, S10, B1 with `--boot 300` (**deviation**: B = 300 instead of 500 for runtime).
* Robust (sign stable across all TRACE variants, primary CI excludes 0): light-layer excess in Atlantic 26.0–27.0, Indian 26.0–27.5, Pacific all layers, Southern Ocean 26.0–27.0; Southern Ocean bottom deficit; Atlantic 28.0–28.1 deficit (but the latter shrinks to −0.006 to −0.02 in S7 post-1990 and S10 phosphate conditioning → treated as **not robust** in interpretation despite meeting the mechanical rule).
* Light-layer excess unchanged when restricted to post-1990 data (S7) → argues against a pre-CRM calibration artefact.
* B1 (no conditioning) roughly doubles Pacific light-layer δ → conditioning on θ, S, AOU is necessary.
* **H4 not supported**: strict vs relaxed QC changes deep β_DIC by ≤ 0.023 vs CI half-widths 0.03–0.10 (0/14 layers QC-sensitive).

## E04 — Leave-one-cruise-out (S8)
* Command: `python -m research.experiments.e04_loco` → `research/results/e04_loco.csv`.
* Result: max change in β_DIC from removing any one cruise ≤ 0.065 (sparse SO surface); ≤ 0.03 elsewhere. No δ changes sign under LOCO in any layer where the primary CI excludes 0. Most influential in deep layers: GEOSECS-era cruises (316N19720718, 318M19730822).

## E03 — Repository hierarchy with/without time (H3), relaxed QC
* Command: `python -m research.experiments.e03_hierarchy_with_time --qc relaxed --covariates none year` (B = 1000 for paired CIs).
* Independent implementation reproduces repository rank-1 validation/external RMSE exactly (20.17/14.48/16.00/10.58; 17.67/18.41/13.82/8.99).
* Linear time term: rank-1 external bias 9.82→1.46, 9.77→2.04, 6.88→1.10, 4.33→−0.76 (Atl/Ind/Pac/SO); external RMSE improves 1.1–2.8 µmol kg⁻¹ in every family and basin.
* Model rankings: Atlantic temporal ordering (rank-1 best) persists with time term; blocked CV ordering differs (rank-2/full-quad better) → ranking depends on validation design, not on the trend.
* Paired cruise-cluster bootstrap: rank-1 vs rank-2 / full-quadratic differences within ±4 µmol kg⁻¹ in all basins and splits. Indian validation contains 3 cruises, external 4.

## E05 — Independent comparison with eMLR(C*) (added analysis; design fixed before computation)
* Purpose: discriminate the explanations for the light-layer excess (δ > 0): (i) TRACE underestimates Cant accumulation, (ii) non-Cant DIC change, (iii) measurement drift.
* Data: Gruber et al. (2019) ΔCant 1994→2007, NCEI accession 0186034 (doi:10.25921/wdn2-pt10), file `dcant_emlr_cstar_gruber_94-07_vs1.nc`, SHA-256 a5dcf2b0f8d58a35652bee6221df54104357705abf739b63febe031b50709a4d, stored read-only in `data/external/gruber2019/`. Standard case DCANT_01; ensemble spread over DCANT_01–14; SIGMA_DCANT_01.
* Design (per basin × γⁿ layer, same samples as E02 but restricted to 1989–2012 to bracket the eMLR period; inclusion rules unchanged):
  1. β_DIC and β_Cant(TRACE) from the identical within-cell estimator on the window.
  2. eMLR rate = mean over the group's samples of DCANT_01/13 at the nearest horizontal grid point, linearly interpolated in depth (ensemble min–max across DCANT_01–14 reported).
  3. TRACE rate over the same period = mean over samples of [Cant(2007) − Cant(1994)]/13 at fixed sample inputs.
* Decision rule (fixed in advance): in layers where δ > 0 robustly (E02), if eMLR rate ≥ β_DIC − 0.5·δ (i.e., closer to observed than to TRACE) → evidence for (i); if eMLR rate ≤ TRACE rate + 0.5·δ → evidence against (i), pointing to (ii)/(iii); otherwise inconclusive. Comparison 3 vs 2 tests TRACE vs eMLR independently of the DIC data.
* Caveats: eMLR(C*) has ~±20 % structural uncertainty (Clement & Gruber, 2018); period 1994–2007 differs from the full-record trend; eMLR is itself derived from GLODAP DIC (not fully independent of our observations, but uses a different estimator and C* rather than AOU conditioning).
* Command: `python -m research.experiments.e05_emlr_comparison --boot 300` → `research/results/e05_emlr_comparison.csv`.
* Result (window 1989–2012): mechanical verdicts — TRACE_low 19, nonCant_or_drift 4, n/a (δ ≤ 0) 11. Among the 16 layers with robust δ > 0 in E02: TRACE_low 13, nonCant_or_drift 3 (Indian 27.0–27.5, Southern Ocean 26.5–27.0, Atlantic 27.0–27.5 is not robust in E02).
* Estimator-independent comparison (eMLR vs TRACE rates, same samples, 1994–2007): eMLR exceeds TRACE by ~10–35 % in γⁿ < 27.0 in all basins and in Pacific/Indian intermediate layers; eMLR is **lower** than TRACE in Southern Ocean 27.0–27.8 (0.46 vs 0.58; 0.26 vs 0.37) and Atlantic deep layers.
* Within the 1989–2012 window, observed δ shrinks and CIs widen (fewer samples); only Pacific 27.0–27.5, Atlantic 26.0–26.5, Indian 26.5–27.0 and SO 26.5–27.0 keep CI excluding 0.
* Interpretation constraint: eMLR(C*) is derived from GLODAP DIC → not independent of the observations; the supportable statement is that two DIC-based estimates agree with each other and exceed the tracer-based (TTD/TRACE) estimate in the upper ocean outside the Southern Ocean, while the Southern Ocean intermediate layers show the opposite sign between eMLR and TRACE.
