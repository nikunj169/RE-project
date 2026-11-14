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
* Result (window 1989–2012): mechanical verdicts — TRACE_low 19, nonCant_or_drift 4, n/a (δ ≤ 0) 11. Among the 15 layers with robust δ > 0 in E02 (16 with CI > 0; Indian γ<26 fails sign stability): TRACE_low 13, nonCant_or_drift 2 (Indian 27.0–27.5, Southern Ocean 26.5–27.0). [Corrected 2026-11-04 after re-check; an earlier version of this entry said 16 and 3.]
* Estimator-independent comparison (eMLR vs TRACE rates, same samples, 1994–2007): eMLR exceeds TRACE by ~10–35 % in γⁿ < 27.0 in all basins and in Pacific/Indian intermediate layers; eMLR is **lower** than TRACE in Southern Ocean 27.0–27.8 (0.46 vs 0.58; 0.26 vs 0.37) and Atlantic deep layers.
* Within the 1989–2012 window, observed δ shrinks and CIs widen (fewer samples); only Pacific 27.0–27.5, Atlantic 26.0–26.5, Indian 26.5–27.0 and SO 26.5–27.0 keep CI excluding 0.
* Interpretation constraint: eMLR(C*) is derived from GLODAP DIC → not independent of the observations; the supportable statement is that two DIC-based estimates agree with each other and exceed the tracer-based (TTD/TRACE) estimate in the upper ocean outside the Southern Ocean, while the Southern Ocean intermediate layers show the opposite sign between eMLR and TRACE.

---
# Round 2 (2026-11-14) — plan in `research_docs/GAP_AUDIT_ROUND2.md`, written before E06–E11
Environment: macOS (Darwin 27.2), fresh `.venv` (Python 3.12, uv), numpy/pandas/scipy/xarray/netCDF4/cartopy, TRACE-Python release v1.1.2 (GitHub zipball SHA-256 `01cf923e…0b51`, identical to Round 1), PyCO2SYS 2.0.0b8. Seeds: 20261104 (E02-type fits), 20261114 (Round-2 bootstraps).

## R0 — Cross-machine reproduction (validation)
* `build_dataset` → 512,384 / 428,264 rows (identical); `e02a_run_trace` → 400,755 finite C_ant (identical); `e02_layer_trends --variant primary` reproduces all committed numbers incl. bootstrap CIs to ≤ 2·10⁻¹⁴. Repository test suite run for the first time: 22/22 pass (after installing sympy); 3 new estimator tests pass.
* `build_dataset.py` gained `nitrate_qc`, `phosphate_qc` columns (additive; all other columns and row order unchanged).

## E06 — Drift-immune differential trends (confirmatory, pre-registered in GAP_AUDIT_ROUND2 §2)
* Command: `python -m research.experiments.e06_differential_trends --boot 300` → `e06_differential_trends.csv`.
* Estimator: two-way FE (cell×layer, cruise); exact solver (`differential_trend.exact_twoway`, Schur complement over cruise effects) — verified against alternating projections to 1.8·10⁻¹⁴. (A first implementation with alternating projections was correct but ~110× slower; replaced before the B = 300 run.)
* Result (Δδ = excess of target layer relative to γⁿ > 27.8 reference, same cruises): Atlantic 26.0–26.5 +0.190 [0.137, 0.268], 26.5–27.0 +0.177 [0.110, 0.237], 27.0–27.5 +0.131 [0.087, 0.186]; Indian +0.126…+0.153 (γ < 27.0), +0.105, +0.064; Pacific +0.124, +0.120, +0.078, +0.070, +0.021; Southern Ocean +0.138 (CI incl. 0), +0.186, +0.122, +0.023, −0.016.
* Pre-registered verdicts: upper layers (γⁿ ≤ 27.0) — 10/12 "excess not explained by cruise offsets", 2 inconclusive (Atlantic and Southern Ocean surface, γ < 26). Intermediate — 4 not explained, 2 attributable to offset/deep signal (Pacific and SO 27.5–27.8), 2 inconclusive.

## E06s — Synthetic drift test (validation)
* Command: `... e06_differential_trends --synthetic --reps 60 --boot 150` → `e06_synthetic_drift.csv`. Cruise offsets = d·(t_cruise − 2000) + N(0, 3²), d ∈ {0, 0.05, 0.1}.
* E06 bias ≤ 0.002 µmol kg⁻¹ yr⁻¹ for all d; 95 % CI coverage 93–100 % → **passed**. The single-layer estimator's bias increases by exactly d (column `e02_target_bias_noiseless` rises 0.05 per 0.05 of drift; its absolute level also contains the real-data cruise component of the synthetic null, so only the increments are interpretable).

## E07 — TRACE structural sensitivity (sensitivity, pre-registered)
* Command: `python -m research.experiments.e07_trace_sensitivity --boot 300` → `e07_trace_sensitivity.csv` (TRACE re-run at all samples for 7 variants).
* Δ/Γ ∈ {1.0, 1.6} and θ-instead-of-T change δ by ≤ 0.015; mean age ×0.8/×1.2 by −0.04…+0.03; disequilibrium scale ×0.9/×1.1 by ±0.03–0.07 in upper layers (dominant knob).
* Rule "excess within TRACE structural uncertainty if any variant's δ CI includes 0": survives all 8 variants in Pacific (all layers ≤ 28.2), Atlantic 26.0–26.5, Indian 26.5–27.0, SO 26.0–26.5; 6–7/8 in other upper layers (the failing variant is always diseq × 1.1 or mean age × 0.8).

## E08 — C* response and control variables (sensitivity, pre-registered)
* Command: `python -m research.experiments.e08_cstar_controls --boot 300` → `e08_cstar_controls.csv`; summary in `summary_round2_cstar.csv`.
* C* = DIC − 117·PO4 − ½(TAlk + 16·PO4) on rows with flag-2/qc-1 PO4 and TAlk (≈ 30–70 % of each group). Upper layers: C* excess CI > 0 in 9/12 (same count as the AOU estimator on the same rows). Pacific 26.0–27.0 (largest AOU trends, +0.16 and +0.43 µmol kg⁻¹ yr⁻¹ at fixed θ,S) drops to +0.08/+0.02 with CIs incl. 0 → excess there depends on the remineralisation tracer.
* Control: TAlk at fixed θ,S has CI excluding 0 in 5/34 groups, all negative (−0.05…−0.08), all deep Pacific/Indian. No physical process lowers deep TAlk; interpreted as TAlk measurement drift (CRMs certified for TAlk only from 1996) → C* is not used at depth.

## E09 — Pooled estimates with a global cruise bootstrap; piecewise time (pre-registered)
* Command: `python -m research.experiments.e09_pooled_and_timeshape --boot 500` → `e09_pooled.csv`, `e09_piecewise_time.csv`.
* Upper (γⁿ ≤ 27.0, 12 groups): δ = 0.130 [0.090, 0.171], R = 1.163 [1.114, 1.213]. Atlantic 1.08 [1.02, 1.15], Indian 1.17 [1.07, 1.32], Pacific 1.19 [1.12, 1.28], SO 1.16 [1.03, 1.26]. Intermediate: δ = 0.042 [0.019, 0.068], R = 1.15 [1.07, 1.25]. Deep: δ = −0.012 [−0.036, 0.013] (Atlantic −0.065, Pacific +0.055).
* Piecewise (knot 2000): post-2000 upper-layer δ > 0 with CI excluding 0 in Atlantic 26.0–27.0, Pacific surface, SO 26.0–26.5 → the excess is present in the CRM era. Pre-2000 slopes are poorly constrained (CIs ±0.2–0.5); post-minus-pre differences significant only in Atlantic 26.5–27.0 (+0.31 [0.11, 0.57]). No resolvable era dependence overall.

## E10 — Two eMLR(C*) products at the same samples (pre-registered)
* Data: Gruber et al. (2019) (SHA-256 `a5dcf2b0…a4d`, identical to Round 1); Müller et al. (2023) NCEI 0279447 `dcant_emlr_cstar_mueller_94-04-14.nc` SHA-256 `6b1669c4dcd2deb7a90bd7b1f03d0cc6385caf7afec2d2421c6fddbd4d9d3aba` (standard case = mask "3", "Standard case"; 11 defined cases). Both read-only in `data/external/`.
* Command: `python -m research.experiments.e10_emlr_products` → `e10_emlr_products.csv`; pooled by `python -m research.evaluation.summarize_round2` → `summary_round2_emlr_pooled.csv`. (First version with lazy xarray point sampling stalled >40 min; replaced by in-memory nearest/linear sampling, numerically identical by construction.)
* Upper-ocean pooled eMLR/TRACE ratios (same samples, same periods, no regression on our side): Gruber 1994–2007 1.13 (ensemble 1.06–1.27); Müller 1994–2004 1.20 (1.01–1.25); Müller 2004–2014 1.13. Observed (E02 pooled) 1.16. Deep Atlantic: eMLR also below TRACE (0.80; 0.35 in 2004–2014), as is observed DIC (0.75). Deep Pacific: eMLR also above TRACE (absolute rates ≤ 0.09 vs 0.01).

## E11 — Proportional form of the drift-immune excess (post hoc summary of E06; logged as such)
* Command: `python -m research.experiments.e11_scale_factor --boot 500` → `e11_scale_factor.csv`. m = LS slope (through origin) of Δδ on the TRACE differential rate; global shared cruise bootstrap.
* m = 0.181 [0.114, 0.232] (all basins); Atlantic 0.20 [0.12, 0.32], Indian 0.19 [0.06, 0.32], Pacific 0.16 [0.06, 0.23], SO 0.18 [0.04, 0.29]. Deviations from proportionality: Indo-Pacific/Atlantic intermediate layers (TRACE rate 0.1–0.25) sit above the line; SO intermediate layers below.

## E12 — TRACE with the unadjusted atmospheric trajectory (post hoc, pre-logged in GAP_AUDIT_ROUND2 §4)
* Observation: TRACE's `CO2TrajectoriesAdjusted.txt` grows at 0.883× the atmospheric rate (1972–2021); its anomaly is 0.89–0.91× atmospheric.
* Commands: `python -m research.experiments.e12_trace_unadjusted --fit --boot 300`; `... e09_pooled_and_timeshape --cant atm`; `... e11_scale_factor --cant atm`.
* TRACE rates rise by 9–10 % (`rate_ratio_atm_default` 1.09–1.13). Upper pooled R = 1.059 [1.014, 1.105] (Atlantic 0.99, Indian 1.07, Pacific 1.09, SO 1.06); m = 0.075 [0.014, 0.122]. Deep: Atlantic R 0.69 [0.53, 0.89], SO 0.74 [0.53, 0.96].
* Pre-logged expectation (R within 0.95–1.05 with CI incl. 1) **not met**: the undamped boundary condition removes ≈ 60 % of the excess; a residual ≈ 6–8 % remains with CIs narrowly excluding 0.

---

# Round 3 (2026-11-14): independent model comparator, reference layers, full-text checks

Design and decision rules were written in `GAP_AUDIT_ROUND3.md` before any model data were sampled (§1–4). E13d was logged in §6 after seeing E13 and before running it.

## R5: Full-text reference checks
* **Bates et al. (2014):** read in full (Oceanography PDF). Table 2 surface pCO₂ trends: BATS 1.69 ± 0.11, HOT 1.72 ± 0.09, ESTOC 1.92 ± 0.24, Iceland Sea 1.29 ± 0.36, Munida 1.28 ± 0.33 µatm yr⁻¹. The text states the subtropical-gyre sites track the atmospheric increase. The manuscript sentence was revised with these values and the high-latitude caveat.
* **Hong & Zemskova (2026):** read from the EarthArXiv preprint v2 (doi:10.31223/X5QX88; the publisher page returns 403). It is an ECCO-Darwin C* storage analysis for 1992–2022, with no TTD comparison, no sample matching and no cruise-offset treatment. No overlap.

## E13: Perfect-model test in MRI-ESM2-0 (pre-registered)
* **Sampling.** Command: `python -m research.experiments.e13_sample_model` → `research/data/model_mri_esm2_0.csv.gz`; log `e13_sample_model.log`; about 35 min, streamed and not stored.
  * Stores: Pangeo CMIP6 Zarr, MRI-ESM2-0 r1i2p1f1 `gr` v20200303. Historical (1972–2014) and ssp585 (2015–2021): Oyr `dissic`, `dissicnat`, `o2`; Omon `thetao`, `so` (annual means).
  * Sampler: nearest column wet at the sample's depth, linear in depth. Depths below the deepest wet level (6150 m) are clamped (first run crashed on this; fixed).
  * Output SHA-256 `6b76b131e3d2e4bb202430459d6cc5c1ceee80204f8c639f7f1ab493b192b01c`. All 343,826 E02 samples are finite. Model vs observed: θ r = 0.975 (mean +0.8 °C); S r = 0.86; AOU r = 0.89 (−14); DIC r = 0.90 (−26 µmol kg⁻¹).
* **Analysis.** Command: `python -m research.experiments.e13_perfect_model --boot 500` → `e13_model_layer_trends.csv`, `e13_model_pooled.csv`, `e13_model_differential.csv`, `e13_model_scale_factor.csv`. TRACE on model S and in-situ T is cached in `research/data/trace_canth_model_mri.csv.gz`.
  * 6,177 rows (1.8 %; warm, salty Pacific/Atlantic thermocline, 0–700 m) return NaN from TRACE for model T, S and are excluded. 5,304 of these are finite for observed T, S.
* **(a) Estimator vs model truth.**
  * Upper pooled R = 0.970 [0.948, 0.982], δ = −0.027 [−0.047, −0.015]; intermediate R = 0.947; deep δ = −0.023.
  * Differential m = −0.011 [−0.045, 0.017]. By basin: Atlantic +0.002, Indian +0.016, Pacific −0.056 [−0.106, −0.030], Southern Ocean −0.009.
  * **H13a met** (point estimates within the band). The Pacific bias is negative, i.e. conservative for the observed excess.
* **(b) Model truth vs TRACE(model T, S).**
  * Upper pooled R = 1.116 [1.103, 1.138]. By basin: Atlantic 1.036, Indian 1.011, Pacific 1.182, Southern Ocean 1.064. Intermediate 1.104; deep 0.757.
  * Differential m = 0.125 [0.102, 0.154]. By basin: Atlantic 0.086, Indian 0.026, Pacific 0.221, Southern Ocean 0.210.
* **(c) Model DIC vs TRACE (analogue of the observed comparison).** Upper R = 1.082 [1.063, 1.099]; m = 0.116 [0.079, 0.153].
* **Reading.** The sign and about two thirds of the size of the observed shortfall (m 0.18) are reproduced without GLODAP DIC. The basin pattern differs from the observations, which show 16–20 % everywhere.

## E13d: TRACE with the unadjusted trajectory on model T, S (post hoc, logged before running)
* Command: `python -m research.experiments.e13_perfect_model --trace-variant atm` → `*_atm.csv`; TRACE cache `trace_canth_model_mri_atm.csv.gz`.
* **(b) result.**
  * Upper pooled R = 1.016 [1.004, 1.036]. By basin: Atlantic 0.944, Indian 0.920, Pacific 1.076, Southern Ocean 0.969. Deep 0.693.
  * m = 0.024 [0.003, 0.050].
* Logged expectation (upper R in [0.95, 1.05]) **met**: the boundary condition accounts for ≈ 85–90 % of the model-world shortfall, against ≈ 60 % of the observed one (E12).

## E14: Alternative reference layers (pre-registered)
* Commands: `python -m research.experiments.e11_scale_factor --ref upperdeep` and `--ref abyssal` → `e11_scale_factor_ref_upperdeep.csv`, `e11_scale_factor_ref_abyssal.csv`.
* Reference γⁿ 27.8–28.0: all-basin m = 0.162 [0.107, 0.206]. By basin: Atlantic 0.18, Indian 0.15, Pacific 0.15, Southern Ocean 0.17.
* Reference γⁿ > 28.0: all-basin m = 0.192 [0.124, 0.249]. By basin: Atlantic 0.19, Indian 0.21, Pacific 0.18, Southern Ocean 0.19.
* Primary (γⁿ > 27.8): 0.181. The decision rule (CI > 0, point within 0.08 of 0.18) is **met**. The Atlantic did not move as much as expected.

## E15: Second model, NorESM2-LM (pre-registered, GAP_AUDIT_ROUND3 §8). **In progress; blocked by bandwidth**

* **Like-for-like MRI control (done).** Command: `python -m research.experiments.e13_perfect_model --max-year 2014 [--trace-variant atm]` → `e13_model_*_to2014[_atm].csv`; 263,118 samples.
  * Upper true/TRACE: R = 1.117 [1.098, 1.144], m = 0.132 [0.093, 0.170]. With the unadjusted boundary: R = 1.020 [1.003, 1.044], m = 0.033.
  * Estimator: R = 0.962, m = −0.023.
  * Dropping 2015–2021 does not change the MRI conclusions.
* **NorESM2-LM sampling** (`python -m research.experiments.e13_sample_model --model noresm2_lm`).
  * First attempt: reached 2008 (about 60 s per year), then hung on a Google Cloud read with no timeout (log `e15_sample_noresm2_lm_attempt0_hung.log`). No checkpoints existed, so that progress was lost.
  * Added per-year checkpoints (`research/data/_cache_model_<name>/`), a local catalogue cache, and a stall watchdog.
  * Attempts 1–5 hung on the first year. Diagnosis (14:40 IST): the whole connection had dropped to about 0.03–0.3 MB/s (Google Cloud, AWS and Cloudflare alike), against about 8 MB/s during E13. Not a code fault.
  * Remaining download is about 45 GB. **Paused** pending bandwidth; the run resumes from checkpoints.
