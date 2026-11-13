# Gap Audit, Round 2 (pre-registered plan for the second iteration)

*Written 2026-11-14 on `research/dev` @ `e3dcab7`, **before** any new experiment was run. Later deviations are logged in `EXPERIMENT_LOG.md` under "Round 2".*

## 0. Reproducibility check on a second machine (done before writing this plan)
* Fresh `.venv` (Python 3.12, uv), TRACE-Python release v1.1.2 (GitHub zipball SHA-256 `01cf923e…0b51`, identical to Round 1).
* `build_dataset` → 512,384 relaxed / 428,264 strict rows (identical); TRACE → 400,755 finite C_ant of 402,747 (identical).
* `e02_layer_trends --variant primary` reproduces every committed number to ≤ 2·10⁻¹⁴ (incl. bootstrap CIs). **Round-1 results are reproducible across machines.**

## 1. Gaps found by auditing code, logs and results

| # | Gap | Evidence | Consequence |
|---|---|---|---|
| G1 | **Depth-uniform, basin-specific deviations in deep water.** Pacific γⁿ 27.5–28.2: δ = +0.047…+0.065 µmol kg⁻¹ yr⁻¹ where TRACE C_ant rates are ≈ 0.01; Atlantic 27.8–28.2: δ ≈ −0.06…−0.08. Stable in every sensitivity variant incl. post-1990. | `e02_layer_trends_*.csv` | A cruise-level calibration drift (offset correlated with year) produces exactly this depth-independent signature, and it would also inflate/deflate the upper-layer excess by the same amount. Round 1 neither modelled nor tested this. W1 is therefore more serious than stated, but also *testable*. |
| G2 | Synthetic test (E01) simulated cruise offsets *uncorrelated* with time. | `e01_synthetic_recovery.py` | Bias from a time-trending calibration drift was never tested; the estimator is not immune to it. |
| G3 | TRACE run only at defaults (TTD shape Δ/Γ = 1.30, disequilibrium scale 1.0, NN mean-age scale factors). | config `trace:` | W3: the excess may sit within TRACE's structural uncertainty. |
| G4 | AOU conditioning: in layers that deoxygenate (N Pacific intermediate), a biased AOU coefficient leaks respiration DIC into β. The phosphate variant (S10) does not remove carbonate-pump (alkalinity) changes. | S5/S10 spread up to 0.2 | Need a C*-type response (DIC corrected for soft-tissue and carbonate pumps using measured PO4 and TAlk) and an explicit "control-variable" check (TAlk trend at fixed θ, S). |
| G5 | 34 basin × layer groups tested separately; bootstrap draws not shared across groups. | `fit_trend` per group | No single pooled estimate of the upper-ocean excess with a valid joint CI; multiplicity not handled. |
| G6 | A single linear trend over 1972–2021; no test of whether the excess sits in early (pre-CRM, GEOSECS) or late data beyond the S7 cut. | config | W4. |
| G7 | Only one DIC-based comparator (Gruber 2019) and no download of it on this machine; no second C_ant estimate. | `data/external/` absent | W2. |
| G8 | Literature search ~15 queries; 9 references unverified ([M]). | `LITERATURE_REVIEW.md` | W6/W9. |
| G9 | Manuscript claims ("excess robust", "18/20 supported") pre-date G1–G6. | `paper/` | Must be revised to the new evidence, whichever way it points. |
| G10 | Study-area figure has no coastlines; no estimator schematic. | `fig1` | Presentation (minor). |

## 2. Pre-registered experiments (Round 2)

**E06 — Drift-immune differential trends (G1, G2; primary new test).**
For each basin, stack each target layer ℓ with a deep reference layer set r (γⁿ > 27.8 for Atlantic/Southern Ocean; > 27.8 for Indian/Pacific) restricted to cruises that sampled both. Fit

  y = α_{cell×layer} + γ_cruise + Σ_layer f_layer(θ, S, AOU) + Δβ_ℓ · (t − 2000)·1[ℓ] + ε

with two-way fixed effects (cell×layer and cruise) by alternating projections. Cruise fixed effects absorb any additive per-cruise offset, *including offsets that drift with time*; Δβ_ℓ is the trend of layer ℓ relative to the deep reference. Apply identically to y = DIC and y = TRACE C_ant; estimand Δδ_ℓ = Δβ_ℓ(DIC) − Δβ_ℓ(C_ant), cruise-cluster bootstrap (B = 300, seed 20261114).
*Interpretation rule (fixed now):* if Δδ_ℓ > 0 with CI excluding 0 in the light layers (γⁿ < 27.0), the upper-ocean excess is not explained by additive cruise-level calibration offsets or drift. If Δδ_ℓ ≈ 0 (CI includes 0 and |Δδ| < half the primary δ), the Round-1 excess is attributed to a basin-wide offset/drift or a deep non-C_ant signal, not to an upper-ocean process.
*Validation:* synthetic test with a time-correlated drift (cruise offset = d·(t_cruise − 2000), d ∈ {0, 0.05, 0.1} µmol kg⁻¹ yr⁻¹, plus random σ_c = 3) applied to both the single-layer (E02) and differential (E06) estimators. Pass: E06 bias < 0.02 for all d; E02 expected to be biased by ≈ d.

**E07 — TRACE structural sensitivity (G3).** Rerun TRACE at all samples with (a) Δ/Γ ∈ {1.0, 1.6} (Waugh et al. use 1.0; default 1.30), (b) disequilibrium scale `canth_diseq` ∈ {0.9, 1.1}, (c) NN mean-age scale factors × {0.8, 1.2}. Refit the E02 estimator for C_ant only (DIC unchanged). Report the range of δ across variants for each layer.
*Rule:* the upper-ocean excess is "within TRACE structural uncertainty" in a layer if any variant gives δ with CI including 0.

**E08 — Pump-corrected response and control variable (G4).**
(a) C* = DIC − r_C:P·PO4 − ½(TAlk + r_N:P·PO4), r_C:P = 117, r_N:P = 16 (Gruber et al., 1996; as in eMLR(C*)), using samples with measured, flag-2 PO4 and TAlk (TAlk secondary-QC 1); conditioned on θ, S only. Expectation: β_C* ≈ β_Cant(TRACE) + Δ(preformed/disequilibrium).
(b) TAlk at fixed θ, S as a control: a non-zero trend indicates alkalinity measurement drift or a real carbonate change.
(c) Trend of AOU at fixed θ, S (deoxygenation context), reported, not tested.

**E09 — Pooled upper-ocean estimate and time shape (G5, G6).**
(a) Global cruise bootstrap (cruises resampled once per replicate across all groups). Pooled sample-weighted R and δ over the 13 layers with γⁿ ≤ 27.0 (the "upper ocean"), and over 27.0–27.8.
(b) Piecewise-linear time (knot 2000) for DIC and C_ant with the same estimator: does the excess differ between 1972–2000 and 2000–2021?

**E10 — Second comparator (G7).** Download Gruber et al. (2019) (already checksummed) and, if publicly available as a gridded product, Müller et al. (2023) ΔC_ant 2004–2014. Compare rates at the same samples as E05. If not obtainable, record this.

**Literature (G8):** systematic search for sample-level/trend-level DIC-vs-TTD comparisons, TRACE validations, deep Pacific DIC change, and GLODAP DIC drift.

## 3. Decision on the paper
The headline claim will be rewritten to whatever E06–E09 support. In particular, if E06 removes the excess, the paper becomes "observed DIC accumulation agrees with TTD expectations once cruise-level drift is controlled; the apparent upper-ocean excess is a calibration artefact": still a publishable, useful result.

## 4. Post-hoc addition (logged before running): E12 — TRACE surface boundary condition
*Added 2026-11-14 after E06–E11, so it is **exploratory/post hoc**, not confirmatory.*
Observation motivating it: TRACE's built-in `CO2TrajectoriesAdjusted.txt` (surface-ocean boundary condition) grows at 0.883× the atmospheric rate over 1972–2021 (anomaly 0.89–0.91× atmospheric), i.e. it encodes a damped surface response. If the outcrop boundary condition actually tracked the atmosphere, TRACE rates would be ≈ 13 % higher — inside the E11 CI for m (0.11–0.23).
Test: rerun TRACE at every sample with the unadjusted atmospheric trajectory (`CO2Trajectories.txt`, column 1), refit E02 (per layer), E09 (pooled upper R) and E11 (drift-immune m).
Expectation if the hypothesis explains the excess: pooled upper-ocean R and (1 + m) fall to within ≈ 0.95–1.05 with CIs including 1.
Caveat fixed in advance: agreement would show the excess is *equivalent* to a ~12 % surface boundary-condition difference, not that the boundary condition is the cause (natural-carbon change or TTD mixing errors could produce the same pattern).
