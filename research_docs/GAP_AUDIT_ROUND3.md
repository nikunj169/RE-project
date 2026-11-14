# Round 3: Independent Model Comparator (Pre-registration)

*Written 2026-11-14, before any model data were sampled. Addresses FINAL_REVIEW (Round 2) R1: every supporting estimate so far uses GLODAP DIC.*

## 1. Design (E13, perfect-model test)

**Model.** MRI-ESM2-0, CMIP6 `historical` r1i2p1f1 (1972–2014) plus `ssp585` r1i2p1f1 (2015–2021). Regular 1° grid (`gr`), Pangeo Google Cloud Zarr copies, version v20200303 (ssp585 version recorded in the log).

**Variables.**
* Annual (`Oyr`): `dissic`, `dissicnat`, `o2`.
* Monthly (`Omon`): `thetao`, `so`, averaged to calendar-year means in memory; nothing large is stored.

The model's true anthropogenic carbon is C_ant^true = `dissic` − `dissicnat`. CMIP6 defines `dissicnat` as DIC under pre-industrial atmospheric CO₂ with the same circulation, so C_ant^true is the model's exact anthropogenic carbon.

**Sampling.** Every strict-QC sample used in E02 (343,826 rows, 1972–2021) is sampled at its own year, position and depth:
* horizontally at the nearest model column that is wet at that depth;
* vertically by linear interpolation.

Unit conversions use GSW:
* mol m⁻³ → µmol kg⁻¹ via in-situ density;
* AOU = O₂sol(S, θ) − O₂;
* in-situ temperature (for TRACE) from θ.

Basin, layer and cell assignments are the observed ones, so the sampling design is identical to the real data.

**Three tests on model data:**

| Test | What it compares | What it tests |
|---|---|---|
| (a) Estimator check | β_DIC (model DIC; θ, S, AOU conditioning; identical estimator, E09 pooled bootstrap) vs β_Cant^true | Do model climate variability, model AOU changes, and the conditioning recover the true C_ant rate? |
| (b) TRACE vs model truth | TRACE run on model T, S at the samples vs β_Cant^true | Is TRACE's rate scale below a model whose surface CO₂ follows the atmosphere? Confounded by model–real ventilation differences (TRACE ages are trained on real-ocean tracers) |
| (c) Differential form | E06-style cruise-FE differential (vs γⁿ > 27.8) of (a) and (b) | Same questions as (a) and (b), with the cruise fixed effects |

## 2. Pre-registered expectations and decision rules

* **H13a (estimator).** Pooled upper-ocean (γⁿ ≤ 27.0) R_model = β_DIC / β_Cant^true lies within [0.95, 1.05], and |δ_model| ≤ 0.03 µmol kg⁻¹ yr⁻¹.
  * If it holds, the observed 16 % excess is not an artefact of the conditioning estimator under realistic climate variability.
  * If it fails (|R_model − 1| > 0.05), the observed excess must be corrected by the model bias, and the headline is revised.
* **H13b (TRACE in a model world).** No directional prediction. We report R_TRACE,model = β_Cant^true / β_TRACE(model T, S), pooled by layer set.
  * A value near the observed 1.16 would be read as consistent with a TRACE rate-scale bias that is independent of GLODAP DIC.
  * A value near 1 would weaken the independent support (but not the observational result).
  * Either way it is confounded by model ventilation, which is stated as a limitation.
* **Secondary.** Deep layers (γⁿ > 27.8): model β_DIC − β_Cant^true is reported as an empirical estimate of non-C_ant DIC change (circulation and biology) at depth.

## 3. Out of scope this round
* A second model (NorESM2-LM, curvilinear, 2.5× larger). It is optional if E13 is informative.
* Model CFC-based TTDs.

## 4. E14: alternative reference layers (R2), pre-registered

The drift-immune excess (E06/E11) is measured relative to all layers with γⁿ > 27.8. Two alternatives re-run E11 (global shared cruise bootstrap, 500 replicates) with the reference restricted to:
* **ref-upperdeep:** (27.8, 28.0] only. Avoids bottom waters, including the deep Atlantic deficit.
* **ref-abyssal:** γⁿ > 28.0 only.

**Decision rule.** The scale-factor finding is reference-robust if the all-basin m has a CI excluding 0 under both alternatives, and its point estimate stays within 0.08 of the primary 0.18. The Atlantic m is expected to move most, because the deep Atlantic deficit sits in the reference.

## 5. Literature items closed this round (R5)
* **Bates et al. (2014):** full text read; Table 2 confirms the manuscript statement for BATS, HOT and ESTOC. Sentence revised with values and the high-latitude caveat.
* **Hong & Zemskova (2026):** read from the EarthArXiv preprint v2. It is an ECCO-Darwin C* storage analysis with no overlap with the sample-matched, drift-immune comparison.

## 6. E13d (post hoc, logged after seeing E13 (a)–(c) and before running)
E13(b) showed the model's true upper-ocean C_ant rate exceeds TRACE-on-model-T/S by ≈ 12 % (pooled R) and m ≈ 0.125 (differential). Test: re-run TRACE on model T, S with the unadjusted atmospheric CO₂ trajectory (E12's `unadjusted_trajectory`), then repeat (b) and (c).

**Expectation.** If the damped boundary condition is the main cause, the pooled upper R_true/TRACE-atm falls within [0.95, 1.05]. Because E12 explained only ≈ 60 % of the observed excess, a partial reduction is the more likely outcome. We report whichever occurs.

## 7. Outcomes (filled in after the runs)

| Item | Outcome |
|---|---|
| H13a (estimator) | **Met.** Upper R = 0.97 [0.95, 0.98]; m = −0.01 [−0.04, 0.02]. The Pacific is biased low (conservative) |
| H13b (TRACE in a model world) | Model truth exceeds TRACE: upper R = 1.12, m = 0.125. Same sign and ≈ ⅔ of the observed size, but a different basin pattern |
| E13d (post hoc) | Expectation **met**: R = 1.016 [1.004, 1.036], m = 0.024 |
| E14 (reference layers) | **Met**: m = 0.16 and 0.19 against the primary 0.18 |
| R5 (literature) | Both references read in full; manuscript wording adjusted for Bates; no overlap from Hong & Zemskova |
