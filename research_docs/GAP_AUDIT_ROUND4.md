# Round 4: Within-cruise depth-dependent artefacts (R4) and closing items (Pre-registration)

*Written 2026-11-14, before E16 was run. Addresses FINAL_REVIEW R4 (within-cruise, depth-dependent DIC artefacts not excluded) and R3 (post hoc boundary-condition mechanism).*

## 1. E16: near reference (R4), pre-registered

The cruise fixed effects remove any offset that is additive and the same at all depths on a cruise. They do not remove an artefact that differs between the upper ocean and the deep reference on the same cruise. Examples:
* a multiplicative (scale) calibration error, which gives a larger absolute offset at higher DIC;
* a depth- or pressure-dependent sampling or storage artefact.

Such an artefact grows with the depth and DIC separation between target and reference. E16 re-runs E11 (global shared cruise bootstrap, 500 replicates, seed 20261114) with a **near reference**:
* reference: the intermediate layers 27.0 < γⁿ ≤ 27.8 (two layers per basin);
* targets: the 12 layers with γⁿ ≤ 27.0.

This cuts the typical target–reference depth separation from several kilometres to a few hundred metres, and the DIC contrast by roughly half.

Under the proportional model (δ_l = m·β_Cant,l for every layer), the scale factor m is the same for any reference. An artefact that grows with separation would shrink m towards 0 when the reference is near.

**Decision rule.** The excess is judged not to come from separation-dependent artefacts if both hold:
* the all-basin m (near reference) has a 95 % CI excluding 0;
* its point estimate is within 0.08 of the primary upper-layer m = 0.18.

If the point estimate falls below 0.10, the manuscript says that part of the differential excess may be separation-dependent. Wider CIs are expected, because the TRACE rate contrast against an intermediate reference is smaller.

**E16b (analytic, same data).** Suppose a uniform scale drift DIC_obs = DIC_true·(1 + κ(t − 2000)) were the cause. The single-layer excess (E02/E09, not drift-immune) needs κ > 0. The differential excess (E06/E11) needs κ·(DIC_target − DIC_ref) > 0, so κ < 0, because upper-ocean DIC is lower than deep DIC. We report the κ each requires, using layer-mean observed DIC. If the two signs disagree, uniform scale drift cannot explain both estimates at once.

## 2. R3: out-of-sample test of the boundary-condition mechanism

E15(d) (GAP_AUDIT_ROUND3 §8) is pre-registered and serves as the held-out test: TRACE with the unadjusted atmospheric boundary condition in a second, independent model (NorESM2-LM). The prediction made before the data were seen is that the unadjusted variant brings R_true/TRACE within [0.95, 1.05].

## 3. Outcomes (filled in after the runs)

| Item | Outcome |
|---|---|
| E16 (near reference) | **Met.** All-basin m = 0.127 [0.065, 0.199]: the CI excludes 0 and the estimate is 0.055 below the primary 0.18 (rule: within 0.08). About 30 % smaller, so a minor separation-dependent part cannot be excluded. Basins: Atlantic 0.09, Indian 0.07, Pacific 0.11, Southern Ocean 0.48 (wide CIs) |
| E16b (scale-drift sign) | 11 of 12 upper layers need κ of opposite sign for the single-layer and differential excess (the exception, Atlantic γⁿ ≤ 26.0, has single-layer excess ≈ 0). The differential excess alone would need κ ≈ −0.04 to −0.46 % yr⁻¹, i.e., a 2–23 % scale drift over 50 years. Uniform scale drift is excluded |
| E15 / R3 | *pending* |
