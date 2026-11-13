# Final Scientific Review: Round 2 (Phase 11)

*Written as a sceptical reviewer for GRL, GBC or Biogeosciences. Covers `paper/MANUSCRIPT.md` (revision 2) and its evidence in `research/`. Date: 2026-11-14. The Round-1 review (2026-11-04) is in git history at `e3dcab7`.*

## 1. What changed since Round 1

| Round-1 weakness | Round-2 response | Status |
|---|---|---|
| W1: excess below an offset "detection limit"; drift not excluded | New estimator with cruise fixed effects. Immune to any additive per-cruise offset incl. time-trending drift; validated synthetically (bias ≤ 0.002). Excess survives in 10/12 upper layers. The old detection-limit argument is superseded. | **Largely resolved.** Residual risk: within-cruise depth-dependent artefacts; reference-layer assumption |
| W2: only one, non-independent comparator | Second eMLR(C*) product (Müller 2023), two periods, 11 cases. All three product-periods exceed TRACE in the upper ocean (1.13–1.20; even the minimum case ≥ 1.01) | **Partly resolved.** Still all DIC-based; no ocean-model or non-GLODAP comparator |
| W3: TRACE at defaults only | 7 structural variants + unadjusted boundary condition (E07, E12). TTD shape and θ are irrelevant; ages ±20 % are minor; the disequilibrium/boundary condition dominates and explains about 60 % of the excess | **Resolved**, with a mechanistic lead |
| W4: single linear time | Piecewise (knot 2000): excess present post-2000; pre-2000 unresolvable | **Partly resolved** |
| W5: AOU conditioning | C* (PO₄ + TAlk) confirms 9/12 upper layers; Pacific thermocline is tracer-dependent; TAlk drift found at depth | **Resolved**, with an honest exception |
| W6/W9: literature and references | All [M] references verified on Crossref; TRACE-Python full text checked (inventory-level validation only); 8 new searches. One new neighbour (Hong 2026) not accessible | **Partly resolved** |
| Multiplicity (34 groups) | Pooled estimates with a global shared cruise bootstrap | **Resolved** |
| Reproducibility | Bit-for-bit reproduction on a second machine; repo tests run (22/22) + 3 new | **Resolved** |

## 2. Overall assessment
The paper now makes a clear, quantitative and falsifiable claim. Upper-ocean DIC at fixed water properties accumulates **16 % [11–21 %]** faster than TRACE predicts at the same samples, and **18 % [11–23 %]** faster relative to the deep ocean once all cruise-level offsets are removed. The excess is consistent across the four basins and matched by two eMLR(C*) products.

The drift-immune design is a genuine methodological contribution: a reusable way to compare any interior C_ant product with observations without mapping and without the calibration-drift confound. The boundary-condition result (E12) offers a concrete, testable mechanism for most of the discrepancy.

The main residual weakness is independence. All corroborating evidence comes from GLODAP DIC. An observation-independent check (ocean-model C_ant sampled identically, or TTDs run on different tracer subsets) would make the paper considerably stronger.

**Publication readiness:** close to submittable to *GBC* or *Biogeosciences* after the items in §5; a *GRL* submission is plausible if items 1–3 hold up. Not ready today. The manuscript is Markdown, needs expert review, and one closely related paper (Hong 2026) has not been read.

## 3. Strongest contributions
1. **Drift-immune, sample-matched comparison** of observed DIC trends with TTD C_ant rates (the cruise-FE differential estimator). It is novel as far as a targeted search can tell, validated, and fast (exact solver).
2. **A quantitative observational constraint** on TRACE's upper-ocean rate scale (1 + m = 1.18 [1.11, 1.23]), consistent across basins.
3. **Attribution of most of the discrepancy** to a specific, inspectable model component (the damped surface CO₂ boundary condition). The pre-logged prediction is reported as only partly met.
4. **Deep-ocean diagnostics.** The deep Atlantic deficit (observations and both eMLR products below steady state) agrees with the literature on non-steady ventilation. The deep Pacific excess is identified as the likely signature of residual offsets.
5. **Data-quality findings useful beyond this paper:**
   - TAlk drift at depth, which matters for C*-based methods;
   - QC-insensitivity of DIC trends;
   - correction of a widespread validation practice for empirical DIC models.

## 4. Remaining weaknesses (likely reviewer objections)

| # | Weakness | Severity | Mitigation |
|---|---|---|---|
| R1 | All supporting estimates use GLODAP DIC; no independent C_ant (ocean model, other tracers) | **High** | Sample RECCAP2/CMIP6 model C_ant (or OCIM) at the same points and dates; apply the identical estimator to model DIC (a perfect-model test of the estimator *and* of TTD) |
| R2 | The reference-layer assumption: Δδ is relative to γⁿ > 27.8; the deep Atlantic deficit inflates Atlantic Δδ | Medium | Already bracketed (E02 vs E06). Add an alternative reference (e.g., 27.8–28.0 only, or deep Indian/Pacific only) as a sensitivity |
| R3 | E11 and E12 are post hoc | Medium | Labelled as such. Ideally re-test the boundary-condition hypothesis on held-out data (e.g., fit on 1972–2005 cruises, test 2006–2021) or with a physically derived surface pCO₂ history (SOCAT-based) instead of the atmosphere |
| R4 | Within-cruise, depth-dependent DIC artefacts not excluded | Medium–low | Argue from coulometric/CRM practice (offsets are scale or additive); a check against cruises with documented CRM corrections |
| R5 | Bates 2014 rate comparison and Hong 2026 content not verified in full text | Medium (must fix) | Read both |
| R6 | Interval coverage ≈ 91 % (single-layer) under large offsets | Low | Disclosed; the differential estimator achieves 93–100 % |
| R7 | Pacific thermocline excess depends on AOU vs PO₄ conditioning | Low–medium | Disclosed. Could add a nitrate-based variant or preformed-nutrient approach |
| R8 | Era dependence unresolved before 2000 | Low | Disclosed |

## 5. Concrete steps before submission (priority order)
1. **Independent comparator (R1):** sample one ocean-model historical run (DIC and C_ant tracers) at GLODAP sample locations and dates; run E02/E06 on model DIC vs model C_ant (estimator check) and vs TRACE applied to model T, S (method check).
2. **Read Hong et al. (2026) and Bates et al. (2014) full texts**; confirm novelty and the surface-pCO₂ statement.
3. **Out-of-sample test of the boundary-condition hypothesis (R3).**
4. **Alternative reference layers (R2).**
5. Expert review by a carbon-cycle oceanographer, ideally a TRACE or eMLR(C*) developer; offer them the diagnostic.
6. Convert to the journal template (LaTeX/Word); renumber figures (7–11 become main-text 3–6; Round-1 Figs. 3, 4 and 6 move to the supplement).
7. Zenodo release of code and derived results.

## 6. Hypothesis status

| Hypothesis | Outcome |
|---|---|
| H1 (R within ±15 %, γⁿ ≤ 27.8) | Supported in 18/20 per layer, but pooled upper-ocean R = 1.16 [1.11, 1.21] sits at the band edge |
| H2 (deviations concentrated in SO / subpolar N Atlantic) | Not as stated. The largest robust excess is in the upper Indo-Pacific and Atlantic thermocline; deficits are in the deep Atlantic and SO bottom/intermediate water |
| H3 (time covariate removes ≥ 70 % of holdout bias) | Supported (36/40) |
| H4 (deep trends QC-sensitive) | Not supported |
| E06 rule (excess not from cruise offsets) | Supported in 10/12 upper layers |
| E12 (post hoc; boundary condition explains the excess) | Partly: ≈ 60 % explained; the pre-logged expectation of full explanation not met |

## 7. Recommended venues
* **Global Biogeochemical Cycles** (best fit): methods plus a global result with uncertainty.
* **Geophysical Research Letters**: if step 1 confirms the excess independently. Short headline: "Tracer-based estimates underestimate upper-ocean carbon accumulation by ~15–20 %".
* **Biogeosciences / Ocean Science** (open access): equivalent to GBC.
* **Conference:** Ocean Sciences Meeting 2028 or EGU 2027 abstract now (results are mature enough for an abstract).
* A separate short note (L&O: Methods) on the cruise-fixed-effects trend estimator and on temporal validation of empirical DIC models.

## 8. Completion status

| Item | Status |
|---|---|
| Round-1 phases 0–11 | Complete (see git history) |
| Round-2 gap audit and plan | Complete; written before experiments |
| E06–E10 (pre-registered) | Complete |
| E11, E12 (post hoc, logged) | Complete |
| Manuscript, abstract, captions, supplement, bibliography | Revised (revision 2) |
| Independent model comparator | **Not done** (R1) |
| Full-text checks of 2 references | **Not done** (R5) |
| Expert review | **Not done** |
