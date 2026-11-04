# Final Scientific Review (Phase 11)

*Written as a skeptical reviewer for an oceanography journal (e.g., GRL, Biogeosciences, GBC). Covers `paper/MANUSCRIPT.md` and its evidence in `research/`. Date: 2026-11-04.*

## 1. Overall assessment
The work asks a sensible, falsifiable question — do observed DIC trends at fixed water properties match tracer-based C_ant accumulation? — and answers it with a transparent, pre-registered, sample-matched design and a validated estimator. The first-order result (tracer- and DIC-based rates agree within ±15 % across most upper and intermediate layers, with a clear ventilation structure) is solid but largely expected. The second-order result (a systematic 10–25 % upper-ocean excess of DIC-based over TRACE rates) is interesting and robust to analysis choices. It is, however, (i) below a conservative measurement-drift detection limit, (ii) supported only by a non-independent DIC-based comparator, and (iii) not yet tested against TRACE's own parameter uncertainty. **In its current form it reads as a careful methods-plus-diagnosis paper with a provocative but not yet conclusive finding.**

**Publication readiness: not ready to submit.** Estimated 2–4 weeks of focused work plus expert input stands between this and a credible submission (§5).

## 2. Strongest contributions
1. **A like-for-like design**: one estimator applied to observed DIC and to TTD C_ant at the same samples. This removes sampling-geometry differences that confound inventory-level comparisons.
2. **Within-cell identification** with density layers and property conditioning, validated by synthetic recovery (bias ≤ 0.012 µmol kg⁻¹ yr⁻¹) and leave-one-cruise-out.
3. **Robustness that is reported honestly**: 11 sensitivity variants, explicit detection limit, a documented coverage shortfall (≈91 %), and claims downgraded where results are fragile (deep Atlantic).
4. **A concrete, quantified correction of a common practice**: static DIC regressions validated on later years are confounded by C_ant uptake, and their model rankings depend on validation design.
5. **Reproducibility**: checksum-verified inputs, logged commands and seeds, pre-registration, versioned results.

## 3. Major weaknesses (likely reviewer objections)
| # | Weakness | Severity | Mitigation available? |
|---|---|---|---|
| W1 | Upper-ocean excess (2.5–4 µmol kg⁻¹ over the record) is ≤ the offset-equivalent detection limit; a coherent DIC calibration drift is not excluded. | **High** | Partial: post-1990 test passes. Add a time-resolved check against cruises with documented CRM use; test drift by era (pre-/post-2000 slopes within cells). |
| W2 | eMLR(C*) is built from the same DIC → "two DIC-based methods agree" is weaker than independent confirmation. | **High** | Compare with Müller et al. (2023) 2004–2014 product; compare with ocean-model C_ant (RECCAP2) sampled at the same points; use TTD with different tracer subsets (CFC-only vs SF₆). |
| W3 | TRACE run with default settings only. The excess could reflect its disequilibrium term or TTD shape (Δ/Γ), which the authors flag as uncertain. | **High** | Easy: rerun TRACE with `canth_diseq` ∈ {0.8, 1.2} and alternative `delta_over_gamma`. Report whether δ survives. |
| W4 | A single linear trend over 1972–2021 conflates accelerating C_ant uptake with era-dependent sampling. R is period-dependent. | Medium | Fit trend with an era-weighted or quadratic time term; compare like-for-like through TRACE (already sample-matched, which partly controls this). |
| W5 | AOU conditioning assumes AOU tracks respiration; biased in Southern Ocean source waters. | Medium | Phosphate variant already gives same sign; add C*-style conditioning with TAlk/nutrients. |
| W6 | Novelty search was limited (~15 queries, some publishers blocked). A sample-matched DIC-vs-TTD comparison may exist in regional studies (e.g., GO-SHIP section papers). | Medium | Systematic Scopus/WoS search; consult a carbon-cycle expert. |
| W7 | Bootstrap intervals ≈91 % coverage under large cruise offsets. | Low–Medium | Already disclosed; could widen to calibrated level. |
| W8 | Arctic, marginal seas, coasts excluded; Indian Ocean coverage thin (29–36 cruises per layer). | Low | Disclose (done). |
| W9 | Several references unverified ([M] tags); Keppler 2023 and Boteler 2023 only at abstract level. | Low (but must fix) | Verify before submission. |

## 4. Status of hypotheses
| Hypothesis | Outcome | Confidence |
|---|---|---|
| H1: R consistent with 1 ± 0.15 (γⁿ ≤ 27.8) | Supported in 18/20; contradicted in 2 Pacific intermediate layers | High (but many "supported" layers sit at the band edge) |
| H2: deviations concentrated in Southern Ocean / subpolar N Atlantic | **Partly**. The largest *robust* deviations are in the upper Indo-Pacific, not where expected. Southern Ocean shows sign reversals between methods. | Medium |
| H3: time/C_ant covariate removes ≥ 70 % of holdout bias | Supported for linear time (36/40 original QC; 30/40 strict); **not** for TRACE C_ant (16/40) — consistent with the main finding | High |
| H4: deep trends QC-sensitive | **Not supported** (0/14) | High |

## 5. Remaining experiments before submission (priority order)
1. **TRACE parameter sensitivity** (W3): disequilibrium and TTD-shape variants, sample-matched. *If δ disappears within TRACE's own uncertainty, the headline becomes "consistent within TRACE uncertainty".*
2. **Second C_ant comparator** (W2): Müller et al. (2023) 2004–2014 field and at least one ocean-model C_ant sampled at the same points.
3. **Drift diagnostics** (W1): within-cell slopes by era (1972–1995 / 1995–2021); cruise-level residual offsets vs year; CRM-era split.
4. **Time-shape sensitivity** (W4): quadratic time or TSS-shaped regressor (∝ atmospheric pCO₂ anomaly).
5. **Expert review and systematic literature search** (W6, W9).

## 6. Recommended target venues
* **Geophysical Research Letters**: if W1–W3 are resolved and the upper-ocean discrepancy survives. Short-format, high-visibility finding: "DIC-based and tracer-based anthropogenic carbon rates diverge in the upper ocean".
* **Biogeosciences** or **Ocean Science** (open access): if the discrepancy is only partly resolved. A full methods-plus-results paper on the sample-matched comparison and its uncertainties.
* **Limnology & Oceanography: Methods / Frontiers in Marine Science**: a separate, shorter note on the H3 finding (temporal validation of empirical DIC models confounded by C_ant).
* Conference: **Ocean Sciences Meeting** or **EGU General Assembly** abstract now — appropriate for preliminary results with stated caveats.

Not recommended: Nature-family journals. The finding is not yet strong or independent enough.

## 7. Concrete improvements required before submission
1. Run §5 items 1–3 and update results honestly, including if the excess vanishes.
2. Verify all [M]/[A] references against full texts; complete the systematic novelty search.
3. Have a carbon-cycle oceanographer review the interpretation (especially preformed DIC, disequilibrium, Southern Ocean).
4. Convert to the target journal's template; renumber figures; add a study-area map with coastlines (cartopy) and a schematic of the estimator.
5. Decide how to relate this work to the repository's original manuscript. It supersedes the "low-dimensional rank-1" claims: the original paper's temporal-validation conclusions should not be submitted as they stand.
6. Release code and derived data (Zenodo DOI), with the research dataset build documented.

## 8. Completion status (honest accounting)
| Phase | Status |
|---|---|
| 0 Safe environment, branch | Complete |
| 1 Repository intelligence | Complete for active pipeline; legacy code partially inspected (documented) |
| 2 GLODAP science | Complete (notes in `SCIENCE_AND_BASELINE.md`) |
| 3 Reconstruction/reproduction | Complete for parametric models; PySR search not re-run; repo test suite not executed (sandbox), replaced by independent re-implementation |
| 4 Literature review | Partial — targeted, not systematic; verification tags in place |
| 5 Proposals | Complete |
| 6 Methodology | Complete, pre-registered |
| 7 Implementation | Complete |
| 8 Validation | Complete for pre-registered set; §5 additions outstanding |
| 9 Interpretation | Complete in manuscript discussion; attribution remains inferential |
| 10 Publication package | Draft complete (`paper/*.md`); not in journal template; not expert-reviewed |
| 11 Final review | This document |
