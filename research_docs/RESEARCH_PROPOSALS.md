# Research Proposals (Phase 5)

All proposals are constrained by what the repository + GLODAPv2.2023 actually provide (Phase 1–2) and by the gap analysis in `LITERATURE_REVIEW.md`. Scores 1 (poor) – 5 (strong). “Leakage/confounding risk” is scored so that 5 = low risk.

## P1 — Sample-matched observational test of tracer-based anthropogenic-carbon accumulation rates
*Question:* At fixed water properties (θ, S, AOU) and within fixed locations, does interior DIC increase at the rate predicted by a transit-time-distribution (TRACE) Cant model, as a function of basin and neutral-density layer, over 1972–2021?
*Design:* identical within-cell fixed-effects estimator applied to (a) observed DIC, (b) TRACE Cant evaluated at the same samples and dates; compare slopes layer-by-layer with cruise-cluster uncertainty. Residual = DIC change not explained by Cant or respiration (natural carbon change, disequilibrium change, TSS/TTD error, residual offsets).
*Builds on repo:* uses the drift the repo mistook for error; reuses the basin framework.
*Main risk:* TRACE is a dependency (installation pending); deep-layer trends are near measurement-offset limits; attribution of residuals is ambiguous.

## P2 — Hidden non-stationarity in empirical DIC models: consequences for validation and model selection
*Question:* How much of the temporal-holdout error and model-ranking in static DIC regressions (S–T–AOU and richer) is caused by the unmodelled anthropogenic trend, and do conclusions about functional form (e.g., rank-1 vs rank-2) survive once it is modelled?
*Design:* re-run the repo hierarchy with/without a time (or TRACE-Cant) covariate under temporal, spatial-block, cruise-block and random splits; correct cluster bootstrap; equivalence vs measurement uncertainty.
*Risk:* methodological, partly known to specialists (ESPER/GLODAP mapped normalise to a reference year) → modest novelty alone.

## P3 — Decadal variation in DIC accumulation (1990s vs 2000s vs 2010s) by density layer
*Question:* Did the within-location DIC accumulation rate at fixed properties change between decades, and is any change consistent with Müller et al. (2023)’s reported 15 % decline in sink sensitivity?
*Risk:* era-specific sampling; few within-cell revisits per decade; overlaps directly with eMLR(C*) literature; low power.

## P4 — Southern Ocean abyssal / AAIW failure: AOU disequilibrium vs water-mass mixing
*Question:* Are the repo’s abyssal and AAIW failures caused by O2 air–sea disequilibrium (AOU ≠ TOU) or by mixing of distinct end-members? Test by replacing AOU with nutrients (PO4, NO3) / preformed-nutrient approaches.
*Risk:* mechanism already known (Ito et al., 2004); incremental; region-specific.

## P5 — Data-quality sensitivity of empirical carbon regressions (measured vs calculated DIC; secondary QC)
*Question:* How do GLODAP primary/secondary QC choices change empirical DIC model skill, biases and trend estimates?
*Risk:* a useful technical note; low scientific significance on its own.

## P6 — Symbolic regression with time/Cant as an input
*Question:* Does SR recover a compact law once time/Cant is included?
*Risk:* fashionable-method-driven; no clear scientific question beyond P2.

---

## Comparative scoring

| Criterion | P1 | P2 | P3 | P4 | P5 | P6 |
|---|---|---|---|---|---|---|
| Scientific significance | 4 | 3 | 4 | 3 | 2 | 2 |
| Novelty (vs literature found) | 4 | 3 | 2 | 2 | 2 | 2 |
| Data availability | 4 | 5 | 3 | 4 | 5 | 5 |
| Methodological feasibility | 4 | 5 | 3 | 4 | 5 | 3 |
| Statistical validity achievable | 4 | 4 | 2 | 3 | 4 | 2 |
| Oceanographic interpretability | 5 | 3 | 4 | 4 | 2 | 2 |
| Publication potential | 4 | 3 | 3 | 3 | 2 | 2 |
| Computational feasibility | 4 | 5 | 4 | 5 | 5 | 2 |
| Reproducibility | 4 | 5 | 4 | 4 | 5 | 3 |
| Leakage/confounding risk (5 = low) | 3 | 4 | 2 | 3 | 4 | 3 |
| **Total (/50)** | **40** | **40** | **31** | **35** | **36** | **26** |

## Selection

**Primary direction: P1, with P2 as an integrated secondary contribution and P5 as a mandatory sensitivity analysis.**

Justification:
* P1 and P2 tie on score, but P1 asks a question about the *ocean* (does observed carbon accumulation match tracer-based expectations, where and by how much does it deviate?) whereas P2 asks a question about *models*. P1 is oceanographically interpretable and falsifiable: TRACE makes an independent, quantitative, depth- and basin-resolved prediction that the DIC data can contradict.
* P2 falls out of P1 at low cost (the same time/Cant covariate) and directly corrects the repository’s published conclusions — this preserves continuity with the existing project.
* P5 is not a paper on its own but must be part of P1 (trend estimates at 0.1 µmol kg⁻¹ yr⁻¹ are vulnerable to offset and calculated-value artefacts).
* The exploratory coverage probe (E00) shows ≥ 10-year within-cell spans for 72–281 cells per basin × layer and a first-look within-cell DIC trend with a physically plausible vertical structure (≈1.0 µmol kg⁻¹ yr⁻¹ in light layers decaying to ≈0.0–0.3 in deep layers). This establishes feasibility only; it is not a result.

Not chosen: P3 (underpowered; would overlap eMLR(C*) studies with weaker data), P4 (known mechanism), P6 (method-driven).

**Dependency flagged:** P1’s primary comparison requires TRACE-Python (MIT; v1.1.2). Its installation was blocked by the session’s sandbox policy and needs user approval. Fallback if unavailable: compare observed rates with published layer/basin rates and with a TSS scaling argument — substantially weaker.
