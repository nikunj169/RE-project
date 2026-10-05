# Final Reviewer Risk Report

## Reviewer 1 — “Is this overfitting or symbolic-regression cherry-picking?”

**Strongest criticism:** The symbolic library contains 599 expressions evaluated on one temporal validation period, and the best candidate was selected from that same period.

**Quantitative answer:** Selection used 2015–2017 validation only; the post-2018 external holdout was excluded programmatically from candidate selection and blocked validation. The symbolic search produced 599 retained candidates across 48 successful searches. Rank-1 recurrence was evaluated cell-by-cell rather than reported as a universal discovery: it appeared in 24/48 basin–operator–seed cells. The frozen symbolic candidates were selected before spatial/cruise evaluation.

**Remaining limitation:** The temporal validation set is still a model-selection set, so validation-selected symbolic RMSE is optimistic. The external holdout is the relevant final evaluation. The family classifier is heuristic and candidates were not fully algebraically deduplicated.

**Manuscript location:** Methods, “Symbolic regression and recurrence” and “Temporal split”; Results, “Rank-1 recurrence and operator-space dependence” and “Temporal, spatial, cruise, and external generalization.”

## Reviewer 2 — “Are GLODAP observations independent?”

**Strongest criticism:** GLODAP observations are clustered by cruise, station, section, location, and depth, so observation-level inference is anticonservative.

**Quantitative answer:** Spatial GroupKFold used 5° latitude × 10° longitude blocks; cruise GroupKFold used the available cruise identifier. Train/test group disjointness was asserted for every fold. Dependence-aware paired-error bootstrap results were reported separately from naive TOST. Fold-level RMSE SDs are retained in `spatial_cruise_validation_summary.csv`.

**Remaining limitation:** Spatial blocks do not guarantee complete independence because nearby hydrographic structure can cross block boundaries. Cruise identifiers group sampling events but do not necessarily capture all section-level dependence. Confidence intervals remain conditional on the chosen grouping design.

**Manuscript location:** Methods, “Blocked validation” and “Dependence-aware inference”; Results, “Temporal, spatial, cruise, and external generalization.”

## Reviewer 3 — “Does the result generalize geographically?”

**Strongest criticism:** Temporal validation can interpolate across repeated geographic sampling patterns and does not establish geographic extrapolation.

**Quantitative answer:** Spatial and cruise blocked validation were performed on pre-2018 observations. For example, Atlantic spatial RMSE was 16.30 for rank-1, 15.46 for rank-2, and 15.55 for full quadratic; Atlantic cruise RMSE was 16.91, 15.58, and 15.63, respectively. The generalization summary reports temporal, spatial, cruise, and locked external metrics side-by-side.

**Remaining limitation:** The blocked designs test transfer across the available GLODAP spatial/cruise groups, not arbitrary global extrapolation. Rankings change across regimes, so geographic universality is not claimed.

**Manuscript location:** Methods, “Blocked validation”; Results, “Temporal, spatial, cruise, and external generalization”; Discussion, “Domain dependence and generalization.”

## Reviewer 4 — “Why should rank-1 matter if more complex equations predict better?”

**Strongest criticism:** The best symbolic candidates achieve lower RMSE, so rank-1 may simply be an inferior approximation.

**Quantitative answer:** The best symbolic candidates improve temporal-validation RMSE relative to rank-1 by approximately 5.8%, 11.6%, 3.2%, and 5.6% across the four basins, but have expression complexities 15–25 versus rank-1 complexity 8 and five fitted parameters. Rank-2 gives no meaningful gain over rank-1: it is worse in Atlantic and Indian and improves by only 0.29 and 0.04 μmol kg⁻¹ in Pacific and Southern Ocean.

**Remaining limitation:** “Efficiency” is descriptive, not a formal decision-theoretic optimum. The study does not establish that rank-1 is the universal best compression under every loss or complexity penalty.

**Manuscript location:** Results, “Rank-1 recurrence and operator-space dependence,” “Accuracy–complexity tradeoff,” and “Rank-1 versus rank-2”; Discussion, “Why rank-1 remains useful despite higher-accuracy models.”

## Reviewer 5 — “Is this ocean physics or an empirical statistical artifact?”

**Strongest criticism:** A squared linear combination discovered from data may have no mechanistic interpretation.

**Quantitative answer:** The rank-1 Hessian establishes rank-1 curvature in predictor space, and derivative signs are consistent across the observed domain. These are mathematical properties of the fitted empirical surface. The analysis does not derive the expression from carbonate chemistry, prove a mixing mechanism, or establish causality.

**Remaining limitation:** Water-mass definitions are threshold-based, AAIW and abyssal failures indicate omitted or non-transferable state information, and physical interpretation remains qualitative.

**Manuscript location:** Methods, “Derivative analysis”; Results, “Derivative structure”; Discussion, “Physical interpretation.”

## Overall reviewer-resistance assessment

The remaining vulnerabilities are transparent rather than hidden: finite operator spaces, heuristic symbolic-family classification, dependence within spatial blocks/cruises, and domain-specific transfer failure. The final manuscript explicitly treats these as limitations and does not claim a universal physical law.
