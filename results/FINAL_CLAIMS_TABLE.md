# Final Claims Table

**Status vocabulary**: SUPPORTED / PARTIALLY SUPPORTED / NOT SUPPORTED / INCONCLUSIVE

| Claim | Evidence | Quantitative result | Status | Appropriate manuscript wording |
|---|---|---|---|---|
| Rank-1 candidates recur | Complete 599-candidate library across 48 runs | 86 rank-1 candidates; 24/48 basin–operator–seed cells | PARTIALLY SUPPORTED | “Rank-1 candidates recurred across all four basins and 24 of 48 search cells.” |
| Rank-1 is invariant to operator space | Four operator spaces were tested | Rank-1 occurred in all spaces in aggregate, but only 1/4 Atlantic spaces and 2/4 Pacific spaces | NOT SUPPORTED | “Recurrence was sensitive to operator space, especially in the Atlantic.” |
| Rank-1 captures substantial predictive structure | Rank-1 versus linear/full quadratic/cubic hierarchy | Rank-1 improved on linear in all basins; full quadratic was better by 0.75, 0.50, 0.11 μmol kg⁻¹ in Indian, Pacific, Southern Ocean | PARTIALLY SUPPORTED | “Rank-1 captures substantial, but not all, predictive structure.” |
| Rank-2 provides little benefit | Same train/validation/external splits | Rank-2 worse than rank-1 by 1.98 and 0.25 μmol kg⁻¹ in Atlantic/Indian; better by only 0.29 and 0.04 in Pacific/SO | SUPPORTED | “A second quadratic direction provided no meaningful additional predictive benefit.” |
| Rank-1 has a strong low-complexity tradeoff | Pareto and model hierarchy | Five fitted parameters, complexity 8; best symbolic RMSE gains 3–12% at complexity 15–25 | SUPPORTED | “Rank-1 is an efficient low-complexity approximation, not the most accurate model.” |
| Rank-1 is the Pareto knee | Pareto threshold analysis | Rank-1 within 5% of best only in Pacific; not within 5% in Atlantic, Indian, SO | NOT SUPPORTED | Do not call rank-1 the Pareto knee. |
| Temporal generalization | Locked post-2018 evaluation | Rank-1 external RMSE 17.67, 18.41, 13.82, 8.99 μmol kg⁻¹ | SUPPORTED | “Temporal external performance remained strong, with basin-specific differences.” |
| Spatial generalization | 5°×10° spatial GroupKFold | Spatial rank-1 RMSE: 16.30, 16.61, 16.03, 9.16 μmol kg⁻¹ | PARTIALLY SUPPORTED | “Performance generalized to spatial blocks, but rankings changed across regimes.” |
| Cruise generalization | Cruise GroupKFold | Cruise rank-1 RMSE: 16.91, 17.39, 16.03, 9.22 μmol kg⁻¹ | PARTIALLY SUPPORTED | “Cruise-blocked results support transfer with basin-specific degradation.” |
| External holdout remained locked | Programmatic year cutoff and frozen candidate manifest | 42,496 post-2018 rows excluded from selection/blocked validation | SUPPORTED | “Post-2018 observations were reserved for final evaluation.” |
| Basin invariance | Four basins analyzed | Rank-1 present in all basins, but recurrence and performance varied | PARTIALLY SUPPORTED | “The structure recurred across basins but was not basin-invariant.” |
| Water-mass invariance | Atlantic water-mass analysis | AAIW rank-1 RMSE ≈28.43 and bias ≈+27.14; full quadratic ≈22.33 | NOT SUPPORTED | “Performance was water-mass dependent, with AAIW as a failure regime.” |
| Southern abyss is universally predictable | Global versus local abyssal fit | Global rank-1 R²=-0.174; local R²=0.786 | NOT SUPPORTED | “Global coefficient transfer fails in the Southern Ocean abyss.” |
| Rank-1 is a physical law | Empirical fit, Hessian, derivatives | No mechanistic derivation or causal identification | NOT SUPPORTED | “Rank-1 curvature is an empirical predictive structure.” |
| Derivative signs are observed-domain properties | Derivative audit | Expected signs and positive core in 100% of observations in all basins | SUPPORTED | “The fitted surface has expected derivative signs over the observed domain.” |
| Derivative signs establish causality | Derivative audit | No causal design | NOT SUPPORTED | Do not make this claim. |
| Rank-1 and full quadratic are equivalent everywhere | Dependence-aware bootstrap | Pacific paired-error CI approximately [-0.63,-0.14] μmol kg⁻¹ | NOT SUPPORTED | “Equivalence was basin- and margin-specific; Pacific showed a small detectable disadvantage.” |

## Final scientific position

The evidence supports substantial empirical low-dimensionality and a useful five-parameter rank-1 approximation. It does not support operator-space invariance, universal accuracy, a Pareto-knee claim, universal equivalence, or physical-law interpretation.
