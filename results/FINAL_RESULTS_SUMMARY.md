# Final Results Summary

## Study Overview

**Central question**: How low-dimensional is the empirical S–T–AOU dependence of ocean TCO2?

**Data**: GLODAP v2.2023, 472,530 QC-passing observations across 4 ocean basins
**Tools**: Julia 1.12.7, PySR 2.3.0, Python 3.14.3
**Experiments**: 48 PySR runs (4 basins × 4 operator sets × 3 seeds), 599 candidates total

---

## What Was Confirmed

1. **The rank-1 quadratic structure is real**: Found in all 4 basins, all 4 search spaces (in aggregate), all 3 seeds. 86 rank-1 candidates out of 599 total.

2. **Rank-2 does NOT improve over rank-1**: In Atlantic and Indian, rank-2 is actually WORSE. In Pacific and Southern Ocean, improvement is negligible (0.04–0.29 μmol/kg). This is the strongest evidence for low-dimensionality.

3. **The model generalizes temporally**: External holdout R² ranges 0.921–0.990. No temporal degradation.

4. **Derivative signs are physically consistent**: dTCO2/dS > 0, dTCO2/dT < 0, dTCO2/dAOU > 0 in 100% of observations across all basins.

5. **Domain dependence is real**: Southern Ocean abyss R² = -0.174 (global), R² = 0.786 (local). AAIW bias = +27.1 μmol/kg. These failures are genuine and informative.

6. **Dependence-aware equivalence**: Cluster bootstrap CIs within ±2 μmol/kg for Atlantic and Southern Ocean.

## What Was NOT Confirmed

1. **Rank-1 is NOT Pareto-optimal at the knee**: Best SR candidates improve RMSE by 3–12% at 2–3× complexity. Rank-1 is at the low-complexity end of the frontier, not the knee.

2. **Rank-1 recurrence is 50% of cells, not 100%**: Atlantic only shows rank-1 in SearchC (explicit squaring operator). Pacific only in SearchC/SearchD.

3. **Full quadratic does NOT universally underperform rank-1**: In Indian, Pacific, and Southern Ocean, full quadratic is slightly better (by 0.12–0.75 μmol/kg).

4. **Equivalence does NOT hold universally**: Pacific cluster bootstrap CI [-0.63, -0.14] is entirely negative — rank-1 is statistically worse than full quadratic (though by a trivial amount).

## Key Numbers

### Model Family Comparison (Validation RMSE, μmol/kg)

| Model | Params | Atlantic | Indian | Pacific | S. Ocean |
|-------|--------|----------|--------|---------|----------|
| Mean baseline | 1 | 66.88 | 122.47 | 131.40 | 59.19 |
| Linear | 4 | 21.59 | 14.96 | 17.77 | 10.66 |
| **Rank-1 quadratic** | **5** | **20.17** | **14.48** | **16.00** | **10.58** |
| Rank-2 quadratic | 9 | 22.16 | 14.74 | 15.71 | 10.54 |
| Full quadratic | 10 | 23.41 | 13.73 | 15.50 | 10.46 |
| Cubic | 20 | 22.60 | 11.83 | 14.24 | 10.04 |
| Best SR | 15–25 | 19.00 | 12.80 | 15.49 | 9.99 |

### Rank-1 vs Rank-2 (the key comparison)

| Basin | Rank-1 | Rank-2 | Δ | Direction |
|-------|--------|--------|---|-----------|
| Atlantic | 20.17 | 22.16 | -1.99 | R1 BETTER |
| Indian | 14.48 | 14.74 | -0.26 | R1 BETTER |
| Pacific | 16.00 | 15.71 | +0.29 | R2 marginally better |
| Southern Ocean | 10.58 | 10.54 | +0.04 | R2 marginally better |

### External Holdout (post-2018, locked)

| Basin | R² | RMSE |
|-------|-----|------|
| Atlantic | 0.921 | 17.67 |
| Indian | 0.968 | 18.41 |
| Pacific | 0.990 | 13.82 |
| Southern Ocean | 0.965 | 8.99 |

### Rank-1 Recurrence

| Metric | Value |
|--------|-------|
| Total rank-1 candidates | 86 / 599 |
| Basin×search×seed cells with rank-1 | 24 / 48 (50%) |
| Basins with rank-1 | 4 / 4 (100%) |
| Search spaces with rank-1 | 4 / 4 (100%) |
| Pareto-optimal rank-1 entries | 45 / 252 (18%) |

---

## Final Assessment

**The central hypothesis is supported with important qualifications:**

The S–T–AOU relationship IS substantially low-dimensional. Rank-1 quadratic (5 params) performs comparably to full quadratic (10 params) and rank-2 does not improve. But rank-1 is not universally best — it's an efficient low-complexity approximation, not the optimal model. The 3–12% RMSE gap from best SR at 2–3× complexity is the quantitative accuracy–complexity tradeoff.

**Final recommendation: B — Manuscript revision required**

Core findings are sound. Manuscript needs to accurately represent:
1. Rank-1's Pareto position (efficient, not at knee)
2. Operator-space dependence (50% cell recurrence)
3. Pacific equivalence result (rank-1 slightly worse)
4. Add spatial-blocked validation
5. Add AIC/BIC comparison

---

## Audit Reports

- `FINAL_CLAIMS_TABLE.md` — All claims with evidence and status
- `FINAL_SYMBOLIC_AUDIT.md` — PySR candidate library audit
- `MANUSCRIPT_NUMERIC_AUDIT.md` — Cross-checked numbers
- `FINAL_RECOMMENDATION.md` — Detailed revision requirements
- `SCIENTIFIC_DECISION_REPORT.md` — Pre-audit decision report
