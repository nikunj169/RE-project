# Final Manuscript Audit Report

**Date**: 2026-09-14  
**Pipeline Status**: Complete & Verified (100% Numeric Match)  
**Manuscript Compilation**: Success (`final_latex/build/main.pdf`, 2.35 MiB)

---

## 1. Scientific Story Summary
The manuscript examines how low-dimensional the empirical dependence of ocean total dissolved inorganic carbon ($\text{TCO}_2$) on salinity ($S$), temperature ($T$), and apparent oxygen utilisation ($\text{AOU}$) is. Analyzing $472{,}530$ quality-controlled GLODAP v2.2023 observations across four basins, the paper evaluates a 599-expression symbolic regression Pareto library against a systematic model hierarchy (mean baseline, linear, rank-1 quadratic, rank-2 quadratic, full quadratic, and cubic models).

## 2. Final Central Claim
"The empirical $S$--$T$--$\text{AOU}$ dependence of ocean $\text{TCO}_2$ is substantially compressible into a low-dimensional nonlinear representation. A five-parameter rank-1 quadratic captures much of the predictive structure of more flexible models, while increasing the quadratic rank from one to two provides little additional predictive benefit. However, the compact structure is not universally invariant: symbolic recovery depends on operator space, predictive performance varies by basin, and specific water-mass and abyssal regimes exhibit substantial domain dependence."

## 3. Evidence Supporting Central Claim
- **Compressibility**: Rank-1 quadratic ($K=5, C=8$) outperforms linear models in all basins, achieving validation RMSEs of 20.17 (Atl), 14.48 (Ind), 16.00 (Pac), 10.58 (SO) $\mu\text{mol kg}^{-1}$.
- **One-Dimensional Quadratic Rank**: Rank-2 ($K=9$) fails to improve upon rank-1 ($\Delta\text{RMSE}$ = $+1.98$, $+0.25$, $-0.29$, $-0.04\,\mu\text{mol kg}^{-1}$), demonstrating that quadratic curvature in state space is effectively one-dimensional.
- **Accuracy vs. Complexity**: Higher-complexity symbolic candidates reduce RMSE by 3--12% but require 2--3× higher expression complexity ($C=15$--$25$ vs. $C=8$).
- **Recurrence & Operator Sensitivity**: Rank-1 candidates recurred in 24/48 search cells across all four basins. In the Atlantic, rank-1 recovery required explicit square operators.
- **Domain Failure Regimes**: Severe systematic bias in Antarctic Intermediate Water (AAIW; rank-1 RMSE 28.43, bias $+27.14\,\mu\text{mol kg}^{-1}$) and global coefficient failure in the Southern Ocean abyss ($>3000$\,m; global $R^2 = -0.174$ vs. local $R^2 = 0.786$) demonstrate clear domain transferability boundaries.

## 4. Claim Classification Table

| Claim | Evidence | Status |
|---|---|---|
| Rank-1 quadratic captures substantial predictive structure | Outperforms linear in all 4 basins; within 0.1-0.75 $\mu\text{mol kg}^{-1}$ of full quad outside Atlantic | **SUPPORTED** |
| Adding rank-2 quadratic rank provides little benefit | Rank-2 worse in Atlantic (+1.98) & Indian (+0.25); marginal in Pacific (-0.29) & SO (-0.04) | **SUPPORTED** |
| Rank-1 recurs across seeds and basins | Present in 24/48 cells; present in all 4 basins | **SUPPORTED** |
| Rank-1 recurrence is operator-space invariant | Atlantic rank-1 recovered ONLY in Search C (with explicit square) | **NOT SUPPORTED** |
| Rank-1 is the Pareto knee across all tolerance thresholds | Rank-1 meets 5% tolerance threshold ONLY in Pacific; fails 1%, 2%, 5% in Atlantic, Indian, SO | **NOT SUPPORTED** |
| Temporal generalization holds | Post-2018 external holdout RMSE: 17.67 (Atl), 18.41 (Ind), 13.82 (Pac), 8.99 (SO) | **SUPPORTED** |
| Spatial and cruise generalization holds | Blocked cross-validation RMSE remains low, but model rankings shift between schemes | **PARTIALLY SUPPORTED** |
| Dependence-aware equivalence holds universally | Pacific paired-error CI is [-0.63, -0.14] $\mu\text{mol kg}^{-1}$ (entirely negative) | **NOT SUPPORTED** |
| AAIW is a major failure domain | Rank-1 RMSE 28.43 $\mu\text{mol kg}^{-1}$, bias $+27.14\,\mu\text{mol kg}^{-1}$ | **SUPPORTED** |
| Southern Ocean abyss is a global transfer failure | Global fit $R^2 = -0.174$; local fit $R^2 = 0.786$ | **SUPPORTED** |
| Rank-1 is a fundamental physical ocean-carbon law | Empirical fit; Hessian $H=2vv^T$ is mathematical; derivative signs describe empirical surface | **NOT SUPPORTED** |

## 5. Major Limitations
1. Predefined operator spaces and complexity bounds ($C_{\max}=30$).
2. Empirical surface properties do not prove physical causality or carbonate equilibrium derivation.
3. Basin-scale models cannot be extrapolated to AAIW or abyssal regimes without regional recalibration.

## 6. Numerical & Figure/Table Consistency Audit
- Total sample N: $472{,}530$ (Train $378{,}507$, Validation $53{,}474$, External $40{,}549$) — **PASS**
- Basin sample Ns: Atlantic $147{,}448$, Indian $40{,}391$, Pacific $187{,}443$, Southern Ocean $97{,}248$ — **PASS**
- Search cells: 48 runs ($4 \times 4 \times 3$), 599 candidates, 252 Pareto, 86 rank-1, 24/48 cells — **PASS**
- Full quadratic parameter count: 10 parameters (not 9) — **PASS**
- Rank-1 parameter count vs complexity: 5 fitted coefficients, PySR complexity 8 — **PASS**
- All 13 TeX tables generated and compiled cleanly — **PASS**
- All 9 figures cross-referenced and rendered in PDF — **PASS**

## 7. Remaining Reviewer Vulnerabilities & Mitigations
- **Vulnerability**: "Is symbolic regression claiming a new physical law?"  
  *Mitigation*: Section 1, Section 5.7, and Section 6 explicitly state rank-1 is an empirical projection and reject physical law claims.
- **Vulnerability**: "Why didn't rank-1 recur in all search cells?"  
  *Mitigation*: Section 4.3 and Section 5.4 document operator-space sensitivity (Search C requirement in Atlantic).
- **Vulnerability**: "Is rank-1 equivalent to full quadratic everywhere?"  
  *Mitigation*: Section 4.15 reports Pacific paired-error CI [-0.63, -0.14] $\mu\text{mol kg}^{-1}$ and rejects universal equivalence.
- **Vulnerability**: "Does the model fail in deep waters?"  
  *Mitigation*: Section 4.12 & 4.13 explicitly detail AAIW bias and Southern Ocean abyssal transfer failure.
