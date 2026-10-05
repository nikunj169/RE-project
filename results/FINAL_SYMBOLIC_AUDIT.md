# Final Symbolic Regression Audit

**Generated**: 2026-09-10
**Auditor**: Automated analysis of 599 PySR candidates across 48 runs

---

## 1. Candidate Library Integrity

- **Total candidates**: 599
- **Runs**: 48 (4 basins × 4 search spaces × 3 seeds) — all successful
- **Evaluation**: All candidates evaluated on the same locked 2015–2018 validation set
- **External holdout**: Post-2018 data was NOT used for candidate selection. Only validation metrics (RMSE_val) appear in the candidate library.
- **Leakage check**: No candidate was selected using external holdout information.

## 2. Complexity Definition

- **PySR expression complexity**: NOT the same as fitted parameter count
- Complexity is a structural measure counting operator nodes in the expression tree
- The rank-1 quadratic `(αS + βT + γA + δ)² + ε` has complexity 8 in PySR but 5 fitted parameters
- This distinction is critical and must be clearly stated in the manuscript

## 3. Family Classifier Assessment

**Potential misclassification risks**:

| Risk | Impact | Mitigation |
|------|--------|------------|
| Variable names differ (sal/tmp/aou100 vs x0/x1/x2) | Fixed — classifier updated | Already addressed |
| Algebraically equivalent rank-1 forms with different syntax | May be missed | Manual inspection of high-complexity candidates recommended |
| Cross-terms that happen to be zero | Misclassified as rank-1 | Low risk — PySR wouldn't retain zero terms |
| `(a*sal + b*tmp)^2 + c` vs `(a*sal)^2 + 2ab*sal*tmp + (b*tmp)^2 + c` | Both rank-1 but second shows cross-terms | Second would be classified as full quadratic — this is a known limitation |

**Net assessment**: The family classifier is adequate but imperfect. The rank-1 recurrence count (86) may include some false positives from expressions that look rank-1 but contain hidden cross-structure, and may miss some true rank-1 expressions that PySR expanded into cross-term form.

## 4. Deduplication

- Expressions were NOT deduplicated across seeds or search spaces
- The same mathematical expression may appear multiple times
- The 86 rank-1 candidates likely include duplicates
- For the Pareto frontier, duplicates at the same complexity don't affect the frontier itself
- For recurrence counting, duplicates inflate the count

**Recommendation**: Report both raw count (86) and unique-basin-search count (number of basin×search cells with at least one rank-1 candidate)

## 5. Leakage Assessment

| Question | Answer |
|----------|--------|
| Were all candidates evaluated on the same validation data? | YES |
| Was the external holdout used for selection? | NO |
| Was the Pareto frontier computed without leakage? | YES — computed on validation metrics only |
| Were any hyperparameters tuned after observing test performance? | NO — PySR hyperparameters were fixed in config |

**Verdict**: No leakage detected.

## 6. Known Limitations

1. PySR subsampled training data to 30,000 observations per run — larger samples might find different structures
2. PySR settings (populations=20, niterations=50) are moderate — more intensive search might find better candidates
3. The deterministic mode ensures reproducibility but may explore less of the search space than stochastic mode
4. The family classifier is heuristic, not algebraic — a proper symbolic equivalence checker would be more reliable
