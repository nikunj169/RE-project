# Verification report for the revised S–T–AOU–TCO₂ manuscript

## Status

**Manuscript readiness: B — revision required before submission.** The scientific narrative is now limited to claims supported by saved results, but the analysis package retains unresolved implementation and artifact issues that should be corrected or documented before release.

## Authoritative sources

The manuscript source is `submission/manuscript_source/main.tex` and its section files. The root `submission/manuscript.tex` and existing PDF were not treated as authoritative because they contain stale paths and superseded claims. Numeric values were checked against `results/processed/`, `results/raw/symbolic_regression/`, `data/processed/`, and the executable analysis scripts.

## Dataset and sample accounting

- Raw GLODAP rows: 1,402,829.
- Rows with finite primary variables: 513,954.
- Final primary-variable QC rows: 512,384.
- Rows assigned to named basins: 472,530.
- Basin-unclassified rows: 39,854.
- Named-basin totals: Atlantic 147,448; Indian 40,391; Pacific 187,443; Southern Ocean 97,248.
- Temporal totals: training 378,507; validation 53,474; external 40,549.

The apparent 472,530 versus 512,384 discrepancy is therefore not a duplicate-observation explanation: it is the difference between all QC-passing rows and the subset assigned to the four named basins. The supplied predictor-distribution table is incomplete/stale because it places the full QC count under `Unclassified` and omits the named-basin rows. It is excluded from the revised main manuscript.

## Figure 1

The stored Figure 1 was empty because the processed CSV previously contained stale/unresolved basin labels while the plotting code filtered for the four named labels. The authoritative preprocessing step was rerun from the supplied raw GLODAP file, and the data-overview figure was regenerated. The current figure contains nonzero observations in all four panels and displays the verified basin counts above. The analysis package should still be made single-source so future figure regeneration cannot silently reuse stale labels.

## Model definitions

Rank-1 is the five-parameter model

\[
\mathrm{TCO}_2=(\alpha S+\beta T+\gamma\,\mathrm{AOU}/100+\delta)^2+\epsilon.
\]

Rank-2 is a sum of two squared affine projections plus an offset. The implementation does not impose orthogonality, normalization, or a canonical ordering. The component representation is rotationally non-identifiable and the resulting quadratic form has rank at most two. The manuscript now uses this wording and does not call the projections orthogonal.

## Symbolic-regression audit

- Completed searches: 48 = 4 basins × 4 operator spaces × 3 seeds.
- Candidate rows: 599, not deduplicated across seeds or operator spaces.
- Heuristic rank-1 classifications: 86 raw rows.
- Cells containing at least one rank-1 candidate: 24/48.
- Effective completed-run settings: 20 populations, population size 40, 50 iterations, maximum size 25, and up to 30,000 training observations per run.
- Candidate selection used 2015–2017 validation RMSE; the post-2018 holdout was locked.

The YAML file lists larger nominal settings for some search spaces. Those nominal values were not the effective settings used by the completed experiment and are not reported as the principal configuration in the manuscript.

## Verified headline results

- Rank-1 temporal-validation RMSE: 20.172, 14.483, 16.002, 10.575 μmol kg⁻¹ for Atlantic, Indian, Pacific, and Southern Ocean.
- Rank-2 minus rank-1 temporal-validation RMSE: +1.984, +0.254, −0.289, −0.035 μmol kg⁻¹.
- Rank-1 spatial-block RMSE: 16.299, 16.609, 16.031, 9.157 μmol kg⁻¹.
- Rank-1 cruise-block RMSE: 16.910, 17.394, 16.031, 9.225 μmol kg⁻¹.
- Rank-1 locked external RMSE: 17.672, 18.413, 13.821, 8.985 μmol kg⁻¹.
- Atlantic AAIW rank-1 RMSE/bias: 28.426/+27.143 μmol kg⁻¹.
- Atlantic NADW rank-1 RMSE/bias: 9.767/−2.701 μmol kg⁻¹.
- Southern abyss global/local rank-1 R²: −0.173896/0.785601; RMSE 7.447/3.182 μmol kg⁻¹.

## Claims revised or omitted

| Issue | Disposition |
|---|---|
| Rank-2 described as orthogonal | Corrected; no orthogonality constraint was implemented. |
| Rank-1 interpreted as physical one-dimensionality | Replaced with a conditional empirical-compression claim. |
| Rank-1 called a universal Pareto knee | Removed; it is described as a low-complexity point on the trade-off. |
| NADW RMSE 13.62 and bias +0.92 | Removed and replaced with the verified 9.767 and −2.701 values. |
| Surface-depth range 22.5–35.1 | Removed because it was not the verified scope of the depth table. |
| Derivative signs 100% in every basin | Not used as a headline result. The stored derivative output gives 99.9864% in Atlantic and 100% in the other basins. |
| Prediction-interval coverage 94.2–95.6% | Omitted; no saved prediction-interval output was available for audit. |
| Pacific equivalence conclusion | Omitted from the main paper and documented below because the inference implementation is not sufficient for the original claim. |

## Equivalence-analysis audit

The stored analysis has three unresolved problems. First, the code contrast and manuscript definition use opposite signs for the paired absolute-error difference. Second, the TOST rows are explicitly observation-level/naive rather than cluster-adjusted. Third, the bootstrap implementation samples cluster labels with replacement but reconstructs rows using membership rather than preserving repeated cluster multiplicities, so it is not a standard cluster bootstrap. The Pacific interval `[-0.626, -0.139]` therefore cannot support the original statement that practical equivalence is rejected. Even under the stored sign convention, the interval lies inside all tested margins of ±1, ±2, ±5, and ±10 μmol kg⁻¹; it may indicate a small directional difference, not a practically important non-equivalence result.

## Remaining release blockers

1. Replace anonymous submission metadata, repository metadata, funding, and acknowledgements as appropriate.
2. Use one asset generator and remove competing legacy/publication table numbering.
3. Correct or remove the stale predictor-distribution table from the release package.
4. Add locked environment files for both Python and Julia where possible.
5. Run the complete pipeline in a clean environment and compare generated-output checksums.
6. Inspect the compiled PDF page by page for figure readability, table overflow, and reference warnings.
7. Add automated text-to-table checks for the headline values and a corrected dependence-aware equivalence analysis if formal equivalence claims are required.
