# Final Recommendation

## Recommendation: **B — MANUSCRIPT REVISION REQUIRED**

The required spatial and cruise validation has now been executed, the external holdout remains locked, and the manuscript compiles after generated-table corrections. The package is scientifically defensible but not yet ready for submission because author/repository metadata remain intentionally anonymous/placeholders and the release still needs a clean-environment reproducibility verification.

## Completed requirements

- Spatial 5° × 10° blocked validation completed for four basins and five model types.
- Cruise-blocked validation completed using the available `cruise` identifier.
- External post-2018 rows were excluded from all model-selection and blocked-validation operations; 42,496 rows were reserved for external evaluation.
- Frozen symbolic candidates were selected using temporal-validation RMSE only.
- Final train/temporal-validation/external model hierarchy generated.
- Parametric AIC/BIC generated only for comparable fitted parametric families; frozen symbolic entries are marked non-comparable.
- Recurrence and Pareto summaries regenerated from the complete candidate library.
- Publication figures and tables regenerated from saved outputs.
- Manuscript rewritten around the evidence and blocked-validation results.
- Tectonic compilation completed after fixing generated LaTeX escaping.
- Old project remained untouched.

## Main quantitative findings

- Rank-1 candidates: 86/599; present in 24/48 basin–operator–seed cells.
- Rank-1 temporal-validation RMSE: 20.17, 14.48, 16.00, 10.58 μmol kg⁻¹ across Atlantic, Indian, Pacific, Southern Ocean.
- Rank-2 versus rank-1: worse by 1.98 and 0.25 μmol kg⁻¹ in Atlantic/Indian; better by only 0.29 and 0.04 μmol kg⁻¹ in Pacific/Southern Ocean.
- Spatial rank-1 RMSE: 16.30, 16.61, 16.03, 9.16 μmol kg⁻¹.
- Cruise rank-1 RMSE: 16.91, 17.39, 16.03, 9.22 μmol kg⁻¹.
- Locked external rank-1 RMSE: 17.67, 18.41, 13.82, 8.99 μmol kg⁻¹.
- AAIW rank-1 RMSE ≈28.43 μmol kg⁻¹, bias ≈+27.14 μmol kg⁻¹; full quadratic RMSE ≈22.33 μmol kg⁻¹.
- Southern abyss global rank-1 R²=-0.174; local R²=0.786 and RMSE≈3.18 μmol kg⁻¹.

## Remaining blockers before “READY FOR SUBMISSION”

1. Replace anonymous author/affiliation metadata with the actual submission metadata.
2. Add final repository URL/DOI and funding/acknowledgement statements.
3. Create and verify exact Python and Julia lockfiles, including Julia `Project.toml`/`Manifest.toml` where possible.
4. Run the complete final-output sequence in a clean environment and compare checksums.
5. Inspect the compiled PDF page by page for table overflow, figure readability, and reference/cross-reference warnings.
6. Add final manuscript numeric cross-check automation and release checksum manifest.

## Publication claim

The defensible claim is that TCO₂’s empirical S–T–AOU dependence is substantially compressible into a low-dimensional nonlinear predictive projection. A five-parameter rank-1 quadratic captures much of the structure of more flexible models, while rank-2 adds little. The approximation is not operator-invariant, universally optimal, universally equivalent to full quadratic, or a physical law; it has basin-, depth-, water-mass-, and validation-regime-dependent failures.
