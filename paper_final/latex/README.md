# Symbolic Regression Discovery of a Compact Quadratic Structure for TCO₂

## LaTeX Project Structure

```
latex/
├── main.tex              # Main manuscript
├── supplementary.tex     # Supplementary material
├── references.bib        # Bibliography
├── Makefile              # Build automation
├── README.md             # This file
└── figures/              # All figures
    ├── final_figure1_scatter.png
    ├── final_figure2_benchmark.png
    ├── final_figure3_depth.png
    ├── final_figure4_uncertainty.png
    ├── COMPLETENESS_FIGURES.png
    └── ...
```

## Compilation

```bash
# Compile main manuscript
make main

# Compile supplementary material
make supp

# Compile everything
make all

# Clean auxiliary files
make clean
```

Or manually:
```bash
pdflatex main
bibtex main
pdflatex main
pdflatex main
```

## Key Corrections Applied

### Critical Fixes
1. **Sample size**: N = 472,154 (was incorrectly stated as 147,387)
2. **Temporal validation**: "Post-2018 chronologically held-out" (was "8-year external holdout")

### Major Fixes
3. **Statistical equivalence**: Removed; replaced with observed RMSE differences + TOST
4. **Bootstrap**: Hybrid bootstrap (case resampling + noise), n=1,000 replicates
5. **Physical interpretation**: "Qualitatively compatible" (was "derives")
6. **AAIW**: "AAIW-consistent residual structure" (was "diagnosis")
7. **Inventory**: "Basin-mean concentration" (was "basin-scale inventory")

### Moderate Fixes
8. **Hessian**: H = 2vv^T (was vv^T)
9. **Coefficient signs**: Derivative analysis added
10. **VIF**: "Consistent with" (was "confirms")
11. **Depth limitations**: Explicitly stated in abstract and conclusions

## New Analyses Added

1. **Robustness testing**: Subsample stability, random seed sensitivity, AOU scaling, complexity comparison
2. **TOST equivalence test**: Formal equivalence testing with δ = 2 and 5 μmol kg⁻¹
3. **Cluster bootstrap**: Spatial cluster resampling (5° lat × 10° lon bins)
4. **Water-mass classification**: T-S-based water-mass identification and AAIW analysis

## Status

**READY FOR SUBMISSION** (pending author/institution information)

All numerical claims are supported by corrected analysis. All placeholder text has been flagged for replacement.
