# Submission Package — Ocean $\text{TCO}_2$ Symbolic Regression Study

This directory contains the submission-ready manuscript package for:

> **How Low-Dimensional Is the S--T--AOU Dependence of Ocean Total Dissolved Inorganic Carbon? A Symbolic Regression Analysis**

---

## 1. Study Objective
This study investigates how low-dimensional the empirical dependence of ocean total dissolved inorganic carbon ($\text{TCO}_2$) on salinity ($S$), temperature ($T$), and apparent oxygen utilisation ($\text{AOU}$) is across global ocean basins. We ask: *How much predictive accuracy is retained as structural constraints become stronger?*

## 2. Main Scientific Contribution
- **Substantial Compressibility**: A five-parameter rank-1 quadratic model, $\text{TCO}_2 = (\alpha S + \beta T + \gamma \text{AOU}/100 + \delta)^2 + \epsilon$, captures substantial predictive structure across $472{,}530$ quality-controlled GLODAP v2.2023 observations.
- **One-Dimensional Quadratic Curvature**: Increasing the quadratic rank from one to two provides no meaningful predictive benefit ($\Delta\text{RMSE} = +1.98, +0.25, -0.29, -0.04\,\mu\text{mol kg}^{-1}$ across Atlantic, Indian, Pacific, Southern Ocean), demonstrating that quadratic curvature in predictor space is effectively one-dimensional.
- **Accuracy vs. Complexity Tradeoff**: Higher-complexity symbolic expressions achieve 3--12% lower RMSE but require 2--3× higher PySR expression complexity. Rank-1 is an efficient low-complexity approximation rather than an invariant Pareto knee.
- **Operator-Space Sensitivity**: Rank-1 candidates recurred in 24 of 48 search cells across four operator spaces and three seeds, but recovery depended on operator set choice (particularly in the Atlantic).
- **Domain Failure Boundaries**: Severe systematic bias in Antarctic Intermediate Water (AAIW; rank-1 RMSE 28.43, bias $+27.14\,\mu\text{mol kg}^{-1}$) and global coefficient failure in the Southern Ocean abyss ($>3000$\,m; global $R^2 = -0.174$ vs. local $R^2 = 0.786$) delineate clear limits of basin-scale empirical parameterizations.

## 3. Dataset
- **Source**: GLODAP v2.2023 ($N = 472{,}530$ QC observations).
- **Basins**: Atlantic ($N=147{,}448$), Indian ($N=40{,}391$), Pacific ($N=187{,}443$), Southern Ocean ($N=97{,}248$).
- **Partitioning**: Training ($<2015$, $N=378{,}507$), Temporal Validation ($2015$--$2017$, $N=53{,}474$), Locked External Holdout ($\ge 2018$, $N=40{,}549$).

## 4. Package Contents
- `final_manuscript.pdf`: Compiled submission PDF.
- `manuscript_source/`: Modular LaTeX source files (`main.tex`, `sections/`).
- `figures/`: High-resolution figures (`fig1_data_overview.png`, `fig2_model_comparison.png`, `fig1_pareto_frontier_final.png`, `fig4_rank1_recurrence.png`, `fig5_generalization.png`, `fig5_depth_stratified.png`, `fig7_watermass_performance.png`, `fig8_southern_abyss.png`, `fig_derivative_analysis.png`).
- `tables/`: Generated TeX table files (`table1` through `table13`).
- `supplementary/`: Supplementary material (`supplementary.tex`).
- `references.bib`: Complete BibTeX citations.
- `cover_letter_draft.md`: Submission cover letter draft.

## 5. Reproduction & Execution
To reproduce the full pipeline and compile the manuscript:

```bash
# 1. Execute analysis pipeline
python3 scripts/run_all.py

# 2. Perform numeric verification audit
python3 scripts/12_final_numeric_audit.py

# 3. Compile LaTeX manuscript
cd final_latex
tectonic main.tex -o build
```

## 6. Environment
Requirements are documented in `environment/environment.yaml`. Key libraries: `PySR` (v0.16+), `Julia` (1.9+), `scikit-learn`, `scipy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `tectonic`.

## 7. Important Limitations
- The rank-1 quadratic is an empirical predictive projection, not a derived physical law.
- Symbolic discovery depends on operator space parameterization.
- Basin-scale models should not be extrapolated to AAIW or abyssal water masses without regional recalibration.
