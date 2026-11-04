# Figure Captions

Figure files are in `research/figures/` (PNG, 300 dpi, and PDF), produced by `python -m research.visualization.figures` from saved results only.

**Figure 1 (`fig1_study_area`).** Study data. (a) Positions of profiles that enter the trend analysis: strict-QC GLODAPv2.2023 samples in 4°×4° cells whose observations span ≥10 years, coloured by basin. (b) Number of these samples per year and basin (thousands), showing the WOCE peak in the 1990s and the GO-SHIP repeat-hydrography era. *Purpose: shows the spatial and temporal sampling that the within-cell design exploits.*

**Figure 2 (`fig2_trend_profiles`).** Rates of change at fixed potential temperature, salinity and AOU (µmol kg⁻¹ yr⁻¹) by neutral-density layer (plotted at the layer midpoint; the lightest layer is drawn at 25.6 and the densest at 28.35 kg m⁻³ for display). Black circles: observed DIC (β_DIC). Coloured squares: TRACE C_ant fitted with the identical estimator at the same samples (β_Cant). Error bars: 95 % cruise-cluster bootstrap intervals (B = 500). *Purpose: shows the vertical ventilation structure of carbon accumulation and its first-order agreement with tracer-based C_ant.*

**Figure 3 (`fig3_ratio_h1`).** Pre-registered test H1. Ratio R = β_DIC/β_Cant (log scale) with 95 % bootstrap intervals for layers with γⁿ ≤ 27.8 kg m⁻³. Grey band: TRACE ±15 % structural uncertainty. Circles: supported (interval overlaps the band). Crosses: contradicted. Basins are offset vertically for readability.

**Figure 4 (`fig4_delta_robustness`).** Difference δ = β_DIC − β_Cant (µmol kg⁻¹ yr⁻¹). Black: primary estimate with 95 % interval. Coloured bars: range of point estimates across the 11 TRACE-based sensitivity variants (cell size, span, conditioning form, conditioning variables, post-1990, seasonal terms, phosphate instead of AOU, unconditioned). Grey band: offset-equivalent detection limit (4 µmol kg⁻¹ divided by the median within-cell time span). *Purpose: separates robust from fragile deviations and shows their size relative to measurement-drift limits.*

**Figure 5 (`fig5_h3_holdout_bias`).** Static empirical DIC models (linear, rank-1 quadratic, rank-2 quadratic, full quadratic, cubic in S, T, AOU), trained on years < 2015 and evaluated on the external set (≥ 2018), strict QC. Top: mean bias (observed − predicted, µmol kg⁻¹). Bottom: RMSE (µmol kg⁻¹). Bars compare static models, models with a linear time term, and models with TRACE C_ant as a covariate. *Purpose: quantifies how the anthropogenic trend confounds temporal validation (H3).*

**Figure 6 (`fig6_emlr_comparison`).** Independent comparison over the eMLR(C*) period. Black: observed within-cell DIC trend using samples from 1989–2012 (95 % interval). Coloured squares: TRACE C_ant change 1994–2007 divided by 13 years, averaged over the same samples. Grey diamonds: eMLR(C*) ΔC_ant 1994–2007 of Gruber et al. (2019) divided by 13 years at the same samples (standard case; bars show the 14-member ensemble range). Units µmol kg⁻¹ yr⁻¹.

**Figure S1 (`figS1_synthetic_recovery`).** Validation of the estimator on synthetic data built from four representative basin × layer groups with imposed trends and cruise offsets (100 replicates per setting). (a) Mean bias (grey band ±0.02 µmol kg⁻¹ yr⁻¹, the pre-registered criterion). (b) Coverage of nominal 95 % bootstrap intervals (grey band 90–99 %).
