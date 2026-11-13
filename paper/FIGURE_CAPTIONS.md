# Figure Captions

Figure files are in `research/figures/` (PNG at 300 dpi, and PDF). They are produced from saved results only, by `python -m research.visualization.figures` (Figs. 2–6, S1) and `python -m research.visualization.figures_round2` (Figs. 1, 7–11).

**Figure 1 (`fig1_study_area`).** Study data.
(a) Positions of the profiles that enter the trend analysis: strict-QC GLODAPv2.2023 samples in 4°×4° cells whose observations span ≥ 10 years, coloured by basin. Robinson projection centred at 150°W.
(b) Number of these samples per year and basin (thousands).

**Figure 2 (`fig2_trend_profiles`).** Rates of change at fixed potential temperature, salinity and AOU (µmol kg⁻¹ yr⁻¹) by neutral-density layer. Layers are plotted at their midpoints; the lightest is drawn at 25.6 and the densest at 28.35 kg m⁻³ for display.
- Black circles: observed DIC (β_DIC).
- Coloured squares: TRACE C_ant fitted with the identical estimator at the same samples (β_Cant).
- Error bars: 95 % cruise-cluster bootstrap intervals (B = 500).

**Figure 3 (`fig3_ratio_h1`).** Pre-registered test H1: the ratio R = β_DIC/β_Cant (log scale) with 95 % intervals, for layers with γⁿ ≤ 27.8.
- Grey band: TRACE's stated ±15 % uncertainty.
- Circles: H1 supported. Crosses: H1 contradicted.

**Figure 4 (`fig4_delta_robustness`).** δ = β_DIC − β_Cant from the single-layer estimator.
- Black: primary estimate with 95 % interval.
- Coloured bars: range across 11 analysis variants.
- Grey band: Round-1 offset-equivalent detection limit. This limit is superseded by the drift-immune estimator (Fig. 7).

**Figure 5 (`fig5_h3_holdout_bias`).** Static empirical DIC models trained on years < 2015 and evaluated on years ≥ 2018 (strict QC).
- Top: mean bias (µmol kg⁻¹).
- Bottom: RMSE (µmol kg⁻¹).
- Bars compare static models, models with a linear time term, and models with a TRACE C_ant covariate.

**Figure 6 (`fig6_emlr_comparison`).** Round-1 comparison with Gruber et al. (2019) over 1989–2012. Rates in µmol kg⁻¹ yr⁻¹.
- Observed within-cell DIC trend (black).
- TRACE 1994–2007 rate (coloured squares).
- eMLR(C*) 1994–2007 rate (grey diamonds; bars show the ensemble range).

**Figure 7 (`fig7_drift_immune_excess`).** Excess of DIC over TRACE (µmol kg⁻¹ yr⁻¹) by layer.
- Grey circles: single-layer δ (E02), sensitive to cruise-level drift.
- Coloured squares: Δδ from the differential estimator with cruise fixed effects (E06). It is relative to the deep reference layers (grey shading, γⁿ > 27.8) and immune to any additive per-cruise offset.
- Error bars: 95 % cruise-cluster bootstrap intervals (B = 300).

*Purpose: shows that the upper-ocean excess survives control of calibration drift, while the depth-uniform deep deviations (Atlantic −, Pacific +) are removed.*

**Figure 8 (`fig8_scale_factor`).** Drift-immune excess Δδ versus the TRACE differential rate for the 20 upper and intermediate basin × layer units, with 95 % intervals.
- Line and shading: least-squares proportionality factor m with its global-bootstrap 95 % interval (B = 500).
- Legend: per-basin m.

This proportional summary is post hoc. Points above the line (intermediate layers outside the Southern Ocean) carry relatively more excess.

**Figure 9 (`fig9_trace_sensitivity`).** δ for TRACE structural variants, from the single-layer estimator with DIC unchanged.
- Default with 95 % interval (black).
- Transit-time-distribution shape Δ/Γ = 1.0 and 1.6.
- Disequilibrium scale ×0.9 and ×1.1.
- Mean ages ×0.8 and ×1.2.
- Potential temperature as input.
- Unadjusted atmospheric CO₂ boundary condition (orange stars; E12, post hoc).

**Figure 10 (`fig10_emlr_products`).** Ratio of C_ant accumulation rate to TRACE at the same samples and periods, for layers where the TRACE rate exceeds 0.05 µmol kg⁻¹ yr⁻¹.
- Gruber et al. (2019), 1994–2007.
- Müller et al. (2023), 1994–2004 and 2004–2014 (standard cases).
- Observed within-cell DIC (stars).
- Grey band: ±15 %.

No regression is involved for the eMLR ratios.

**Figure 11 (`fig11_cstar_controls`).** Sensitivity to the remineralisation correction, on samples with quality-controlled PO₄ and TAlk.
- Excess of DIC over TRACE conditioned on θ, S and AOU (black).
- Excess of the pump-corrected tracer C* conditioned on θ and S (coloured).
- −½ × the TAlk trend at fixed θ, S (crosses): the alkalinity contribution to C*. Its negative deep values indicate TAlk measurement drift.

**Figure S1 (`figS1_synthetic_recovery`).** Validation of the single-layer estimator on synthetic data: (a) mean bias; (b) coverage of nominal 95 % bootstrap intervals.
