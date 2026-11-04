# Phase 2–3: GLODAP Science Notes and Reconstruction of the Existing Study

## Part A — GLODAPv2.2023 as used here (Phase 2)

**Product.** GLODAPv2.2023 (Lauvset et al., 2024, ESSD 16:2047, doi:10.5194/essd-16-2047-2024) is a synthesis of ship-based interior-ocean bottle data (1,402,829 samples in the merged master file, 1972–2021 in our QC subset). It is a *bottle* product (discrete samples on cruise/station/cast/bottle), distinct from the *mapped* GLODAPv2 climatology (Lauvset et al., 2016), which normalised DIC to 2002 using TTD-based anthropogenic-carbon rates.

**File acquired.** `data/raw/GLODAPv2.2023_Merged_Master_File.mat` from NOAA NCEI accession 0283442; SHA-256 `002881fa…8766` (see `data/raw/SHA256SUMS`); read-only. The repository’s processed CSV is reproduced exactly from it (512,384 rows; identical values).

**Quality control layers.**
* *Primary QC* — WOCE flags per sample (`G2xxxf`): 2 = good, 0 = value calculated/interpolated (for DIC: computed from other carbonate parameters), 9 = missing.
* *Secondary QC* — cruise-level flags (`G2xxxqc`): 1 = the cruise’s data for that variable passed crossover-based secondary QC (adjustments applied for systematic offsets, judged on deep data); 0 = not (yet) subjected / not passed. GLODAP recommends qc = 1 for most uses.
* The repository used **neither**: 35,819 (7.0 %) of its rows are calculated DIC, 67,839 (13.2 %) lack secondary QC. Strict QC (DIC f=2 & qc=1, S and O2 f=2 & qc=1) retains 428,264 rows.
* Stated internal-consistency targets for adjustments are ~4 µmol kg⁻¹ for DIC (Olsen et al., 2016; Lauvset et al., 2024). **Implication:** a trend of 0.1 µmol kg⁻¹ yr⁻¹ sustained over 40 years is a 4 µmol kg⁻¹ change — at the level of residual inter-cruise offsets. Deep-ocean trends are therefore at the edge of detectability, and this must be treated explicitly (cruise-level random effects / cruise-cluster inference, leave-one-cruise-out).

**Variables available (finite counts in raw file).** DIC 548,917; TAlk 496,638; O2 1,252,778; NO3 1,123,248; PO4 1,070,157; SiO4 1,122,518; θ 1,363,780; pressure 1,394,969; neutral density γⁿ; CFC-12 400,054; SF6 103,791; δ¹³C 37,511; Δ¹⁴C 42,331; pH 492,023. Units: µmol kg⁻¹ for chemistry, °C (ITS-90) temperature, PSS-78 salinity, dbar pressure.

**Coverage and sampling bias.** Sampling is dominated by WOCE (1990–1998) and GO-SHIP repeat sections (2003–present); pre-1990 data are sparse (< 3,500 yr⁻¹); 2019–2021 are sparse. Coverage is boreal-summer-biased at high latitudes and section-based (dense along lines, absent between). Time and space are therefore strongly confounded: a naive regression on year partly measures *where* cruises went in each era. Any trend estimate must identify time effects *within* locations.

**Key oceanographic relationships for this study.**
* DIC = preformed DIC + regenerated (soft tissue ∝ AOU via Redfield-like ratio; carbonate pump) + anthropogenic CO2. Preformed DIC depends on outcrop θ, S, alkalinity and air–sea disequilibrium.
* AOU = O2,sat(θ,S) − O2. AOU is **not** independent of θ and S, and AOU ≠ true oxygen utilisation where water forms undersaturated (Southern Ocean deep-water formation: disequilibrium up to 50–100 µmol kg⁻¹; Ito et al., 2004). Conditioning DIC on AOU removes most respiration signal but with a known Southern-Ocean bias.
* Anthropogenic CO2 enters at the surface and is carried along ventilation pathways; under transient steady state (TSS), the local accumulation rate is proportional to the local Cant inventory (Tanhua et al., 2007; Gruber et al., 2019). Expected vertical structure: ~1 µmol kg⁻¹ yr⁻¹ in surface/mode waters, decaying with density; larger in ventilated North Atlantic deep water than in old Pacific/Indian deep water.
* Correlation ≠ mechanism: a residual DIC trend at fixed (θ, S, AOU) is consistent with, but not proof of, Cant invasion; it also contains changes in preformed DIC (disequilibrium, alkalinity), CaCO3 dissolution, AOU-disequilibrium changes, and residual measurement offsets.

## Part B — Reconstruction of the existing research (Phase 3)

| # | Question | Answer |
|---|---|---|
| 1 | Central question | How low-dimensional is TCO2(S, T, AOU)? How much accuracy survives rank-1 / rank-2 / full-quadratic constraints? |
| 2 | Hypotheses | H1 rank-1 ≈ full quadratic; H2 rank-2 adds nothing (curvature is 1-D); H3 PySR repeatedly rediscovers rank-1; H4 (implicit) temporal stability |
| 3 | Methodology | Basin-wise fits on pre-2015 data; 2015–17 validation; ≥2018 external; PySR (48 runs); blocked CV; TOST; threshold water masses |
| 4 | Assumptions | Stationarity in time; exchangeable observations; in-situ T adequate; AOU = respiration; no QC flags needed |
| 5 | Results | Rank-1 val RMSE 20.2/14.5/16.0/10.6; rank-2 ΔRMSE +2.0…−0.04; SR 3–12 % better at 2–3× complexity; failure in SO abyss and Atlantic AAIW |
| 6 | What they establish | Within basins, DIC is very well predicted (R² 0.92–0.99) by smooth functions of S, T, AOU; the choice among smooth low-order forms changes RMSE by ≲ 2 µmol kg⁻¹ — below DIC measurement consistency (~4). |
| 7 | What they fail to establish | (a) that curvature is “one-dimensional” — rank-1 is nested in rank-2/full-quad, which have *lower* training error; out-of-sample reversals arise from distribution shift. (b) temporal stability — residuals drift ≈ +0.3–0.4 µmol kg⁻¹ yr⁻¹ (depth-averaged), with the vertical structure of anthropogenic CO2. (c) equivalence — the TOST is anticonservative and the cluster bootstrap is implemented incorrectly. (d) “rank-1 recurrence” — family labels come from a string heuristic. |
| 8 | Strongest limitations | Unmodelled time trend confounds every temporal-validation comparison; differences between model families are below measurement uncertainty; 15 % of data are calculated or non-QC’d DIC. |
| 9 | Already known in literature | Empirical DIC/TAlk regressions on hydrographic predictors are mature (LIR, ESPER, CANYON-B, eMLR). Modern tools include time or a Cant adjustment explicitly (CANYON-B, ESPER; the mapped GLODAPv2 normalises to 2002). That the static S–T–AOU form drifts in time is expected from Cant uptake. |
| 10 | Most promising new direction | Turn the confounder into the object of study: a sample-matched, within-location test of whether the observed DIC change at fixed (θ, S, AOU) agrees with independent tracer-based (TTD/TRACE) Cant accumulation, by basin and density layer, over 1972–2021 — and quantify how ignoring it corrupts temporal validation of empirical DIC models. |

**Reproduction status.** Rank-1 results reproduced exactly with independent code; full-quadratic within 0.2 µmol kg⁻¹ (conditioning). PySR library not re-run (cost; and not needed for the new direction). Repository tests not executed (sandbox policy blocked running repository code); independent re-implementation used instead.
