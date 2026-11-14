# Literature Review and Research-Gap Analysis

**Method.** Targeted web searches (2026-11-04) of publisher pages, abstracts and data-centre pages; full text was read only where open access allowed it. Each reference carries a verification tag:
**[V]** bibliographic details and the cited claim checked against publisher/abstract/data-centre page this session;
**[A]** checked only at abstract/secondary-source level (claim should be confirmed against full text);
**[M]** cited from prior knowledge, *not* re-verified this session — must be checked before submission.

---

## 1. Data product

| Ref | Content relevant here | Tag |
|---|---|---|
| Lauvset, S. K., et al. (2024). The annual update GLODAPv2.2023: the global interior ocean biogeochemical data product. *ESSD* 16, 2047–2072. doi:10.5194/essd-16-2047-2024 | Primary WOCE flags + cruise-level secondary-QC flags (1 = passed secondary QC) for 13 core variables; secondary QC targets systematic offsets using deep data. | [V] |
| Olsen, A., et al. (2016). The Global Ocean Data Analysis Project version 2 (GLODAPv2). *ESSD* 8, 297–323. doi:10.5194/essd-8-297-2016 | Original v2 synthesis; adjustment thresholds (DIC ~4 µmol kg⁻¹). | [M] |
| Lauvset, S. K., et al. (2016). A new global interior ocean mapped climatology: the 1°×1° GLODAP version 2. *ESSD* 8, 325–340. doi:10.5194/essd-8-325-2016 | Mapped product normalised DIC/pH to 2002 using TTD-based Cant rates — i.e., the community already treats the anthropogenic trend as a confounder for static mapping. | [A] |

## 2. Anthropogenic carbon (Cant) in the ocean interior

| Ref | Method / data | Finding / limitation | Tag |
|---|---|---|---|
| Sabine, C. L., et al. (2004). The oceanic sink for anthropogenic CO2. *Science* 305, 367–371. doi:10.1126/science.1097403 | ΔC* on WOCE/JGOFS | 118 ± 19 Pg C by 1994 | [M] |
| Friis, K., Körtzinger, A., Pätsch, J., Wallace, D. W. R. (2005). On the temporal increase of anthropogenic CO2 in the subpolar North Atlantic. *Deep-Sea Res. I* 52, 681–698. doi:10.1016/j.dsr.2004.11.017 | introduces eMLR | regression-difference between occupations | [M] |
| Tanhua, T., Körtzinger, A., Friis, K., Waugh, D. W., Wallace, D. W. R. (2007). An estimate of anthropogenic CO2 inventory from decadal changes in oceanic carbon content. *PNAS* 104, 3037–3042. doi:10.1073/pnas.0606574104 | TSS framework | Cant change over a period ∝ Cant at its start under invariant circulation | [A] (PMC1802019) |
| Clement, D., & Gruber, N. (2018). The eMLR(C*) method to determine decadal changes in the global ocean storage of anthropogenic CO2. *GBC* 32, 654–679. doi:10.1002/2017GB005819 | eMLR on C*; synthetic-data tests | ~±10 % (synthetic), ~±20 % expected for real data; dominant uncertainty is structural (configuration choices) | [V] |
| Gruber, N., et al. (2019). The oceanic sink for anthropogenic CO2 from 1994 to 2007. *Science* 363, 1193–1199. doi:10.1126/science.aau5153 | eMLR(C*) on GLODAPv2 | +34 ± 4 Pg C (1994–2007); globally consistent with TSS but North Atlantic ~20 % below TSS expectation | [V] |
| Müller, J. D., et al. (2023). Decadal trends in the oceanic storage of anthropogenic carbon from 1994 to 2014. *AGU Advances* 4, e2023AV000875. doi:10.1029/2023AV000875 | eMLR(C*), two decades | 29 ± 3 vs 27 ± 3 Pg C dec⁻¹; sensitivity to atmospheric CO2 fell ~15 ± 11 %; storage vs flux mismatch → natural carbon loss; authors caution limited ability to resolve decadal variability | [V] |
| Steinfeldt, R., Rhein, M., Kieke, D. (2024). Anthropogenic carbon storage and its decadal changes in the Atlantic between 1990–2020. *Biogeosciences* 21, 3839–3867. doi:10.5194/bg-21-3839-2024 | modified TTD with dilution; CFCs, SF6, tritium | Atlantic inventory 43.0 → 68.2 Pg C (1990→2020); deviations from steady state strongest in subpolar N. Atlantic and Southern Ocean; deep S. Atlantic changes small | [V] |
| Carter, B. R., et al. (2025). Tracer-based Rapid Anthropogenic Carbon Estimation (TRACE). *ESSD* 17, 3073–3088. doi:10.5194/essd-17-3073-2025 | inverse-Gaussian TTD; α fitted to CFC-11/12, SF6 (GLODAPv2.2023) + OCIM ideal age; NN for preformed properties; inputs lon, lat, depth, S, T, year | Not trained on observed DIC; steady circulation; ±4.4 µmol kg⁻¹ or 15 % (1σ) reconstruction error; likely overestimates in upwelling surface waters; lower deep inventories than GLODAPv2 mapped below ~500 m; no direct comparison with observed DIC trends reported | [V] |
| Sandborn, D. E., Carter, B. R., Warner, M. J., Dias, L. M. (2026). TRACE-Python. *GMD* 19, 5961–5978. doi:10.5194/gmd-19-5961-2026; software doi:10.5281/zenodo.17822675 | Python port (MIT) | — | [V] |
| Keppler, L., et al. (2023). Recent trends and variability in the oceanic storage of dissolved inorganic carbon. *GBC* 37, e2022GB007677. doi:10.1029/2022GB007677 | MOBO-DIC ML mapping, 2004–2019, upper 1500 m | DIC pool +3.2 ± 0.7 Pg C yr⁻¹ (extrapolated to 4000 m), statistically indistinguishable from Cant uptake; interannual variability largest in tropical Pacific (ENSO) | [A] (publisher page blocked; abstract/secondary sources) |
| Boteler, C., et al. (2023). Trends of anthropogenic DIC in the Northwest Atlantic … state-space model. *JGR Oceans*. | statistical decomposition, regional | regime change ~2000; strong summer sampling bias | [A] (title/author from preprint PDF; DOI not verified) |
| DeVries, T., et al. (2023). Magnitude, trends, and variability of the global ocean carbon sink from 1985 to 2018. *GBC* 37, e2023GB007780. doi:10.1029/2023GB007780 | RECCAP2 synthesis | context | [A] |
| Gruber, N., Bakker, D. C. E., DeVries, T., et al. (2023). Trends and variability in the ocean carbon sink. *Nat. Rev. Earth Environ.* 4, 119–134. doi:10.1038/s43017-022-00381-x | review | context | [M] |

## 3. Empirical algorithms for interior DIC

| Ref | Predictors / time handling | Tag |
|---|---|---|
| Carter, B. R., et al. (2018). Updated methods for global locally interpolated estimation of alkalinity, pH, and nitrate. *L&O Methods* 16, 119–131. doi:10.1002/lom3.10232 | LIR; location-varying linear regressions | [M] |
| Carter, B. R., et al. (2021). New and updated global empirical seawater property estimation routines. *L&O Methods* 19, 785–809. doi:10.1002/lom3.10461 | ESPER (LIR, NN, mixed); includes a first-principles adjustment of DIC/pH for anthropogenic change to a reference year | [V] |
| Bittig, H. C., et al. (2018). An alternative to static climatologies: Robust estimation of open ocean CO2 variables and nutrient concentrations from T, S, and O2 data using Bayesian neural networks. *Front. Mar. Sci.* 5, 328. doi:10.3389/fmars.2018.00328 | CANYON-B; inputs include time, position, pressure, T, S, O2 | [V] |
| Broullón, D., et al. (2020). A global monthly climatology of oceanic total dissolved inorganic carbon: a neural network approach. *ESSD* 12, 1725–1743. doi:10.5194/essd-12-1725-2020 | NN; includes year/CO2 term | [M] |
| Cranmer, M. (2023). Interpretable machine learning for science with PySR and SymbolicRegression.jl. arXiv:2305.01582 | symbolic regression tool used in the repository | [M] |

No published symbolic-regression study of ocean DIC was found (one targeted search) — the repository’s SR framing is novel in tooling, but its substantive finding (smooth low-order forms predict DIC within basins) is long established by the algorithms above.

## 4. AOU as a respiration tracer

| Ref | Content | Tag |
|---|---|---|
| Ito, T., Follows, M. J., Boyle, E. A. (2004). Is AOU a good measure of respiration in the oceans? *GRL* 31, L17305. doi:10.1029/2004GL020900 | Model-based: surface O2 disequilibrium in deep-water formation regions decouples AOU from true O2 utilisation; Southern Ocean largest | [V] (magnitude 50–100 µmol kg⁻¹ is from citing literature, [A]) |

## 5. Validation of models on clustered geoscience data

| Ref | Content | Tag |
|---|---|---|
| Roberts, D. R., et al. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography* 40, 913–929. doi:10.1111/ecog.02881 | blocked CV rationale | [M] |
| Ploton, P., et al. (2020). Spatial validation reveals poor predictive performance of large-scale ecological mapping models. *Nat. Commun.* 11, 4540. doi:10.1038/s41467-020-18321-y | random CV over-optimistic under spatial autocorrelation | [V] |
| Kapoor, S., & Narayanan, A. (2023). Leakage and the reproducibility crisis in machine-learning-based science. *Patterns* 4, 100804. doi:10.1016/j.patter.2023.100804 | leakage taxonomy | [M] |

---

## 6. Gap analysis

| Direction | What exists | What remains | Verdict |
|---|---|---|---|
| Static S–T–AOU DIC compression (existing repo) | LIR/ESPER/CANYON-B achieve ≲ 10 µmol kg⁻¹; Cant time handling built in | Differences among smooth forms are below measurement consistency | **Low novelty**; also confounded by Cant trend |
| Cant decadal change | eMLR(C*) (Gruber 2019; Müller 2023), TTD (Steinfeldt 2024; TRACE 2025), ML DIC (Keppler 2023) | eMLR differences two occupations/eras; TTD is tracer-based and assumes steady circulation; *no study found that applies one identical, within-location estimator to observed DIC and to TTD-predicted Cant at the same samples over the full 1972–2021 record* | **Credible gap** — a direct, like-for-like observational consistency test of TTD-based accumulation rates, by basin × density layer |
| Non-Cant DIC change at fixed AOU | Müller 2023 infer natural carbon loss from storage vs flux mismatch; Keppler 2023 attribute trend mainly to Cant | Residual (observed − TTD) DIC trend at fixed θ, S, AOU is a direct, if model-dependent, diagnostic of non-Cant DIC change or of TTD/TSS bias | **Gap**, but attribution is ambiguous (TRACE error vs real natural change) |
| Validation pitfalls for empirical DIC models | Community normalises to a reference year (GLODAPv2 mapped, ESPER) — so experts know the trend matters | Quantifying how an unmodelled trend changes temporal-holdout skill and model *rankings* is not documented for DIC regressions | **Modest methodological gap** (secondary contribution) |
| Southern Ocean AOU disequilibrium effects on DIC regressions | Ito 2004 (model); later work | Abyssal failure of the repo’s regressions is plausibly AOU disequilibrium + water-mass mixing | Known mechanism; incremental |

**Novelty caveat.** “Not found” is based on a limited search (≈15 queries; some publishers blocked). Before submission, a systematic search (Web of Science / Scopus: “TRACE” OR “transit time distribution” AND “DIC trend”; “eMLR” AND “TTD” comparison; Tanhua 2007 follow-ups; Carter et al. 2017 Pacific DIC change; Woosley et al. 2016 Indian Ocean; Pérez et al. Atlantic) is required. Known closest neighbours: the eMLR(C*) vs TTD comparisons inside Steinfeldt et al. (2024) and the TRACE paper (inventory-level, not trend-level, not sample-matched).

---
## Round-2 additions (2026-11-14)

**Verification upgrade.** All nine Round-1 [M] references were checked against Crossref metadata (title, journal, volume, pages, year): Olsen 2016, Sabine 2004, Friis 2005, Gruber 2023, Carter 2018, Broullón 2020, Roberts 2017, Kapoor 2023, Ito 2004. All match; they are retagged **[V-meta]**. Bibliographic details are verified; the cited claims rest on prior knowledge or abstracts.

| Ref | Relevance | Tag |
|---|---|---|
| Sandborn, D. E., et al. (2026). TRACE-Python … *GMD* 19, 5961–5978. doi:10.5194/gmd-19-5961-2026 | Full text read (§4.1). Comparison with DIC-based estimates (Sabine 2004; Gruber 2019; Müller 2023) is **inventory-level** and qualitative ("good"). The accumulation-rate comparison is with Davila et al. (2022) and DeVries (2014), not with decadal DIC-based rates. **No sample-level or rate-level comparison with observed DIC.** | [V] |
| Müller, J. D., et al. (2023) data, NCEI 0279447 | Gridded ΔC_ant 1994–2004 and 2004–2014; 11 cases; standard = mask "3", "Standard case". The README was read in full. | [V] |
| Carter, B. R., et al. (2017). Two decades of Pacific anthropogenic carbon storage … P16 and P02. *GBC* 31, 306–327. doi:10.1002/2016GB005485 | eMLR along repeat sections; C_ant increase reaches ~1000 m | [V-meta] + [A] |
| Sabine, C. L., et al. (2008). Decadal changes in Pacific carbon. *JGR Oceans* 113. doi:10.1029/2007JC004577 | AOU-related circulation term can dominate DIC change in the North Pacific thermocline. This is relevant to our tracer-dependent Pacific-thermocline result (E08). | [V-meta] + [A] |
| Carter, B. R., et al. (2024). Random and systematic uncertainty in ship-based seawater carbonate chemistry observations. *L&O* 69, 2473–2488. doi:10.1002/lno.12674 | Systematic (cruise-level) uncertainty in GLODAP carbonate data; GLODAP adjustments improve consistency for some subsets and worsen it for others. Motivates the drift-immune design. | [A] |
| Gruber, N., Sarmiento, J. L., Stocker, T. F. (1996). *GBC* 10, 809–837. doi:10.1029/96GB01608 | C* definition used in E08 | [V-meta] |
| Bates, N. R., et al. (2014). *Oceanography* 27(1), 126–141. doi:10.5670/oceanog.2014.16 | Seven surface time series (Table 2). Surface pCO₂ trends: BATS 1.69 ± 0.11, HOT 1.72 ± 0.09, ESTOC 1.92 ± 0.24 µatm yr⁻¹; Iceland Sea 1.29 ± 0.36; Munida 1.28 ± 0.33; Irminger 2.37 ± 0.49; CARIACO 2.95 ± 0.43. Text (p. 135): at ESTOC, HOT and BATS "the increase in surface pCO2 was similar to the increase observed in atmospheric pCO2", and ΔpCO₂ "has remained fairly constant across the subtropical gyres"; Iceland Sea and Munida lag the atmosphere. **Full text read 2026-11-14 (Round 3).** | [V-full] |
| Hong, F. & Zemskova, V. E. (2026). Changes in oceanic carbon storage due to anthropogenic carbon input over the past three decades. *Atmosphere-Ocean* 64, 296–319. doi:10.1080/07055900.2026.2660361 | Annual C* (Gruber 1996 tracer) from the **ECCO-Darwin** ocean state estimate on a regular grid, 1992–2022: +60 Pg C (≈28 % of emissions); 36 Pg C over 2004–2019 vs 44–49 Pg C from Gruber 2019. Accelerating accumulation; EOF modes linked to Pacific climate modes. **No** TTD/TRACE comparison, no sample-matched observational trend estimation, and no treatment of cruise offsets, so no overlap with our contribution. **Read in full from EarthArXiv preprint v2 (doi:10.31223/X5QX88, 16 Dec 2025); the published version is paywalled (403), and differences from the preprint are not checked.** | [V-full preprint] |
| Boteler, C., et al. (2023). *JGR Oceans* 128. doi:10.1029/2022JC019483 | Isopycnal-layer C_ant trends in the NW Atlantic (state-space); ~0.57 µmol kg⁻¹ yr⁻¹ across layers 1993–2015 | [V-meta] + [A] |

**Novelty re-check (Round 2).** Eight further targeted searches were run, covering TRACE validation against DIC trends, TTD vs repeat-hydrography rate comparisons, fixed-effects or isopycnal DIC trends, and GLODAP drift. None found:
- a sample-matched comparison of DIC trends with TTD C_ant rates; or
- an estimator that removes cruise-level offsets (cruise fixed effects) from interior DIC trend estimates.

The closest are:
- the inventory-level comparisons in the TRACE papers;
- Steinfeldt et al. (2024): TTD vs eMLR, Atlantic only, inventory-level;
- Hong & Zemskova (2026): model-based (ECCO-Darwin) C* storage; read in preprint form (Round 3); no overlap.

The novelty claim stands at "not found in a targeted search". A Web of Science/Scopus search and expert consultation are still required.
