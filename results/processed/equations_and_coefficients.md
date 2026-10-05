## Equation forms

Let $Y=\mathrm{TCO}_2$, $S$ be salinity, $T$ be in-situ temperature, and $A$ be AOU. For Models 0–3, $a=A/100$; Models 4–5 use raw AOU ($A$), matching the implementation.

- **Mean baseline:** $Y=m$
- **Linear:** $Y=b_0+b_S S+b_T T+b_A a$
- **Rank-1 quadratic:** $Y=(\alpha S+\beta T+\gamma a+\delta)^2+\epsilon$
- **Rank-2 quadratic:** $Y=(a_1S+b_1T+c_1a+d_1)^2+(a_2S+b_2T+c_2a+d_2)^2+\epsilon$
- **Full second-order polynomial:** $Y=b_0+b_S S+b_T T+b_A A+b_{S2}S^2+b_{T2}T^2+b_{A2}A^2+b_{ST}ST+b_{SA}SA+b_{TA}TA$
- **Cubic polynomial:** $Y=b_0+b_0' S+b_1' T+b_2' A+b_3'S^2+b_4'T^2+b_5'A^2+b_6'ST+b_7' SA+b_8'TA+b_9'S^3+b_{10}'T^3+b_{11}'A^3+b_{12}'S^2T+b_{13}'S^2A+b_{14}'T^2S+b_{15}'T^2A+b_{16}'A^2S+b_{17}'A^2T+b_{18}'STA$

For the cubic model, the coefficient labels `coef_0`–`coef_18` follow this term order: $S,T,A,S^2,T^2,A^2,ST,SA,TA,S^3,T^3,A^3,S^2T,S^2A,T^2S,T^2A,A^2S,A^2T,STA$.

## Fitting and dataset summary

- **Dataset:** GLODAP v2.2023 observations.
- **Predictors:** salinity, in-situ temperature, and apparent oxygen utilisation (AOU). **Target:** TCO₂.
- **Quality control:** finite values with 25 < salinity < 42, −2.5 < temperature < 35 °C, 1700 < TCO₂ < 2600 μmol kg⁻¹, and AOU > −50 μmol kg⁻¹.
- **Basin fits:** Atlantic, Indian, Pacific, and Southern Ocean; models were fitted separately within each basin.
- **Temporal split:** training years < 2015; validation years 2015–2017; external holdout years ≥ 2018.
- **Fitting:** linear and polynomial models use least-squares linear regression; rank-1 uses nonlinear least squares; rank-2 uses nonlinear optimization with multiple starts.
- **Interpretation:** coefficients are empirical, basin-specific fitted parameters, not universal physical constants.

## Coefficient values

### Atlantic

#### Mean baseline

| Parameter | Coefficient |
|---|---:|
| `mean` | `2152.61537722` |

#### Linear

| Parameter | Coefficient |
|---|---:|
| `intercept` | `901.83859811` |
| `coef_S` | `35.8575138152` |
| `coef_T` | `-6.10881812257` |
| `coef_A` | `65.9414573712` |

#### Rank-1 quadratic

| Parameter | Coefficient |
|---|---:|
| `alpha` | `1.48552371271` |
| `beta` | `-0.237665343548` |
| `gamma` | `2.21964093907` |
| `delta` | `-36.610677187` |
| `epsilon` | `1923.45373189` |

#### Rank-2 quadratic

| Parameter | Coefficient |
|---|---:|
| `a1` | `1.76340316668` |
| `b1` | `-0.0939230769796` |
| `c1` | `-1.00050960825` |
| `d1` | `-82.8531109439` |
| `a2` | `1.30292543394` |
| `b2` | `-0.119524859662` |
| `c2` | `0.243061486251` |
| `d2` | `-0.130045531898` |
| `epsilon` | `-355.022439492` |

#### Full second-order polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `6510.48538479` |
| `coef_S` | `-289.458932222` |
| `coef_T` | `10.7280861014` |
| `coef_A` | `0.970342825742` |
| `coef_S2` | `4.71066260265` |
| `coef_T2` | `-0.0506885667229` |
| `coef_A2` | `4.72039243766e-05` |
| `coef_ST` | `-0.446054705175` |
| `coef_SA` | `-0.0073119372594` |
| `coef_TA` | `-0.00694172440054` |

#### Cubic polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1937.29844319` |
| `coef_0` | `-0.00643289262686` |
| `coef_1` | `-0.0336834502529` |
| `coef_2` | `0.016944028983` |
| `coef_3` | `-0.226012160969` |
| `coef_4` | `-0.384449576482` |
| `coef_5` | `-0.0157831897178` |
| `coef_6` | `-0.625986233753` |
| `coef_7` | `0.160922795735` |
| `coef_8` | `-0.789845950847` |
| `coef_9` | `0.0108197549799` |
| `coef_10` | `0.0137540590558` |
| `coef_11` | `-1.99269958875e-07` |
| `coef_12` | `0.0210619115343` |
| `coef_13` | `-0.00386195285313` |
| `coef_14` | `-0.00943319569205` |
| `coef_15` | `0.00167218447252` |
| `coef_16` | `0.000450502903651` |
| `coef_17` | `7.98244087982e-05` |
| `coef_18` | `0.020699158751` |

### Indian

#### Mean baseline

| Parameter | Coefficient |
|---|---:|
| `mean` | `2194.23658071` |

#### Linear

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1305.42684569` |
| `coef_S` | `25.3481755827` |
| `coef_T` | `-8.8108561039` |
| `coef_A` | `69.6481121197` |

#### Rank-1 quadratic

| Parameter | Coefficient |
|---|---:|
| `alpha` | `0.501000054944` |
| `beta` | `-0.169608163253` |
| `gamma` | `1.28779343415` |
| `delta` | `9.95496355353` |
| `epsilon` | `1435.88910714` |

#### Rank-2 quadratic

| Parameter | Coefficient |
|---|---:|
| `a1` | `0.0696608295062` |
| `b1` | `-0.0596366723258` |
| `c1` | `0.669162508271` |
| `d1` | `105.944237171` |
| `a2` | `-0.34186692224` |
| `b2` | `-0.121373496848` |
| `c2` | `2.62386946091` |
| `d2` | `-4.25493904205` |
| `epsilon` | `-9817.70608309` |

#### Full second-order polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `-5787.20641324` |
| `coef_S` | `383.117939265` |
| `coef_T` | `28.5586939788` |
| `coef_A` | `3.39757783859` |
| `coef_S2` | `-4.41187045147` |
| `coef_T2` | `0.149358593742` |
| `coef_A2` | `0.000353644914068` |
| `coef_ST` | `-1.19712907912` |
| `coef_SA` | `-0.0793531137247` |
| `coef_TA` | `0.00126938955762` |

#### Cubic polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1694.29020115` |
| `coef_0` | `0.000298277330054` |
| `coef_1` | `-0.000743582055941` |
| `coef_2` | `0.0116477952901` |
| `coef_3` | `0.011015135494` |
| `coef_4` | `-0.0197757332807` |
| `coef_5` | `-0.0381225675205` |
| `coef_6` | `-0.010211048133` |
| `coef_7` | `0.22793308103` |
| `coef_8` | `0.158701815312` |
| `coef_9` | `0.010774024426` |
| `coef_10` | `0.00333142575763` |
| `coef_11` | `-5.57728133792e-06` |
| `coef_12` | `-0.0030260692558` |
| `coef_13` | `-0.00609874335649` |
| `coef_14` | `-0.00605447961991` |
| `coef_15` | `0.00272439549644` |
| `coef_16` | `0.00119094591325` |
| `coef_17` | `-3.55308955942e-05` |
| `coef_18` | `-0.00625535711323` |

### Pacific

#### Mean baseline

| Parameter | Coefficient |
|---|---:|
| `mean` | `2208.37655223` |

#### Linear

| Parameter | Coefficient |
|---|---:|
| `intercept` | `787.754412883` |
| `coef_S` | `40.2917212793` |
| `coef_T` | `-9.09816412525` |
| `coef_A` | `74.8288037433` |

#### Rank-1 quadratic

| Parameter | Coefficient |
|---|---:|
| `alpha` | `1.1593663523` |
| `beta` | `-0.269559477126` |
| `gamma` | `1.77552781742` |
| `delta` | `-20.152634285` |
| `epsilon` | `1790.28281032` |

#### Rank-2 quadratic

| Parameter | Coefficient |
|---|---:|
| `a1` | `-0.468819135966` |
| `b1` | `-0.315860989207` |
| `c1` | `1.16526853579` |
| `d1` | `29.2547755481` |
| `a2` | `0.889263424329` |
| `b2` | `-0.0593906037728` |
| `c2` | `0.759025380668` |
| `d2` | `-2.27691127406` |
| `epsilon` | `1207.32900882` |

#### Full second-order polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1776.32011457` |
| `coef_S` | `-3.1866377046` |
| `coef_T` | `-31.4396050975` |
| `coef_A` | `-1.66722765119` |
| `coef_S2` | `0.407139713008` |
| `coef_T2` | `0.0318484811347` |
| `coef_A2` | `-0.000353032591619` |
| `coef_ST` | `0.667595518211` |
| `coef_SA` | `0.0778619001323` |
| `coef_TA` | `-0.0210817627636` |

#### Cubic polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1588.58902345` |
| `coef_0` | `-5.77559973679e-05` |
| `coef_1` | `-0.00126306009848` |
| `coef_2` | `0.0162244107781` |
| `coef_3` | `-0.00176443087186` |
| `coef_4` | `-0.0188861148472` |
| `coef_5` | `-0.0267629645037` |
| `coef_6` | `-0.0221546663735` |
| `coef_7` | `0.266293445028` |
| `coef_8` | `0.088348621114` |
| `coef_9` | `0.0148743768008` |
| `coef_10` | `0.00706972233102` |
| `coef_11` | `-1.01114895316e-05` |
| `coef_12` | `-0.00650468675659` |
| `coef_13` | `-0.00767452238836` |
| `coef_14` | `-0.00587114747578` |
| `coef_15` | `0.00328859619244` |
| `coef_16` | `0.000937859177583` |
| `coef_17` | `-6.39632749395e-05` |
| `coef_18` | `-0.00367900964791` |

### Southern Ocean

#### Mean baseline

| Parameter | Coefficient |
|---|---:|
| `mean` | `2211.66109208` |

#### Linear

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1112.65262289` |
| `coef_S` | `30.7144327097` |
| `coef_T` | `-8.51800923798` |
| `coef_A` | `62.8712207251` |

#### Rank-1 quadratic

| Parameter | Coefficient |
|---|---:|
| `alpha` | `0.728905150259` |
| `beta` | `-0.230706263923` |
| `gamma` | `1.63784629084` |
| `delta` | `-6.20133788418` |
| `epsilon` | `1811.14938636` |

#### Rank-2 quadratic

| Parameter | Coefficient |
|---|---:|
| `a1` | `0.0685899882312` |
| `b1` | `-0.0951308378944` |
| `c1` | `2.62410765975` |
| `d1` | `2.4715254753` |
| `a2` | `-0.627555257278` |
| `b2` | `0.17033579823` |
| `c2` | `-0.611571640387` |
| `d2` | `-1.47915150915` |
| `epsilon` | `1614.62848998` |

#### Full second-order polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1659.46690988` |
| `coef_S` | `0.018205563555` |
| `coef_T` | `-8.03225759653` |
| `coef_A` | `0.268896963701` |
| `coef_S2` | `0.430344843807` |
| `coef_T2` | `-0.0321040590401` |
| `coef_A2` | `0.000742502142152` |
| `coef_ST` | `0.00572943892866` |
| `coef_SA` | `0.00794907267637` |
| `coef_TA` | `-0.0148189624888` |

#### Cubic polynomial

| Parameter | Coefficient |
|---|---:|
| `intercept` | `1754.12885308` |
| `coef_0` | `0.000123287926758` |
| `coef_1` | `-0.00293753145708` |
| `coef_2` | `0.000773656963118` |
| `coef_3` | `0.00386694578322` |
| `coef_4` | `-0.0374000328238` |
| `coef_5` | `0.0117464023005` |
| `coef_6` | `-0.0489852041431` |
| `coef_7` | `0.0378576409006` |
| `coef_8` | `-0.222530974493` |
| `coef_9` | `0.0101464116067` |
| `coef_10` | `0.0135893393152` |
| `coef_11` | `1.7751876149e-05` |
| `coef_12` | `-0.00393202333989` |
| `coef_13` | `-0.00047907377833` |
| `coef_14` | `-0.00943133615123` |
| `coef_15` | `-0.000500360549495` |
| `coef_16` | `-0.000435028122919` |
| `coef_17` | `-9.41455309576e-05` |
| `coef_18` | `0.00616977515977` |

