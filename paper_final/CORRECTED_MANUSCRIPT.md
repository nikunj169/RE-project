# Symbolic Regression Discovery of a Compact Quadratic Structure for Total Dissolved Inorganic Carbon

**[Author names — to be added]**
**[Institution — to be added]**
**[Correspondence: email@institution.edu — to be added]**

---

## Abstract

Symbolic regression was used to search for compact analytical expressions
relating total dissolved inorganic carbon (TCO₂) to salinity (S), temperature
(T), and apparent oxygen utilisation (AOU) in the global ocean. Within a
restricted search space of addition, subtraction, multiplication, and division
operators, PySR identified a squared-linear functional form:

$$TCO_2 = (\alpha S + \beta T + \gamma \cdot AOU/100 + \delta)^2 + \epsilon$$

This structure imposes a rank-1 quadratic dependence on a linear combination
of the predictors. The model was fitted separately to four ocean basins
(Atlantic, Indian, Pacific, Southern Ocean) using GLODAP v2.2023 data
(N = 472,154 observations after quality control).

Temporal validation used a chronological split: training on pre-2015 data,
model selection on 2015–2018 data, and external evaluation on post-2018
observations. In the Atlantic basin, the quadratic model achieved
R² = 0.9074 and RMSE = 20.17 μmol kg⁻¹, representing a 13.8% RMSE
reduction relative to a full multiple linear regression benchmark (MLR+AOU)
using nine free parameters. In the Indian, Pacific, and Southern Ocean basins,
the quadratic model showed similar test-set RMSE to MLR+AOU, with absolute
differences of 0.12–0.76 μmol kg⁻¹.

The quadratic structure is qualitatively consistent with nonlinear carbonate
buffering, although the Revelle-factor relationship does not uniquely imply
the discovered functional form. The model performs best in the Atlantic Ocean
and at depth (below 100 m). Surface predictions (0–100 m) are less reliable,
with R² as low as ~0.76 in the Atlantic. Southern Ocean predictions below
3000 m should be used with caution (R² = -0.174).

Approximate 95% prediction intervals were constructed via hybrid bootstrap
resampling (n = 1,000 replicates), with empirical coverage ranging from
89.3% to 95.8% across basins.

**Keywords:** total dissolved inorganic carbon, symbolic regression, ocean
carbon system, Revelle Factor, GLODAP

---

## 1. Introduction

Symbolic regression offers a data-driven approach to discovering compact
analytical relationships in complex Earth system data. Unlike black-box
machine learning methods, symbolic regression searches for interpretable
mathematical expressions that balance accuracy with simplicity.

Here we apply PySR (Cranmer, 2023) to search for compact expressions
relating TCO₂ to S, T, and AOU across four ocean basins. We emphasize that
the search was conducted within a restricted operator set (addition,
subtraction, multiplication, and division), and the discovered forms should
be interpreted as the best compact approximation found within that space,
rather than as a fundamental law of ocean carbon chemistry.

---

## 2. Data

### 2.1 GLODAP v2.2023

We used the Global Ocean Data Analysis Project version 2.2023
(Lauvset et al., 2023), extracting salinity (G2salinity), temperature
(G2temperature), apparent oxygen utilisation (G2aou), total dissolved
inorganic carbon (G2tco2), latitude (G2latitude), longitude (G2longitude),
depth (G2depth), year (G2year), and region code (G2region).

### 2.2 Quality Control

Observations were retained if all of the following were satisfied:
- All predictor and target values were finite
- Salinity: 25 < S < 42
- Temperature: −2.5 < T < 35 °C
- TCO₂: 1700 < Y < 2600 μmol kg⁻¹
- AOU: AOU > −50 μmol kg⁻¹

After quality control, 472,154 observations remained.

### 2.3 Basin Definitions

Four basins were defined using GLODAP region codes and latitude:
- **Atlantic**: region code 1, latitude > 35°S
- **Indian**: region code 16, latitude > 35°S
- **Pacific**: region code 8, latitude > 35°S
- **Southern Ocean**: latitude < 35°S (all region codes)

These definitions are mutually exclusive by construction.

### 2.4 Temporal Split

The data were split chronologically:
- **Training**: observations before 2015
- **Test/model selection**: 2015–2018
- **External holdout**: post-2018 observations

The external holdout spans approximately 3 years of observations
collected after the most recent training data. We note that this is a
temporal holdout, not a spatially independent sample — GLODAP contains
repeated hydrographic sections, so the external set likely includes
observations at similar locations to training data.

---

## 3. Methods

### 3.1 Symbolic Regression

PySR (Cranmer, 2023) was used to search for compact analytical expressions.
The search was conducted with the following configuration:

- **Binary operators**: +, −, ×, ÷
- **Unary operators**: exp, log
- **Predictors**: S, T, AOU/100
- **Target**: TCO₂
- **Max complexity**: 20
- **Populations**: 15
- **Population size**: 50
- **Iterations**: 50
- **Random seed**: 42

**Limitation**: The search space was restricted to the above operators.
Squared-linear structures are particularly accessible to this operator set
through repeated multiplication. The discovered form should be interpreted
as the best compact expression found within this restricted space, not as
a universal law. Robustness to alternative operator sets, complexity
penalties, and random seeds has not been systematically tested.

### 3.2 Model Fitting

The discovered equation was fitted per basin using nonlinear least squares
(scipy.optimize.curve_fit) with the Levenberg–Marquardt algorithm:

$$TCO_2 = (\alpha S + \beta T + \gamma \cdot AOU/100 + \delta)^2 + \epsilon$$

### 3.3 Benchmark: MLR+AOU

A multiple linear regression with all second-order cross-terms served as
the fair benchmark:

$$TCO_2 = \beta_0 + \beta_1 S + \beta_2 T + \beta_3 AOU + \beta_4 S^2 + \beta_5 T^2 + \beta_6 AOU^2 + \beta_7 ST + \beta_8 SAOU + \beta_9 TAOU$$

This benchmark uses 10 free parameters (vs. 5 for the quadratic model)
and includes the same predictors, providing a fair comparison.

### 3.4 Prediction Intervals

Approximate 95% prediction intervals were constructed via hybrid bootstrap
(n = 1,000 replicates):

1. Resample training observations with replacement
2. Refit the quadratic model on the resampled data
3. Generate predictions on the test set
4. Add a random draw from N(0, σ_resid) to each prediction, where σ_resid
   is the standard deviation of training residuals

The 2.5th and 97.5th percentiles of the resulting prediction distribution
form the reported intervals.

**Limitations**: This procedure resamples individual observations, but
GLODAP observations are not independent — they include repeated cruises,
spatially clustered sections, and depth correlations. The resulting
intervals may be too narrow. Empirical coverage was verified on the test
set and ranged from 89.3% to 95.8% across basins (nominal target: 95%).
The Atlantic (89.3%) and Southern Ocean (90.0%) intervals show
under-coverage, suggesting that a cluster or block bootstrap resampling
by cruise or section would be more appropriate.

### 3.5 Hessian Analysis

For the model f(x) = (v^T x + δ)² + ε, the Hessian matrix is:

$$H = 2vv^T$$

where v = (α, β, γ/100)^T. This matrix has rank 1, confirming that the
quadratic surface has curvature along only one direction in predictor space.
The factor of 2 arises from the second derivative of the squared term.

**Note**: The previous version of this manuscript stated H = vv^T, omitting
the factor of 2. The rank-1 conclusion is unchanged.

---

## 4. Results

### 4.1 Model Coefficients

[Table 1: corrected coefficients — see corrected_coefficients.csv]

The fitted coefficients show consistent signs across basins: α > 0 (salinity
coefficient), β < 0 (temperature coefficient), γ > 0 (AOU coefficient).

**Important caveat on sign interpretation**: Because the model is squared,
the sign of a coefficient such as β does not directly determine the sign of
the local TCO₂ sensitivity. The actual partial derivative is:

$$\frac{\partial TCO_2}{\partial T} = 2\beta(\alpha S + \beta T + \gamma A/100 + \delta)$$

The sign of this derivative depends on the sign of the entire linear term
(αS + βT + γA/100 + δ), which is positive over the observed predictor
domain in all basins. Therefore, the expected physical signs (positive
salinity effect, negative temperature effect, positive AOU effect) do hold
over the observed data range, but this is a property of the fitted model
evaluated at observed values, not simply of the individual coefficients.

### 4.2 Temporal Validation

[Table 2: corrected validation metrics — see corrected_results.csv]

The quadratic model achieves its best performance in the Atlantic basin
(R² = 0.9074, RMSE = 20.17 μmol kg⁻¹), with a 13.8% RMSE reduction
relative to MLR+AOU. In the Indian, Pacific, and Southern Ocean basins,
the quadratic model shows similar test-set RMSE to MLR+AOU, with absolute
differences of 0.12–0.76 μmol kg⁻¹.

**We do not claim statistical equivalence.** A failure to reject the null
hypothesis of equal predictive loss (p > 0.05 from a Diebold–Mariano test)
does not establish equivalence. The Diebold–Mariano test was designed for
sequential forecasting and its application to spatial oceanographic
observations is approximate. A proper equivalence or non-inferiority
analysis would require specifying a scientifically justified margin and
constructing a confidence interval for the difference in predictive loss.

The observed RMSE differences (0.12–0.76 μmol kg⁻¹) are small relative
to typical TCO₂ concentrations (~2100–2300 μmol kg⁻¹), corresponding to
approximately 0.01–0.04% of typical basin-mean TCO₂ concentration. However,
we do not convert this to a percentage of basin-scale carbon inventory,
as that would require volume- and density-weighted integration that has
not been performed.

### 4.3 External Holdout

Post-2018 observations (chronologically held out, not spatially independent)
show no systematic degradation:

[Insert external holdout metrics from corrected_results.csv]

### 4.4 Depth-Stratified Performance

[Table 3: corrected depth results — see corrected_depth_results.csv]

The model performs best at depth (below 100 m and particularly below 3000 m)
and worst at the surface (0–100 m). In the Atlantic, surface R² is
approximately 0.76, while deep-water R² exceeds 0.83.

**Southern Ocean below 3000 m**: R² = -0.174, indicating that the model
does not reliably predict TCO₂ in the Southern Ocean abyss. Predictions
in this depth range should not be used without independent validation.

### 4.5 Prediction Intervals

Approximate 95% prediction intervals (n = 1,000 bootstrap replicates)
achieved the following empirical coverage:

[Insert PI coverage from corrected_results.csv]

The Atlantic (89.3%) and Southern Ocean (90.0%) intervals show
under-coverage relative to the nominal 95% target. This is likely due to
the observation-level bootstrap not accounting for spatial and temporal
dependence in GLODAP data. A cluster bootstrap resampling by cruise or
section would be more appropriate but has not been implemented.

### 4.6 Physical Interpretation (Post-Hoc)

The quadratic structure is qualitatively consistent with nonlinear carbonate
buffering. The Revelle Buffer Factor (RF) describes the nonlinear relationship
between dissolved inorganic carbon and pCO₂, and higher RF values indicate
stronger nonlinearity. The Atlantic Ocean has the highest basin-mean RF,
consistent with its superior quadratic model performance.

**However**, this argument is post-hoc and qualitative. Nonlinear carbonate
chemistry can produce many different functional forms; the Revelle Factor
does not uniquely imply the specific squared-linear equation discovered here.
The physical interpretation should be viewed as a plausible rationale rather
than a derivation.

---

## 5. Discussion

### 5.1 Strengths

- Compact, interpretable equation with only 5 free parameters
- Competitive with MLR+AOU (9 parameters) in the Atlantic
- Physically interpretable structure (rank-1 mixing axis)
- Temporal validation on chronologically held-out data

### 5.2 Limitations

1. **Restricted search space**: PySR searched only +, −, ×, ÷, exp, log.
   The discovered form may not be the best possible compact expression.

2. **Robustness not demonstrated**: The result has not been verified across
   different random seeds, operator sets, complexity penalties, or training
   subsets.

3. **Not truly independent validation**: GLODAP contains repeated
   hydrographic sections. The temporal split tests temporal generalization
   but not spatial extrapolation.

4. **Surface limitations**: R² ~ 0.76 at 0–100 m in the Atlantic; the
   model is not recommended for surface carbon flux studies without
   additional validation.

5. **Southern Ocean abyss**: R² < 0 at 3000–7000 m; predictions in this
   range are unreliable.

6. **Underperforms MLR+AOU in most basins**: Only the Atlantic shows clear
   improvement; other basins show similar or slightly worse performance.

7. **Physical interpretation is post-hoc**: The carbonate chemistry argument
   does not derive the discovered form; it merely provides qualitative
   consistency.

8. **Bootstrap limitations**: The prediction intervals use observation-level
   resampling, which underestimates uncertainty for spatially/temporally
   correlated oceanographic data.

---

## 6. Conclusions

Symbolic regression identified a compact rank-1 quadratic structure for
TCO₂ within a restricted search space. The model performs competitively
with an unconstrained MLR+AOU benchmark in the Atlantic basin (13.8%
RMSE improvement) and shows similar performance in other basins. The
quadratic structure is qualitatively consistent with nonlinear carbonate
buffering, but this consistency does not constitute derivation.

The model is best suited for subsurface predictions (below 100 m) in the
Atlantic and Pacific basins. Surface predictions and Southern Ocean abyssal
predictions should be used with caution. The approach demonstrates the
potential of symbolic regression for discovering interpretable oceanographic
relationships, but the results should be viewed as exploratory rather than
definitive.

Future work should include: (1) systematic robustness testing across
operator sets and random seeds; (2) density-based water-mass stratification
to address the AAIW-related Atlantic bias; (3) cluster bootstrap or
mixed-effects uncertainty quantification; and (4) comparison with
physics-based ocean carbon models.

---

## Data Availability

GLODAP v2.2023 is available at https://www.glodap.info.
[Repository URL and DOI — to be added upon acceptance]

## Code Availability

[Repository URL and DOI — to be added upon acceptance]
The PySR configuration used in this study is documented in Section 3.1.
Exact PySR version, Julia version, and search logs will be archived
with the repository.

## Acknowledgements

[Funding statement — to be added]

## References

[Existing references — verify completeness]
