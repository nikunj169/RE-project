"""
Generate all additional data needed for the revised manuscript.
"""
import numpy as np
import scipy.io as sio
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression
import os, warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Load data
for path in (os.path.join(SCRIPT_DIR, 'data/raw/GLODAPv2.2023_Merged_Master_File.mat'),
             os.path.join(SCRIPT_DIR, 'GLODAPv2.2023_Merged_Master_File.mat')):
    try:
        mat = sio.loadmat(path, squeeze_me=True); break
    except FileNotFoundError:
        continue

S_all = mat['G2salinity'].astype(float)
T_all = mat['G2temperature'].astype(float)
A_all = mat['G2aou'].astype(float)
Y_all = mat['G2tco2'].astype(float)
Reg   = mat['G2region'].astype(float)
Lat   = mat['G2latitude'].astype(float)
Lon   = mat['G2longitude'].astype(float)
Dep   = mat['G2depth'].astype(float)
Yr    = mat['G2year'].astype(float)

qc = (np.isfinite(S_all) & np.isfinite(T_all) & np.isfinite(A_all)
      & np.isfinite(Y_all) & np.isfinite(Dep) & np.isfinite(Yr)
      & (S_all > 25) & (S_all < 42) & (T_all > -2.5) & (T_all < 35)
      & (Y_all > 1700) & (Y_all < 2600) & (A_all > -50))

SO_LAT = -35; TEMPORAL = 2015; EXTERNAL = 2018

BASINS = {
    'Atlantic':      {'codes': [1],    'lat_min': SO_LAT},
    'Indian':        {'codes': [16],   'lat_min': SO_LAT},
    'Pacific':       {'codes': [8],    'lat_min': SO_LAT},
    'Southern Ocean':{'codes': None,   'lat_max': SO_LAT},
}

def get_mask(name, cfg):
    if name == 'Southern Ocean':
        return qc & (Lat < cfg['lat_max'])
    return qc & np.isin(Reg, cfg['codes']) & (Lat > cfg['lat_min'])

def quad_model(X, a, b, g, d, e):
    S, T, A = X
    return (a*S + b*T + g*(A/100.) + d)**2 + e

def mlr_feats(S, T, A):
    return np.column_stack([S, T, A, S**2, T**2, A**2, S*T, S*A, T*A])

# ══════════════════════════════════════════════════════════════════════════════
# 1. Basin predictor summary statistics
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("1. Basin predictor summary statistics")
print("=" * 70)

pred_stats = []
for bname, cfg in BASINS.items():
    m = get_mask(bname, cfg)
    S, T, A, Y, dp = S_all[m], T_all[m], A_all[m], Y_all[m], Dep[m]
    for var_name, var in [('S', S), ('T', T), ('AOU', A), ('TCO2', Y), ('Depth', dp)]:
        pred_stats.append({
            'Basin': bname, 'Variable': var_name, 'N': len(var),
            'Mean': round(np.mean(var), 2), 'SD': round(np.std(var), 2),
            'Min': round(np.min(var), 2), 'Q25': round(np.percentile(var, 25), 2),
            'Median': round(np.median(var), 2), 'Q75': round(np.percentile(var, 75), 2),
            'Max': round(np.max(var), 2),
        })

df_pred = pd.DataFrame(pred_stats)
df_pred.to_csv(os.path.join(SCRIPT_DIR, 'revision_predictor_stats.csv'), index=False)
print(df_pred.to_string(index=False))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Predictor correlations per basin
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("2. Predictor correlations per basin")
print("=" * 70)

corr_stats = []
for bname, cfg in BASINS.items():
    m = get_mask(bname, cfg)
    S, T, A = S_all[m], T_all[m], A_all[m]
    r_ST, _ = stats.pearsonr(S, T)
    r_SA, _ = stats.pearsonr(S, A)
    r_TA, _ = stats.pearsonr(T, A)
    corr_stats.append({
        'Basin': bname,
        'r(S,T)': round(r_ST, 4),
        'r(S,AOU)': round(r_SA, 4),
        'r(T,AOU)': round(r_TA, 4),
    })
    print(f"  {bname}: r(S,T)={r_ST:.4f}, r(S,AOU)={r_SA:.4f}, r(T,AOU)={r_TA:.4f}")

df_corr = pd.DataFrame(corr_stats)
df_corr.to_csv(os.path.join(SCRIPT_DIR, 'revision_predictor_correlations.csv'), index=False)

# ══════════════════════════════════════════════════════════════════════════════
# 3. Derivative analysis across predictor domain
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("3. Derivative analysis across predictor domain")
print("=" * 70)

deriv_results = []
p0 = [1.0, -0.3, 2.5, 20.0, 2000.0]

for bname, cfg in BASINS.items():
    m = get_mask(bname, cfg)
    S, T, A, Y, yr = S_all[m], T_all[m], A_all[m], Y_all[m], Yr[m]
    tr = yr < TEMPORAL

    Xtr = (S[tr], T[tr], A[tr]); Ytr = Y[tr]
    popt, _ = curve_fit(quad_model, Xtr, Ytr, p0=p0, maxfev=30000)

    alpha, beta, gamma, delta, eps = popt
    core = alpha*S + beta*T + gamma*(A/100) + delta

    dS = 2 * alpha * core
    dT = 2 * beta * core
    dA = 2 * gamma * (1/100) * core

    # Fraction of domain where derivative has expected sign
    frac_dS_pos = np.mean(dS > 0) * 100
    frac_dT_neg = np.mean(dT < 0) * 100
    frac_dA_pos = np.mean(dA > 0) * 100

    deriv_results.append({
        'Basin': bname,
        'frac_dS_positive': round(frac_dS_pos, 2),
        'frac_dT_negative': round(frac_dT_neg, 2),
        'frac_dA_positive': round(frac_dA_pos, 2),
        'dS_median': round(np.median(dS), 2),
        'dT_median': round(np.median(dT), 2),
        'dA_median': round(np.median(dA), 4),
        'core_median': round(np.median(core), 2),
        'core_min': round(np.min(core), 2),
        'core_max': round(np.max(core), 2),
    })
    print(f"  {bname}: core=[{np.min(core):.1f}, {np.median(core):.1f}, {np.max(core):.1f}]")
    print(f"    dS>0: {frac_dS_pos:.1f}%, dT<0: {frac_dT_neg:.1f}%, dA>0: {frac_dA_pos:.1f}%")

df_deriv = pd.DataFrame(deriv_results)
df_deriv.to_csv(os.path.join(SCRIPT_DIR, 'revision_derivative_analysis.csv'), index=False)

# ══════════════════════════════════════════════════════════════════════════════
# 4. TOST on per-sample loss differentials (correct approach)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("4. Corrected TOST on per-sample absolute errors")
print("=" * 70)

tost_corrected = []
for bname, cfg in BASINS.items():
    m = get_mask(bname, cfg)
    S, T, A, Y, yr = S_all[m], T_all[m], A_all[m], Y_all[m], Yr[m]
    tr = yr < TEMPORAL; te = (yr >= TEMPORAL) & (yr < EXTERNAL)
    if te.sum() < 100: continue

    Xtr = (S[tr], T[tr], A[tr]); Ytr = Y[tr]
    Xte = (S[te], T[te], A[te]); Yte = Y[te]

    popt, _ = curve_fit(quad_model, Xtr, Ytr, p0=p0, maxfev=30000)
    Yp_q = quad_model(Xte, *popt)
    Yp_m = LinearRegression().fit(mlr_feats(*Xtr), Ytr).predict(mlr_feats(*Xte))

    # Per-sample absolute errors
    ae_quad = np.abs(Yte - Yp_q)
    ae_mlr = np.abs(Yte - Yp_m)

    # Paired difference: AE_MLR - AE_Quad (positive = quad better)
    d = ae_mlr - ae_quad
    n = len(d)
    d_bar = np.mean(d)
    d_sd = np.std(d, ddof=1)
    se = d_sd / np.sqrt(n)

    for delta_name, delta in [('2umol', 2.0), ('5umol', 5.0)]:
        # TOST on paired differences
        # H0: |mean(d)| >= delta  vs  H1: |mean(d)| < delta
        # Test 1: d_bar > -delta
        t1 = (d_bar - (-delta)) / se
        p1 = stats.t.cdf(t1, df=n-1)
        # Test 2: d_bar < +delta
        t2 = (d_bar - delta) / se
        p2 = 1 - stats.t.cdf(t2, df=n-1)
        # TOST p-value
        tost_p = max(p1, p2)
        equiv = tost_p < 0.05

        # CI
        tcrit = stats.t.ppf(0.975, df=n-1)
        ci_lo = d_bar - tcrit * se
        ci_hi = d_bar + tcrit * se

        ci_within = (ci_lo > -delta) and (ci_hi < delta)

        tost_corrected.append({
            'Basin': bname, 'Margin': delta_name, 'Delta': delta,
            'Mean_AE_diff': round(d_bar, 4), 'SE': round(se, 4),
            'TOST_p': round(tost_p, 6), 'Equiv_by_test': 'YES' if equiv else 'NO',
            'CI_lo': round(ci_lo, 4), 'CI_hi': round(ci_hi, 4),
            'CI_within_bounds': 'YES' if ci_within else 'NO',
            'N': n,
        })
        print(f"  {bname} (δ={delta}): d_bar={d_bar:+.4f}, CI=[{ci_lo:+.4f}, {ci_hi:+.4f}], "
              f"TOST p={tost_p:.6f}, test={'EQUIV' if equiv else 'NOT'}, "
              f"CI within={'YES' if ci_within else 'NO'}")

df_tost = pd.DataFrame(tost_corrected)
df_tost.to_csv(os.path.join(SCRIPT_DIR, 'revision_tost_corrected.csv'), index=False)

# ══════════════════════════════════════════════════════════════════════════════
# 5. Southern Ocean abyssal analysis
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("5. Southern Ocean abyssal analysis (3000-7000m)")
print("=" * 70)

m_so = get_mask('Southern Ocean', BASINS['Southern Ocean'])
S_so = S_all[m_so]; T_so = T_all[m_so]; A_so = A_all[m_so]
Y_so = Y_all[m_so]; yr_so = Yr[m_so]; dp_so = Dep[m_so]

abyss = dp_so >= 3000
tr_abyss = abyss & (yr_so < TEMPORAL)
te_abyss = abyss & (yr_so >= TEMPORAL) & (yr_so < EXTERNAL)

print(f"  Abyssal SO: train={tr_abyss.sum()}, test={te_abyss.sum()}")
print(f"  Test predictor ranges:")
print(f"    S: [{S_so[te_abyss].min():.2f}, {S_so[te_abyss].max():.2f}], mean={S_so[te_abyss].mean():.2f}")
print(f"    T: [{T_so[te_abyss].min():.2f}, {T_so[te_abyss].max():.2f}], mean={T_so[te_abyss].mean():.2f}")
print(f"    AOU: [{A_so[te_abyss].min():.2f}, {A_so[te_abyss].max():.2f}], mean={A_so[te_abyss].mean():.2f}")
print(f"    TCO2: [{Y_so[te_abyss].min():.2f}, {Y_so[te_abyss].max():.2f}], mean={Y_so[te_abyss].mean():.2f}")

# Fit model
Xtr_a = (S_so[tr_abyss], T_so[tr_abyss], A_so[tr_abyss])
Ytr_a = Y_so[tr_abyss]
Xte_a = (S_so[te_abyss], T_so[te_abyss], A_so[te_abyss])
Yte_a = Y_so[te_abyss]

popt_a, _ = curve_fit(quad_model, Xtr_a, Ytr_a, p0=p0, maxfev=30000)
Yp_a = quad_model(Xte_a, *popt_a)
res_a = Yte_a - Yp_a

# Also fit MLR for comparison
mlr_a = LinearRegression().fit(mlr_feats(*Xtr_a), Ytr_a)
Yp_mlr_a = mlr_a.predict(mlr_feats(*Xte_a))
res_mlr_a = Yte_a - Yp_mlr_a

print(f"  Quad: R²={r2_score(Yte_a, Yp_a):.4f}, RMSE={np.sqrt(mean_squared_error(Yte_a, Yp_a)):.2f}")
print(f"  MLR:  R²={r2_score(Yte_a, Yp_mlr_a):.4f}, RMSE={np.sqrt(mean_squared_error(Yte_a, Yp_mlr_a)):.2f}")

# Mean reference (predict mean)
mean_pred = np.full_like(Yte_a, np.mean(Ytr_a))
rmse_mean = np.sqrt(mean_squared_error(Yte_a, mean_pred))
print(f"  Mean reference RMSE: {rmse_mean:.2f}")
print(f"  R² negative means model is worse than predicting the mean")

# Residual patterns
print(f"\n  Residual statistics (Quad):")
print(f"    Mean bias: {np.mean(res_a):.2f}")
print(f"    SD: {np.std(res_a):.2f}")
print(f"    Median: {np.median(res_a):.2f}")

# Correlation of residuals with predictors
r_res_T, _ = stats.pearsonr(T_so[te_abyss], res_a)
r_res_S, _ = stats.pearsonr(S_so[te_abyss], res_a)
r_res_A, _ = stats.pearsonr(A_so[te_abyss], res_a)
print(f"    Corr(residual, T): {r_res_T:.4f}")
print(f"    Corr(residual, S): {r_res_S:.4f}")
print(f"    Corr(residual, AOU): {r_res_A:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# 6. MLR parameter count verification
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("6. MLR parameter count verification")
print("=" * 70)
print("  mlr_feats produces: S, T, AOU, S², T², AOU², ST, SAOU, TAOU")
print("  = 9 features + 1 intercept = 10 parameters total")

# Verify
test_S = np.array([35.0, 36.0])
test_T = np.array([10.0, 15.0])
test_A = np.array([200.0, 250.0])
feats = mlr_feats(test_S, test_T, test_A)
print(f"  mlr_feats shape: {feats.shape} = (2 samples, {feats.shape[1]} features)")
print(f"  With intercept: {feats.shape[1] + 1} parameters")

# ══════════════════════════════════════════════════════════════════════════════
# 7. Rank-1 geometry data (core variable)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("7. Rank-1 geometry: TCO2 vs core variable")
print("=" * 70)

for bname, cfg in BASINS.items():
    m = get_mask(bname, cfg)
    S, T, A, Y, yr = S_all[m], T_all[m], A_all[m], Y_all[m], Yr[m]
    tr = yr < TEMPORAL; te = (yr >= TEMPORAL) & (yr < EXTERNAL)
    if te.sum() < 100: continue

    Xtr = (S[tr], T[tr], A[tr]); Ytr = Y[tr]
    popt, _ = curve_fit(quad_model, Xtr, Ytr, p0=p0, maxfev=30000)

    core = popt[0]*S[te] + popt[1]*T[te] + popt[2]*(A[te]/100) + popt[3]
    r_core, _ = stats.pearsonr(core, Y[te])
    print(f"  {bname}: r(core, TCO2)={r_core:.4f}")

print("\nAll revision data generated.")
