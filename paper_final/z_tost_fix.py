"""Corrected TOST calculation on per-sample absolute errors."""
import numpy as np
import scipy.io as sio
import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
import os, warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

for path in (os.path.join(SCRIPT_DIR, 'data/raw/GLODAPv2.2023_Merged_Master_File.mat'),):
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
    if name == 'Southern Ocean': return qc & (Lat < cfg['lat_max'])
    return qc & np.isin(Reg, cfg['codes']) & (Lat > cfg['lat_min'])

def quad_model(X, a, b, g, d, e):
    S, T, A = X; return (a*S + b*T + g*(A/100.) + d)**2 + e

def mlr_feats(S, T, A):
    return np.column_stack([S, T, A, S**2, T**2, A**2, S*T, S*A, T*A])

p0 = [1.0, -0.3, 2.5, 20.0, 2000.0]
results = []

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

    # Paired difference of absolute errors: positive = quad better
    d = np.abs(Yte - Yp_m) - np.abs(Yte - Yp_q)
    n = len(d)
    d_bar = np.mean(d)
    d_sd = np.std(d, ddof=1)
    se = d_sd / np.sqrt(n)

    print(f"\n{bname}: d_bar={d_bar:+.4f}, se={se:.4f}, n={n}")

    for delta_name, delta in [('2umol', 2.0), ('5umol', 5.0)]:
        # TOST: reject H0 of non-equivalence if BOTH one-sided tests reject
        # H01: d_bar <= -delta  vs  H11: d_bar > -delta
        t1 = (d_bar - (-delta)) / se
        p1 = 1 - stats.t.cdf(t1, df=n-1)  # P(T >= t1)

        # H02: d_bar >= +delta  vs  H12: d_bar < +delta
        t2 = (d_bar - delta) / se
        p2 = stats.t.cdf(t2, df=n-1)  # P(T <= t2)

        tost_p = max(p1, p2)
        equiv = tost_p < 0.05

        tcrit = stats.t.ppf(0.975, df=n-1)
        ci_lo = d_bar - tcrit * se
        ci_hi = d_bar + tcrit * se
        ci_within = (ci_lo > -delta) and (ci_hi < delta)

        print(f"  δ={delta}: t1={t1:.2f} p1={p1:.6f}, t2={t2:.2f} p2={p2:.6f}, "
              f"TOST p={tost_p:.6f} -> {'EQUIV' if equiv else 'NOT'}, "
              f"CI=[{ci_lo:+.4f}, {ci_hi:+.4f}] within=[{-delta},{+delta}]={ci_within}")

        results.append({
            'Basin': bname, 'Margin': delta_name, 'Delta': delta,
            'Mean_AE_diff': round(d_bar, 4), 'SE': round(se, 4),
            't1': round(t1, 2), 'p1': round(p1, 6),
            't2': round(t2, 2), 'p2': round(p2, 6),
            'TOST_p': round(tost_p, 6), 'Equivalent': 'YES' if equiv else 'NO',
            'CI_lo': round(ci_lo, 4), 'CI_hi': round(ci_hi, 4),
            'CI_within': 'YES' if ci_within else 'NO', 'N': n,
        })

pd.DataFrame(results).to_csv(os.path.join(SCRIPT_DIR, 'revision_tost_final.csv'), index=False)
print("\nDone.")
