"""
================================================================================
CORRECTED CANYON-B Comparison Script
================================================================================

Fixes all errors in the previous canyonb.py and canyonb2.py scripts:

  ERROR 1 (canyonb2.py): Fixed date datetime(2020,1,1) for ALL observations
           → FIXED: uses actual decimal year per observation from G2year/G2month

  ERROR 2 (canyonb2.py): No temporal split — ran on entire GLODAP dataset
           → FIXED: restricts to post-2018 external holdout only

  ERROR 3 (canyonb2.py): No basin separation — one global RMSE
           → FIXED: separate results per basin matching paper Table 3

  ERROR 4 (canyonb2.py): No QC filter — included Arctic estuarine obs (S≈0.1)
           → FIXED: same QC as all other scripts (S>25, T>-2.5, TCO2>1700, AOU>-50)

  ERROR 5 (canyonb2.py): Used G2pressure (many NaNs) instead of G2depth
           → FIXED: uses G2depth, converts to pressure (dbar ≈ depth in m)

  ERROR 6 (canyonb.py): CANYONB import failed silently → placeholder output
           → FIXED: uses canyonbpy package (the one you have installed)

PACKAGE REQUIRED:
  pip install canyonbpy
  (This is the package used in canyonb2.py — from PyPI, not GitHub)

OUTPUT:
  table3_corrected.csv   — basin RMSE values to fill Table 3
  table3_corrected.png   — bar chart for the paper
================================================================================
"""

import numpy as np
import scipy.io as sio
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ── Try importing CANYON-B (canyonbpy package) ────────────────────────────────
try:
    from canyonbpy import canyonb as CANYONB
    CANYONB_AVAILABLE = True
    print("✓ canyonbpy imported successfully")
except ImportError:
    CANYONB_AVAILABLE = False
    print("✗ canyonbpy not found. Install with: pip install canyonbpy")

# ── Config ────────────────────────────────────────────────────────────────────
EXTERNAL = 2018   # post-2018 = external holdout
TEMPORAL = 2015   # pre-2015 = training
SO_LAT   = -35

BASIN_CFG = {
    'Atlantic':      {'codes': [1],    'lat_min': SO_LAT},
    'Indian':        {'codes': [16],   'lat_min': SO_LAT},   # code 16 only (verified)
    'Pacific':       {'codes': [8],    'lat_min': SO_LAT},   # code 8 only (verified)
    'Southern Ocean':{'codes': None,   'lat_max': SO_LAT},
}

# ── Model definitions ─────────────────────────────────────────────────────────
def quadratic_model(X, alpha, beta, gamma, delta, epsilon):
    S, T, A = X
    core = (alpha * S) + (beta * T) + (gamma * (A / 100.0)) + delta
    return core**2 + epsilon

def mlr_feats(S, T, A):
    a = A / 100.0
    return np.column_stack([S, T, a, S**2, T**2, a**2, S*T, S*a, T*a])

# ── Load GLODAP ───────────────────────────────────────────────────────────────
print("Loading GLODAP data...")
mat = None
for path in ('GLODAPv2.2023_Merged_Master_File.mat',
             'data/raw/GLODAPv2.2023_Merged_Master_File.mat'):
    try:
        mat = sio.loadmat(path, squeeze_me=True)
        print(f"  Loaded: {path}")
        break
    except FileNotFoundError:
        continue

if mat is None:
    raise FileNotFoundError("GLODAP .mat not found.")

S_all   = mat['G2salinity'].astype(float)
T_all   = mat['G2temperature'].astype(float)
A_all   = mat['G2aou'].astype(float)
Y_all   = mat['G2tco2'].astype(float)
O2_all  = mat['G2oxygen'].astype(float)       # dissolved O2 for CANYON-B
Reg_all = mat['G2region'].astype(float)
Lat_all = mat['G2latitude'].astype(float)
Lon_all = mat['G2longitude'].astype(float)
Dep_all = mat['G2depth'].astype(float)        # FIX 5: depth not pressure
Yr_all  = mat['G2year'].astype(float)

# FIX 1: Build decimal year per observation using actual month
try:
    Mo_all = mat['G2month'].astype(float)
    DecYr_all = Yr_all + (Mo_all - 0.5) / 12.0
    print("  Decimal year built from G2year + G2month")
except KeyError:
    DecYr_all = Yr_all + 0.5
    print("  ⚠ G2month not found — using mid-year approximation")

# FIX 1 continued: build datetime objects per observation for canyonbpy
def decyr_to_datetime(decyr):
    """Convert decimal year array to list of datetime objects."""
    dts = []
    for dy in decyr:
        year = int(dy)
        remainder = dy - year
        # days in year
        import calendar
        days_in_year = 366 if calendar.isleap(year) else 365
        day_of_year = int(remainder * days_in_year) + 1
        try:
            dt = datetime(year, 1, 1) + pd.Timedelta(days=day_of_year - 1)
            dts.append(dt)
        except Exception:
            dts.append(datetime(year, 7, 1))
    return dts

# ── QC mask ───────────────────────────────────────────────────────────────────
# FIX 4: same QC filter as all other scripts
base_qc = (np.isfinite(S_all) & np.isfinite(T_all) & np.isfinite(A_all)
           & np.isfinite(Y_all) & np.isfinite(O2_all)
           & np.isfinite(Dep_all) & np.isfinite(Yr_all)
           & (S_all > 25) & (S_all < 42)
           & (T_all > -2.5) & (T_all < 35)
           & (Y_all > 1700) & (Y_all < 2600)
           & (A_all > -50)
           & (O2_all > 0) & (O2_all < 600))

def basin_mask(name, cfg):
    if name == 'Southern Ocean':
        return base_qc & (Lat_all < cfg['lat_max'])
    return base_qc & np.isin(Reg_all, cfg['codes']) & (Lat_all > cfg['lat_min'])

# ── Main loop ─────────────────────────────────────────────────────────────────
results = []

for bname, cfg in BASIN_CFG.items():
    print(f"\n{'='*60}")
    print(f"  Basin: {bname}")

    m   = basin_mask(bname, cfg)
    S   = S_all[m];    T   = T_all[m];    A   = A_all[m]
    Y   = Y_all[m];    O2  = O2_all[m];   dp  = Dep_all[m]
    yr  = Yr_all[m];   lat = Lat_all[m];  lon = Lon_all[m]
    dyr = DecYr_all[m]

    # FIX 2: restrict to correct temporal splits
    tr = yr < TEMPORAL
    ex = yr >= EXTERNAL      # post-2018 holdout only

    print(f"  Train (pre-{TEMPORAL}): {tr.sum():,}  |  Holdout (post-{EXTERNAL}): {ex.sum():,}")

    if ex.sum() < 50:
        print(f"  ⚠ Too few holdout obs — skipping")
        continue

    # Fit quadratic on training data
    popt, _ = curve_fit(quadratic_model,
                        (S[tr], T[tr], A[tr]), Y[tr],
                        p0=[1.0, -0.3, 2.5, 20.0, 2000.0], maxfev=30000)
    Yp_quad = quadratic_model((S[ex], T[ex], A[ex]), *popt)
    rmse_quad = np.sqrt(mean_squared_error(Y[ex], Yp_quad))
    r2_quad   = r2_score(Y[ex], Yp_quad)
    bias_quad = float(np.mean(Y[ex] - Yp_quad))
    print(f"  Quadratic — R²={r2_quad:.4f}  RMSE={rmse_quad:.3f}  Bias={bias_quad:.3f}")

    # MLR+AOU benchmark
    mlr = LinearRegression().fit(mlr_feats(S[tr], T[tr], A[tr]), Y[tr])
    Yp_mlr  = mlr.predict(mlr_feats(S[ex], T[ex], A[ex]))
    rmse_mlr = np.sqrt(mean_squared_error(Y[ex], Yp_mlr))
    r2_mlr   = r2_score(Y[ex], Yp_mlr)
    bias_mlr = float(np.mean(Y[ex] - Yp_mlr))
    print(f"  MLR+AOU   — R²={r2_mlr:.4f}  RMSE={rmse_mlr:.3f}  Bias={bias_mlr:.3f}")

    # CANYON-B on holdout
    rmse_cb = np.nan; r2_cb = np.nan; bias_cb = np.nan; n_cb = 0

    if CANYONB_AVAILABLE:
        print(f"  Running CANYON-B on {ex.sum():,} holdout observations...")
        try:
            # FIX 1: use actual per-observation datetimes (NOT a fixed date)
            dates_ex = decyr_to_datetime(dyr[ex])

            # FIX 5: use depth (dbar ≈ depth in m for most purposes)
            cb_out = CANYONB(
                gtime = dates_ex,
                lat   = lat[ex],
                lon   = lon[ex],
                pres  = dp[ex],       # depth in m ≈ pressure in dbar
                temp  = T[ex],
                psal  = S[ex],
                doxy  = O2[ex]        # FIX: dissolved O2, NOT AOU
            )

            Yp_cb = cb_out['CT']
            valid = np.isfinite(Yp_cb) & np.isfinite(Y[ex])
            n_cb  = valid.sum()
            print(f"  CANYON-B valid: {n_cb:,} / {ex.sum():,}")

            rmse_cb = np.sqrt(mean_squared_error(Y[ex][valid], Yp_cb[valid]))
            r2_cb   = r2_score(Y[ex][valid], Yp_cb[valid])
            bias_cb = float(np.mean(Y[ex][valid] - Yp_cb[valid]))
            print(f"  CANYON-B  — R²={r2_cb:.4f}  RMSE={rmse_cb:.3f}  Bias={bias_cb:.3f}")

        except Exception as e:
            print(f"  ✗ CANYON-B error: {e}")
    else:
        print("  CANYON-B not available")

    results.append({
        'Basin':       bname,
        'N_holdout':   int(ex.sum()),
        'Quad_R2':     round(r2_quad, 4),
        'Quad_RMSE':   round(rmse_quad, 3),
        'Quad_Bias':   round(bias_quad, 3),
        'MLR_R2':      round(r2_mlr, 4),
        'MLR_RMSE':    round(rmse_mlr, 3),
        'MLR_Bias':    round(bias_mlr, 3),
        'CB_R2':       round(r2_cb, 4) if not np.isnan(r2_cb) else None,
        'CB_RMSE':     round(rmse_cb, 3) if not np.isnan(rmse_cb) else None,
        'CB_Bias':     round(bias_cb, 3) if not np.isnan(bias_cb) else None,
        'N_CB_valid':  int(n_cb),
    })

# ── Save CSV ──────────────────────────────────────────────────────────────────
df = pd.DataFrame(results)
df.to_csv('table3_corrected.csv', index=False)
print("\n" + "="*60)
print("TABLE 3 RESULTS:")
print(df[['Basin','N_holdout','Quad_RMSE','MLR_RMSE','CB_RMSE']].to_string(index=False))

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
basins = df['Basin'].tolist()
x = np.arange(len(basins)); w = 0.25

quad_rmse = df['Quad_RMSE'].tolist()
mlr_rmse  = df['MLR_RMSE'].tolist()
cb_rmse   = [v if v is not None else np.nan for v in df['CB_RMSE'].tolist()]

b1 = ax.bar(x-w, quad_rmse, w, label='Quadratic (this work)',
            color='#1f77b4', edgecolor='black', lw=0.8)
b2 = ax.bar(x,   mlr_rmse,  w, label='MLR+AOU (10-param)',
            color='#aec7e8', edgecolor='black', lw=0.8)
b3 = ax.bar(x+w, cb_rmse,   w, label='CANYON-B (Bittig et al. 2018)',
            color='#d62728', edgecolor='black', lw=0.8, hatch='//')

for bar, val in zip(b1, quad_rmse):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.2,
            f'{val:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
for bar, val in zip(b2, mlr_rmse):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.2,
            f'{val:.1f}', ha='center', va='bottom', fontsize=9)
for bar, val in zip(b3, cb_rmse):
    if not np.isnan(val):
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.2,
                f'{val:.1f}', ha='center', va='bottom', fontsize=9, color='#d62728')

ax.set_xticks(x); ax.set_xticklabels(basins, fontsize=11)
ax.set_ylabel('RMSE (μmol kg⁻¹)', fontsize=12)
ax.set_title('Three-Way RMSE Comparison: Post-2018 External Holdout\n'
             'Quadratic Model  |  MLR+AOU Benchmark  |  CANYON-B',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=10); ax.grid(axis='y', linestyle='--', alpha=0.4)

plt.tight_layout()
fig.savefig('table3_corrected.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print("✓ Saved: table3_corrected.csv  and  table3_corrected.png")
print()
print("NOTE FOR METHODS SECTION:")
print("  CANYON-B (canyonbpy implementation of Bittig et al. 2018) was applied")
print("  to the post-2018 external holdout using dissolved oxygen (G2oxygen),")
print("  in-situ temperature, practical salinity, depth, latitude, longitude,")
print("  and actual observation date as inputs. Published weights were used")
print("  without retraining. Basin definitions match those of the quadratic model.")