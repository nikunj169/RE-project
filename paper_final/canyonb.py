"""
================================================================================
CANYON-B vs Quadratic Model vs MLR+AOU — Fill Table 3
================================================================================

PURPOSE:
  Run CANYON-B on the post-2018 external holdout (same data used for your
  quadratic model external RMSE) and compute a fair three-way comparison:
    - Quadratic model (your model)
    - MLR+AOU benchmark (10-parameter OLS)
    - CANYON-B (Bittig et al. 2018)

OUTPUT:
  table3_canyonb_comparison.csv   — numbers to paste into Table 3 of the paper
  table3_canyonb_comparison.png   — bar chart ready to insert as a figure

CANYON-B INSTALLATION:
  pip install canyonb
  OR clone from: https://github.com/HCBScienceProducts/CANYON-B
  The Python port is: https://github.com/BottomRedoxModel/canyonb
  Install with: pip install git+https://github.com/BottomRedoxModel/canyonb.git

  CANYON-B inputs required per observation:
    - date (decimal year)
    - latitude (degrees N)
    - longitude (degrees E)
    - depth (m)
    - temperature (deg C, in-situ)
    - salinity (practical)
    - oxygen (umol/kg)     ← this is O2, NOT AOU

  CANYON-B predicts: AT, CT (=TCO2), pH, pCO2, NO3, PO4, SiOH4
  We use the CT output.

NOTES:
  - CANYON-B was trained on GLODAP data. Using post-2018 GLODAP as holdout
    means some overlap is possible — flag this in the paper.
  - Report the CANYON-B version used (important for reproducibility).
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
import warnings
warnings.filterwarnings('ignore')

# ── Try importing CANYON-B ────────────────────────────────────────────────────
try:
    from canyonb import CANYONB
    CANYONB_AVAILABLE = True
    print("✓ CANYON-B imported successfully")
except ImportError:
    CANYONB_AVAILABLE = False
    print("✗ CANYON-B not found. Install with:")
    print("  pip install git+https://github.com/BottomRedoxModel/canyonb.git")
    print("  or: pip install canyonb")
    print("\nScript will run but CANYON-B column will be NaN.")
    print("Fill those cells manually once you have CANYON-B installed.\n")

# ── Config ────────────────────────────────────────────────────────────────────
EXTERNAL = 2018          # post-2018 is the holdout
TEMPORAL = 2015          # pre-2015 is training
SO_LAT   = -35

BASIN_CFG = {
    'Atlantic':      {'codes': [1],    'lat_min': SO_LAT},
    'Indian':        {'codes': [3,16], 'lat_min': SO_LAT},
    'Pacific':       {'codes': [2, 8], 'lat_min': SO_LAT},
    'Southern Ocean':{'codes': None,   'lat_max': SO_LAT},
}

# Your quadratic model external RMSE from your existing results
# (these are already known — from your final_summary_for_paper.txt)
QUAD_KNOWN = {
    'Atlantic':       {'r2_ext': 0.9210, 'rmse_ext': 17.68},
    'Indian':         {'r2_ext': 0.9682, 'rmse_ext': 18.44},
    'Pacific':        {'r2_ext': 0.9899, 'rmse_ext': 13.82},
    'Southern Ocean': {'r2_ext': 0.9651, 'rmse_ext':  8.99},
}

# ── Model definitions ─────────────────────────────────────────────────────────
def quadratic_model(X, alpha, beta, gamma, delta, epsilon):
    S, T, A = X
    core = (alpha * S) + (beta * T) + (gamma * (A / 100.0)) + delta
    return core**2 + epsilon

def mlr_feats(S, T, A):
    """10-parameter MLR+AOU feature matrix."""
    a = A / 100.0
    return np.column_stack([
        S, T, a,
        S**2, T**2, a**2,
        S*T, S*a, T*a
    ])

# ── Load GLODAP ───────────────────────────────────────────────────────────────
print("Loading GLODAP data...")
mat = None
for path in ('GLODAPv2.2023_Merged_Master_File.mat',
             'data/raw/GLODAPv2.2023_Merged_Master_File.mat'):
    try:
        mat = sio.loadmat(path, squeeze_me=True)
        print(f"  Loaded from: {path}")
        break
    except FileNotFoundError:
        continue

if mat is None:
    raise FileNotFoundError(
        "GLODAP .mat file not found. Place it in the working directory "
        "or data/raw/ subdirectory."
    )

# Extract variables
S_all   = mat['G2salinity'].astype(float)
T_all   = mat['G2temperature'].astype(float)
A_all   = mat['G2aou'].astype(float)
Y_all   = mat['G2tco2'].astype(float)
Reg_all = mat['G2region'].astype(float)
Lat_all = mat['G2latitude'].astype(float)
Lon_all = mat['G2longitude'].astype(float)
Dep_all = mat['G2depth'].astype(float)
Yr_all  = mat['G2year'].astype(float)

# Try to get dissolved O2 for CANYON-B (not AOU)
try:
    O2_all = mat['G2oxygen'].astype(float)
    has_o2 = True
    print("  Dissolved O2 (G2oxygen) available for CANYON-B")
except KeyError:
    has_o2 = False
    print("  ⚠ G2oxygen not found — will derive O2 from AOU + sat")
    print("    (CANYON-B uses O2, not AOU)")

# Try to get decimal year (for CANYON-B date input)
try:
    Month_all = mat['G2month'].astype(float)
    Day_all   = mat['G2day'].astype(float)
    DecYr_all = Yr_all + (Month_all - 0.5) / 12.0
except KeyError:
    DecYr_all = Yr_all + 0.5   # mid-year approximation
    print("  ⚠ Month/day not found — using mid-year approximation for CANYON-B")

# ── QC function ───────────────────────────────────────────────────────────────
def qc_mask(extra_mask):
    return (extra_mask
            & np.isfinite(S_all) & np.isfinite(T_all) & np.isfinite(A_all)
            & np.isfinite(Y_all) & np.isfinite(Dep_all) & np.isfinite(Yr_all)
            & (S_all > 25) & (S_all < 42)
            & (T_all > -2.5) & (T_all < 35)
            & (Y_all > 1700) & (Y_all < 2600)
            & (A_all > -50))

def basin_mask(name, cfg):
    if name == 'Southern Ocean':
        return Lat_all < cfg['lat_max']
    return np.isin(Reg_all, cfg['codes']) & (Lat_all > cfg['lat_min'])

# ── Main loop ─────────────────────────────────────────────────────────────────
results = []

for bname, cfg in BASIN_CFG.items():
    print(f"\n{'='*60}")
    print(f"  Processing: {bname}")
    print(f"{'='*60}")

    m = qc_mask(basin_mask(bname, cfg))
    S  = S_all[m];   T  = T_all[m];   A  = A_all[m]
    Y  = Y_all[m];   dp = Dep_all[m]; yr = Yr_all[m]
    lat= Lat_all[m]; lon= Lon_all[m]
    dyr= DecYr_all[m]

    if has_o2:
        o2 = O2_all[m]
    else:
        # Approximate: O2 = O2_sat - AOU
        # O2_sat from Garcia & Gordon 1992 (simplified)
        Ts = np.log((298.15 - T) / (273.15 + T))
        lnO2sat = (2.00907 + 3.22014*Ts + 4.05010*Ts**2
                   + 4.94457*Ts**3 - 0.256847*Ts**4 + 3.88767*Ts**5
                   + S*(-0.00624523 - 0.00737614*Ts - 0.010341*Ts**2
                        - 0.00817083*Ts**3)
                   - 4.88682e-7 * S**2)
        o2sat = np.exp(lnO2sat) * 44.614   # convert to umol/kg
        o2 = o2sat - A

    # Split into train / holdout
    tr = yr < TEMPORAL
    ex = yr >= EXTERNAL

    n_train = tr.sum()
    n_ext   = ex.sum()
    print(f"  Train (pre-{TEMPORAL}): {n_train:,}  |  Holdout (post-{EXTERNAL}): {n_ext:,}")

    if n_ext < 50:
        print(f"  ⚠ Too few external observations ({n_ext}) — skipping basin")
        continue

    # ── Fit quadratic model on training data ──────────────────────────────────
    p0 = [1.0, -0.3, 2.5, 20.0, 2000.0]
    popt, _ = curve_fit(quadratic_model,
                        (S[tr], T[tr], A[tr]), Y[tr],
                        p0=p0, maxfev=20000)

    Yp_quad_ext = quadratic_model((S[ex], T[ex], A[ex]), *popt)
    rmse_quad   = np.sqrt(mean_squared_error(Y[ex], Yp_quad_ext))
    r2_quad     = r2_score(Y[ex], Yp_quad_ext)
    bias_quad   = np.mean(Y[ex] - Yp_quad_ext)
    print(f"  Quadratic — R²={r2_quad:.4f}  RMSE={rmse_quad:.3f}  Bias={bias_quad:.3f}")

    # Verify against known values (sanity check)
    known_rmse = QUAD_KNOWN[bname]['rmse_ext']
    if abs(rmse_quad - known_rmse) > 1.0:
        print(f"  ⚠ RMSE mismatch: computed={rmse_quad:.2f}, "
              f"expected≈{known_rmse:.2f} — check QC or basin definition")
    else:
        print(f"  ✓ RMSE consistent with known value ({known_rmse:.2f})")

    # ── Fit MLR+AOU benchmark on training data ────────────────────────────────
    mlr = LinearRegression().fit(
        mlr_feats(S[tr], T[tr], A[tr]), Y[tr]
    )
    Yp_mlr_ext = mlr.predict(mlr_feats(S[ex], T[ex], A[ex]))
    rmse_mlr   = np.sqrt(mean_squared_error(Y[ex], Yp_mlr_ext))
    r2_mlr     = r2_score(Y[ex], Yp_mlr_ext)
    bias_mlr   = np.mean(Y[ex] - Yp_mlr_ext)
    print(f"  MLR+AOU  — R²={r2_mlr:.4f}  RMSE={rmse_mlr:.3f}  Bias={bias_mlr:.3f}")

    # ── CANYON-B on holdout ───────────────────────────────────────────────────
    rmse_cb = np.nan;  r2_cb = np.nan;  bias_cb = np.nan
    n_cb_valid = 0

    if CANYONB_AVAILABLE:
        print(f"  Running CANYON-B on {n_ext:,} holdout observations...")
        try:
            # CANYON-B call — adjust argument names to match your installed version
            # Standard Python port uses: CANYONB(date, lat, lon, pres, temp, psal, doxy)
            # where pres is in dbar (≈ depth in m for surface waters)
            cb_out = CANYONB(
                date = dyr[ex],
                lat  = lat[ex],
                lon  = lon[ex],
                pres = dp[ex],          # depth in m ≈ pressure in dbar
                temp = T[ex],
                psal = S[ex],
                doxy = o2[ex]
            )

            # CANYON-B returns a dict; TCO2 / CT key varies by version
            for key in ['CT', 'TCO2', 'tco2', 'DIC', 'dic']:
                if key in cb_out:
                    Yp_cb = cb_out[key]
                    break
            else:
                raise KeyError(f"TCO2 key not found in CANYON-B output. "
                               f"Available keys: {list(cb_out.keys())}")

            # Remove NaN predictions
            valid = np.isfinite(Yp_cb) & np.isfinite(Y[ex])
            n_cb_valid = valid.sum()
            print(f"  CANYON-B valid predictions: {n_cb_valid:,} / {n_ext:,}")

            rmse_cb = np.sqrt(mean_squared_error(Y[ex][valid], Yp_cb[valid]))
            r2_cb   = r2_score(Y[ex][valid], Yp_cb[valid])
            bias_cb = np.mean(Y[ex][valid] - Yp_cb[valid])
            print(f"  CANYON-B — R²={r2_cb:.4f}  RMSE={rmse_cb:.3f}  "
                  f"Bias={bias_cb:.3f}  (n={n_cb_valid:,})")

        except Exception as e:
            print(f"  ✗ CANYON-B failed: {e}")
            print("  Fill CANYON-B column manually in output CSV.")
    else:
        print("  CANYON-B not available — skipping (fill manually)")

    # ── Improvement calculations ──────────────────────────────────────────────
    quad_vs_mlr_pct  = (rmse_mlr - rmse_quad) / rmse_mlr * 100
    quad_vs_cb_pct   = (rmse_cb  - rmse_quad) / rmse_cb  * 100 \
                       if not np.isnan(rmse_cb) else np.nan

    results.append({
        'Basin':             bname,
        'N_holdout':         int(n_ext),
        'Quad_R2':           round(r2_quad, 4),
        'Quad_RMSE':         round(rmse_quad, 3),
        'Quad_Bias':         round(bias_quad, 3),
        'MLR_R2':            round(r2_mlr, 4),
        'MLR_RMSE':          round(rmse_mlr, 3),
        'MLR_Bias':          round(bias_mlr, 3),
        'CANYONB_R2':        round(r2_cb, 4) if not np.isnan(r2_cb) else '[run CANYON-B]',
        'CANYONB_RMSE':      round(rmse_cb, 3) if not np.isnan(rmse_cb) else '[run CANYON-B]',
        'CANYONB_Bias':      round(bias_cb, 3) if not np.isnan(bias_cb) else '[run CANYON-B]',
        'N_CANYONB_valid':   int(n_cb_valid),
        'Quad_vs_MLR_pct':   round(quad_vs_mlr_pct, 2),
        'Quad_vs_CB_pct':    round(quad_vs_cb_pct, 2) if not np.isnan(quad_vs_cb_pct) else '[run CANYON-B]',
    })

# ── Save CSV ──────────────────────────────────────────────────────────────────
df = pd.DataFrame(results)
df.to_csv('table3_canyonb_comparison.csv', index=False)

print("\n" + "="*60)
print("RESULTS SUMMARY — Copy into LaTeX Table 3")
print("="*60)
print(df[['Basin','N_holdout','Quad_RMSE','MLR_RMSE','CANYONB_RMSE']].to_string(index=False))

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
basins = df['Basin'].tolist()
x = np.arange(len(basins))
w = 0.25

quad_rmse = df['Quad_RMSE'].tolist()
mlr_rmse  = df['MLR_RMSE'].tolist()
cb_rmse_raw = df['CANYONB_RMSE'].tolist()
cb_rmse = [v if isinstance(v, float) else np.nan for v in cb_rmse_raw]

b1 = ax.bar(x - w, quad_rmse, w, label='Quadratic Model (this work)',
            color='#1f77b4', edgecolor='black', linewidth=0.8)
b2 = ax.bar(x,     mlr_rmse,  w, label='MLR+AOU Benchmark (10-param)',
            color='#aec7e8', edgecolor='black', linewidth=0.8)
b3 = ax.bar(x + w, cb_rmse,   w, label='CANYON-B (Bittig et al. 2018)',
            color='#d62728', edgecolor='black', linewidth=0.8,
            hatch='//')

# Value labels
for rect, val in zip(b1, quad_rmse):
    ax.text(rect.get_x() + rect.get_width()/2, rect.get_height() + 0.3,
            f'{val:.1f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
for rect, val in zip(b2, mlr_rmse):
    ax.text(rect.get_x() + rect.get_width()/2, rect.get_height() + 0.3,
            f'{val:.1f}', ha='center', va='bottom', fontsize=9)
for rect, val in zip(b3, cb_rmse):
    if not np.isnan(val):
        ax.text(rect.get_x() + rect.get_width()/2, rect.get_height() + 0.3,
                f'{val:.1f}', ha='center', va='bottom', fontsize=9, color='#d62728')
    else:
        ax.text(rect.get_x() + rect.get_width()/2, 1.5,
                '[run\nCANYON-B]', ha='center', va='bottom',
                fontsize=7, color='gray', style='italic')

ax.set_xticks(x)
ax.set_xticklabels(basins, fontsize=11)
ax.set_ylabel('RMSE ($\\mu$mol kg$^{-1}$)', fontsize=12)
ax.set_title('Three-Way RMSE Comparison: Post-2018 External Holdout\n'
             'Quadratic Model  |  MLR+AOU Benchmark  |  CANYON-B',
             fontsize=12, fontweight='bold')
ax.legend(fontsize=10, loc='upper right')
ax.set_ylim(0, max([v for v in quad_rmse + mlr_rmse if v])*1.3)
ax.grid(axis='y', linestyle='--', alpha=0.4)
ax.set_axisbelow(True)

plt.tight_layout()
fig.savefig('table3_canyonb_comparison.png', dpi=300, bbox_inches='tight')
plt.close(fig)

print("\n✓ Saved: table3_canyonb_comparison.csv")
print("✓ Saved: table3_canyonb_comparison.png")

# ── LaTeX table snippet ───────────────────────────────────────────────────────
print("\n" + "="*60)
print("LATEX TABLE 3 SNIPPET — paste into your .tex file")
print("="*60)
print(r"""
\begin{table}[H]
\centering
\caption{Three-way RMSE comparison on the post-2018 external holdout.
CANYON-B uses the full input vector ($S$, $T$, $O_2$, depth, latitude,
longitude, date). All models trained/fitted on pre-2015 data only.}
\label{tab:canyonb}
\begin{tabular}{lcccccc}
\toprule
Basin & $N_\text{holdout}$ &
\multicolumn{2}{c}{Quadratic (this work)} &
\multicolumn{1}{c}{MLR+AOU} &
\multicolumn{1}{c}{CANYON-B} \\
\cmidrule(lr){3-4}
& & $R^2$ & RMSE & RMSE & RMSE \\
& & & ($\mu$mol~kg$^{-1}$) & ($\mu$mol~kg$^{-1}$) &
($\mu$mol~kg$^{-1}$) \\
\midrule""")

for row in results:
    cb_rmse_str = (f"{row['CANYONB_RMSE']:.2f}"
                   if isinstance(row['CANYONB_RMSE'], float)
                   else r'$[\ast]$')
    mlr_str = f"{row['MLR_RMSE']:.2f}"
    print(f"{row['Basin']:15s} & {row['N_holdout']:>6,} & "
          f"{row['Quad_R2']:.4f} & {row['Quad_RMSE']:.2f} & "
          f"{mlr_str} & {cb_rmse_str} \\\\")

print(r"""\bottomrule
\end{tabular}
\end{table}""")

print("\n" + "="*60)
print("NEXT STEPS:")
print("="*60)
print("""
1. If CANYON-B was not installed, install it:
      pip install git+https://github.com/BottomRedoxModel/canyonb.git

2. Re-run this script — CANYON-B RMSE values will fill automatically.

3. Copy the LaTeX snippet above into your paper replacing Table 3.

4. Note for your Methods section:
   'CANYON-B v[version] was applied to the post-2018 external holdout
   using dissolved oxygen (G2oxygen), in-situ temperature, practical
   salinity, depth, latitude, longitude, and decimal year as inputs.
   CANYON-B was not retrained on our dataset; we use the published
   weights to ensure a fair comparison of out-of-the-box performance.'

5. If CANYON-B RMSE is ~10-14 umol/kg (expected from literature),
   your quadratic model will show a 5-8 umol/kg gap — be transparent
   about this and use the 3 use-case framing already in your paper.
""")