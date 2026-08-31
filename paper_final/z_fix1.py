"""
================================================================================
FIX 1 — Corrected Bootstrap Prediction Intervals
================================================================================

PROBLEM: Previous bootstrap gave 1–2% coverage instead of 95%.
CAUSE:   Old code only propagated PARAMETER uncertainty. With 100k+ training
         points, parameters are pinned so tightly that the band collapsed to
         ~1 μmol/kg wide — missing all the irreducible residual noise.

SOLUTION: A proper 95% PREDICTION interval = parameter uncertainty + residual
          noise. After each bootstrap resample-and-refit, we add a draw from
          N(0, σ_resid) to each prediction. This gives the correct interval
          that should contain ~95% of future observations.

          PI_width² = parameter_uncertainty² + residual_variance²

OUTPUTS:
  fix1_prediction_intervals.png   — corrected PI plot (target: ≥90% coverage)
  fix1_pi_summary.csv             — coverage and width per basin
================================================================================
"""

import numpy as np
import scipy.io as sio
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error, r2_score
import warnings
warnings.filterwarnings('ignore')

# ── Config ────────────────────────────────────────────────────────────────────
SO_LAT   = -35
N_BOOT   = 300      # more resamples for stable intervals
TEMPORAL = 2015
EXTERNAL = 2018
PLOT_N   = 600      # points shown in figure (sorted subsample)

BASIN_CFG = {
    'Atlantic':      {'codes': [1],    'lat_min': SO_LAT},
    'Pacific':       {'codes': [2, 8], 'lat_min': SO_LAT},
    'Southern Ocean':{'codes': None,   'lat_max': SO_LAT},
}
BASIN_COLORS = {'Atlantic': '#1f77b4', 'Pacific': '#ff7f0e', 'Southern Ocean': '#9467bd'}

# ── Model ─────────────────────────────────────────────────────────────────────
def quadratic_model(X, alpha, beta, gamma, delta, epsilon):
    S, T, A = X
    core = (alpha * S) + (beta * T) + (gamma * (A / 100.0)) + delta
    return core**2 + epsilon

# ── Data loading ──────────────────────────────────────────────────────────────
print("Loading GLODAP data...")
for path in ('GLODAPv2.2023_Merged_Master_File.mat',
             'data/raw/GLODAPv2.2023_Merged_Master_File.mat'):
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

def qc(mask):
    return (mask
            & np.isfinite(S_all) & np.isfinite(T_all) & np.isfinite(A_all)
            & np.isfinite(Y_all) & np.isfinite(Dep)   & np.isfinite(Yr)
            & (S_all > 25) & (S_all < 42)
            & (T_all > -2.5) & (T_all < 35)
            & (Y_all > 1700) & (Y_all < 2600)
            & (A_all > -50))

def basin_mask(name, cfg):
    if name == 'Southern Ocean':
        return Lat < cfg['lat_max']
    return np.isin(Reg, cfg['codes']) & (Lat > cfg['lat_min'])

# ── Main loop ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
summary = []

for ax, (bname, cfg) in zip(axes, BASIN_CFG.items()):
    print(f"\n  {bname}")
    m   = qc(basin_mask(bname, cfg))
    S,T,A,Y,dp,yr = S_all[m],T_all[m],A_all[m],Y_all[m],Dep[m],Yr[m]

    tr = yr < TEMPORAL
    te = (yr >= TEMPORAL) & (yr < EXTERNAL)
    if te.sum() < 100:
        print(f"    ⚠ Only {te.sum()} test samples — skipping")
        ax.set_visible(False); continue

    Xtr = (S[tr], T[tr], A[tr]);  Ytr = Y[tr]
    Xte = (S[te], T[te], A[te]);  Yte = Y[te]
    n_tr = Ytr.shape[0]

    # ── OLS fit on full training set ──────────────────────────────────────────
    p0 = [1.0, -0.3, 2.5, 20.0, 2000.0]
    popt, _ = curve_fit(quadratic_model, Xtr, Ytr, p0=p0, maxfev=20000)

    Ypred_te = quadratic_model(Xte, *popt)

    # ── CORRECT residual σ from training set ─────────────────────────────────
    train_resid = Ytr - quadratic_model(Xtr, *popt)
    sigma_resid = np.std(train_resid)          # irreducible noise level
    print(f"    σ_resid (train) = {sigma_resid:.3f} μmol/kg")

    # ── Bootstrap: resample → refit → predict → ADD residual noise ───────────
    n_te = Yte.shape[0]
    boot_preds = np.full((N_BOOT, n_te), np.nan)
    ok = 0

    for b in range(N_BOOT):
        idx = np.random.choice(n_tr, n_tr, replace=True)
        Xb  = tuple(x[idx] for x in Xtr)
        Yb  = Ytr[idx]
        try:
            pb, _ = curve_fit(quadratic_model, Xb, Yb, p0=popt, maxfev=8000)
            # ── KEY FIX: add irreducible noise draw to each prediction ────────
            noise = np.random.normal(0, sigma_resid, n_te)
            boot_preds[b] = quadratic_model(Xte, *pb) + noise
            ok += 1
        except RuntimeError:
            pass

    print(f"    Bootstrap: {ok}/{N_BOOT} successful")

    valid = boot_preds[~np.any(np.isnan(boot_preds), axis=1)]
    ci_lo = np.percentile(valid,  2.5, axis=0)
    ci_hi = np.percentile(valid, 97.5, axis=0)

    coverage     = np.mean((Yte >= ci_lo) & (Yte <= ci_hi)) * 100
    mean_width   = np.mean(ci_hi - ci_lo)
    rmse_te      = np.sqrt(mean_squared_error(Yte, Ypred_te))
    r2_te        = r2_score(Yte, Ypred_te)

    print(f"    Coverage: {coverage:.1f}%  |  Mean width: {mean_width:.1f} μmol/kg")
    print(f"    R²={r2_te:.4f}  RMSE={rmse_te:.2f}")

    summary.append({
        'Basin': bname,
        'R2_test': round(r2_te, 4),
        'RMSE_test': round(rmse_te, 3),
        'sigma_resid': round(float(sigma_resid), 3),
        'N_boot_ok': ok,
        'PI_coverage_pct': round(coverage, 1),
        'PI_mean_width': round(mean_width, 2),
        'Alpha': round(float(popt[0]), 4),
        'Beta':  round(float(popt[1]), 4),
        'Gamma': round(float(popt[2]), 4),
        'Delta': round(float(popt[3]), 4),
        'Epsilon': round(float(popt[4]), 2),
    })

    # ── Plot ──────────────────────────────────────────────────────────────────
    idx_sort = np.argsort(Yte)[:PLOT_N]
    Ys  = Yte[idx_sort];  Yps = Ypred_te[idx_sort]
    los = ci_lo[idx_sort]; his = ci_hi[idx_sort]
    xs  = np.arange(len(Ys))
    col = BASIN_COLORS[bname]

    ax.fill_between(xs, los, his, alpha=0.35, color=col,
                    label=f'95% PI (n={ok} boots)')
    ax.plot(xs, Ys,  'k.', ms=3, alpha=0.5, label='Measured')
    ax.plot(xs, Yps, '-',  color=col, lw=1.5, alpha=0.9, label='Predicted (OLS)')

    cov_color = 'darkgreen' if abs(coverage - 95) < 5 else \
                'darkorange' if abs(coverage - 95) < 10 else 'darkred'
    ax.text(0.97, 0.04,
            f'Coverage: {coverage:.1f}%\n(target: 95%)\nWidth: {mean_width:.1f} μmol/kg',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=10,
            color=cov_color, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.9, edgecolor=cov_color, lw=1.5))

    ax.set_title(f'{bname}\n$R^2$={r2_te:.4f}, RMSE={rmse_te:.2f} μmol kg$^{{-1}}$',
                 fontweight='bold')
    ax.set_xlabel('Sample index (sorted by measured TCO$_2$)')
    ax.set_ylabel('TCO$_2$ (μmol kg$^{-1}$)')
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.3, linestyle=':')

fig.suptitle(
    'FIX 1 — Corrected 95% Bootstrap Prediction Intervals\n'
    'PI = parameter uncertainty + residual noise  '
    '(correct method: boot_pred + N(0, σ_resid))',
    fontsize=13, fontweight='bold'
)
plt.tight_layout()
fig.savefig('fix1_prediction_intervals.png', dpi=300, bbox_inches='tight')
plt.close(fig)

df = pd.DataFrame(summary)
df.to_csv('fix1_pi_summary.csv', index=False)

print("\n" + "="*60)
print("FIX 1 COMPLETE")
print("  fix1_prediction_intervals.png")
print("  fix1_pi_summary.csv")
print("\nSummary:")
print(df[['Basin','PI_coverage_pct','PI_mean_width','R2_test','RMSE_test']].to_string(index=False))
print("="*60)
print("""
EXPLANATION FOR METHODS SECTION:
  95% prediction intervals were constructed via parametric bootstrap
  (n=300 resamples). Each bootstrap iteration resampled the training
  set with replacement, refitted the quadratic model, generated
  predictions on the test set, and added a draw from N(0, σ) where
  σ is the training-set residual standard deviation (capturing
  irreducible observation noise). The 2.5th and 97.5th percentiles
  of the resulting distribution form the reported prediction interval.
""")