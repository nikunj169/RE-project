import numpy as np
import scipy.io as sio
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.linear_model import LinearRegression
import warnings

# Ignore calculation warnings during fitting
warnings.filterwarnings('ignore')

# =============================================================================
# 1. LOAD DATA & PRE-PROCESSING
# =============================================================================
print("Loading GLODAP Data...")
filename = 'GLODAPv2.2023_Merged_Master_File.mat'

# Path handling for different environments
try:
    mat = sio.loadmat(filename, squeeze_me=True)
except FileNotFoundError:
    mat = sio.loadmat('data/raw/' + filename, squeeze_me=True)

# Variable extraction
S_all = mat['G2salinity']
T_all = mat['G2temperature']
AOU_all = mat['G2aou']
Y_all = mat['G2tco2']
Region = mat['G2region']
Lat = mat['G2latitude']
Depth = mat['G2depth']
Year = mat['G2year']

# =============================================================================
# 2. MODEL DEFINITIONS
# =============================================================================
def your_quadratic_model(X, alpha, beta, gamma, delta, epsilon):
    """The proposed Ocean Carbon Model (Symbolic Regression Derivative)"""
    S, T, A = X
    A_norm = A / 100.0
    core = (alpha * S) + (beta * T) + (gamma * A_norm) + delta
    # The quadratic form approximates the non-linear CO2 solubility/buffer system
    return (core**2) + epsilon

# =============================================================================
# 3. SETTINGS & INITIALIZATION
# =============================================================================
basin_map = {
    'Atlantic': [1], 
    'Pacific': [8, 2], 
    'Indian': [16, 3], 
    'Southern Ocean': 'Lat'
}
temporal_threshold = 2015
depth_bins = [0, 100, 500, 1000, 3000, 6000]
depth_labels = ['0-100m', '100-500m', '500-1000m', '1000-3000m', '>3000m']

final_results = []
depth_stats = []

print(f"\n{'Basin':<15} | {'Year Split':<12} | {'Your RMSE':<10} | {'Poly RMSE':<10} | {'Improvement'}")
print("-" * 80)

# =============================================================================
# 4. PROCESSING LOOP
# =============================================================================
for name, code in basin_map.items():
    # Region masking
    if name == 'Southern Ocean':
        mask = (Lat < -50)
    else:
        mask = np.isin(Region, code)
    
    # Cleaning data for valid entries only
    valid = mask & ~np.isnan(S_all) & ~np.isnan(T_all) & ~np.isnan(AOU_all) & ~np.isnan(Y_all) & ~np.isnan(Year) & ~np.isnan(Depth)
    S, T, A, Y, yr, dep = S_all[valid], T_all[valid], AOU_all[valid], Y_all[valid], Year[valid], Depth[valid]
    
    # Temporal Split (Train on Historical, Test on Recent)
    train_mask = yr < temporal_threshold
    test_mask = yr >= temporal_threshold
    
    # Safety check for sample size
    if np.sum(test_mask) < 50:
        continue

    X_train = (S[train_mask], T[train_mask], A[train_mask])
    Y_train = Y[train_mask]
    
    X_test = (S[test_mask], T[test_mask], A[test_mask])
    Y_test = Y[test_mask]
    Dep_test = dep[test_mask]

    # --- 1. FIT YOUR PROPOSED MODEL ---
    # FIX: Added maxfev=10000 to resolve Southern Ocean convergence error
    try:
        popt, _ = curve_fit(your_quadratic_model, X_train, Y_train, 
                            p0=[1.0, -0.3, 2.5, 20, 2000], 
                            maxfev=10000)
        
        Y_pred_your = your_quadratic_model(X_test, *popt)
        rmse_your = np.sqrt(mean_squared_error(Y_test, Y_pred_your))
    except RuntimeError:
        print(f"{name:<15} | Fit failed to converge.")
        continue

    # --- 2. FIT LITERATURE BENCHMARK (Polynomial S,T) ---
    # This represents standard physical empirical relationships
    poly_train = np.column_stack([S[train_mask], T[train_mask], S[train_mask]**2, T[train_mask]**2])
    poly_test = np.column_stack([S[test_mask], T[test_mask], S[test_mask]**2, T[test_mask]**2])
    
    poly_model = LinearRegression().fit(poly_train, Y_train)
    Y_pred_poly = poly_model.predict(poly_test)
    rmse_poly = np.sqrt(mean_squared_error(Y_test, Y_pred_poly))

    improvement = ((rmse_poly - rmse_your) / rmse_poly) * 100
    
    print(f"{name:<15} | <{temporal_threshold} vs >= | {rmse_your:<10.2f} | {rmse_poly:<10.2f} | {improvement:.1f}%")

    final_results.append({
        'Basin': name, 
        'Test_Year_Range': f'{temporal_threshold}-2023',
        'Your_RMSE': rmse_your, 
        'Literature_RMSE': rmse_poly, 
        'Pct_Better': improvement
    })

    # --- 3. DEPTH-STRATIFIED ANALYSIS ---
    for i in range(len(depth_bins)-1):
        d_mask = (Dep_test >= depth_bins[i]) & (Dep_test < depth_bins[i+1])
        if np.sum(d_mask) > 10:
            d_rmse = np.sqrt(mean_squared_error(Y_test[d_mask], Y_pred_your[d_mask]))
            depth_stats.append({
                'Basin': name, 
                'Layer': depth_labels[i], 
                'RMSE': d_rmse
            })

# =============================================================================
# 5. GENERATE PUBLICATION PLOTS & EXPORTS
# =============================================================================
df_depth = pd.DataFrame(depth_stats)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

# Plot A: Temporal Validation Summary
basins = [r['Basin'] for r in final_results]
y_rmse = [r['Your_RMSE'] for r in final_results]
l_rmse = [r['Literature_RMSE'] for r in final_results]
x = np.arange(len(basins))

ax1.bar(x - 0.2, l_rmse, 0.4, label='Literature Benchmark (Poly S,T)', color='lightgrey')
ax1.bar(x + 0.2, y_rmse, 0.4, label='Your Proposed Model', color='navy')
ax1.set_ylabel('RMSE (μmol/kg)', fontsize=12)
ax1.set_title('A. Temporal Stability (Predictive Power on 2015-2023 Data)', fontsize=14, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(basins)
ax1.legend()
ax1.grid(axis='y', linestyle='--', alpha=0.6)

# Plot B: Vertical Robustness
for basin in df_depth['Basin'].unique():
    subset = df_depth[df_depth['Basin'] == basin]
    ax2.plot(subset['Layer'], subset['RMSE'], marker='s', label=basin, linewidth=2.5)

ax2.set_ylabel('RMSE (μmol/kg)', fontsize=12)
ax2.set_title('B. Depth Robustness (Accuracy by Water Column Layer)', fontsize=14, fontweight='bold')
ax2.grid(linestyle='--', alpha=0.6)
ax2.legend()

plt.tight_layout()
plt.savefig('Final_Validation_Results_Corrected.png', dpi=300)

# Export Data for Manuscript Tables
pd.DataFrame(final_results).to_csv('results_temporal_summary.csv', index=False)
df_depth.pivot(index='Layer', columns='Basin', values='RMSE').to_csv('results_depth_breakdown.csv')

print("\n" + "="*60)
print("✓ VALIDATION COMPLETE")
print("1. Plot saved: Final_Validation_Results_Corrected.png")
print("2. Summary Table: results_temporal_summary.csv")
print("3. Depth Table: results_depth_breakdown.csv")
print("="*60)