"""
Residual Analysis and Diagnostic Plots
======================================

Creates publication-quality diagnostic plots to assess model performance:
1. Residual vs. Measured TCO2
2. Residual vs. Depth (if available)
3. Residual vs. Latitude
4. Residual distribution histogram
5. Q-Q plot for normality
6. Spatial residual map

Run AFTER fitting your model to identify systematic biases.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.optimize import curve_fit
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# MODEL DEFINITION
# =============================================================================
def stable_model(X, alpha, beta, gamma, delta, epsilon):
    S, T, A = X
    A_norm = A / 100.0
    core = (alpha * S) + (beta * T) + (gamma * A_norm) + delta
    return (core**2) + epsilon

# =============================================================================
# LOAD DATA
# =============================================================================
print("Loading data and fitting models...")
mat = sio.loadmat('data/raw/GLODAPv2.2023_Merged_Master_File.mat', squeeze_me=True)
S_all, T_all, AOU_all, Y_all = mat['G2salinity'], mat['G2temperature'], mat['G2aou'], mat['G2tco2']
Region, Lat = mat['G2region'], mat['G2latitude']

# Try to load additional variables for diagnostic plots
try:
    Depth_all = mat['G2depth']
    has_depth = True
except:
    has_depth = False
    print("Warning: Depth not available")

try:
    Lon_all = mat['G2longitude']
    has_lon = True
except:
    has_lon = False
    print("Warning: Longitude not available for spatial maps")

basin_map = {'Atlantic': [1], 'Pacific': [8, 2], 'Indian': [16, 3], 'Southern Ocean': 'Lat'}

# =============================================================================
# CREATE FIGURE
# =============================================================================
fig = plt.figure(figsize=(20, 14))

for idx, (basin_name, b_id) in enumerate(basin_map.items()):
    print(f"\nProcessing {basin_name}...")
    
    # Apply masks
    if basin_name == 'Southern Ocean':
        mask = (Lat < -30) & (Region != 4)
    else:
        mask = (np.isin(Region, b_id)) & (Lat > -30)
    
    mask &= (S_all > 30) & (Y_all > 1800) & (AOU_all > 0)
    S, T, A, Y = S_all[mask], T_all[mask], AOU_all[mask], Y_all[mask]
    Lat_basin = Lat[mask]
    
    if has_depth:
        Depth = Depth_all[mask]
    if has_lon:
        Lon_basin = Lon_all[mask]
    
    # Fit model
    popt, _ = curve_fit(stable_model, (S, T, A), Y, p0=[0.1, -1.0, 5.0, 10.0, 1500])
    Y_pred = stable_model((S, T, A), *popt)
    residuals = Y - Y_pred
    
    # Calculate statistics
    rmse = np.sqrt(np.mean(residuals**2))
    mae = np.mean(np.abs(residuals))
    bias = np.mean(residuals)
    
    # =========================================================================
    # PLOT 1: Residual vs. Measured TCO2
    # =========================================================================
    ax1 = plt.subplot(4, 5, idx*5 + 1)
    
    # Hexbin for density visualization
    hb = ax1.hexbin(Y, residuals, gridsize=50, cmap='YlOrRd', bins='log', mincnt=1)
    ax1.axhline(0, color='black', linestyle='--', linewidth=2, label='Zero Line')
    ax1.axhline(rmse, color='red', linestyle=':', linewidth=1.5, alpha=0.7, label=f'±RMSE ({rmse:.1f})')
    ax1.axhline(-rmse, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
    
    ax1.set_xlabel('Measured TCO$_2$ (μmol/kg)', fontsize=10)
    ax1.set_ylabel('Residual (μmol/kg)', fontsize=10)
    ax1.set_title(f'{basin_name}\nResidual vs. Measured', fontsize=11, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':')
    ax1.legend(fontsize=8, loc='upper left')
    
    # Add text box with statistics
    textstr = f'RMSE={rmse:.2f}\nBias={bias:.2f}\nMAE={mae:.2f}'
    ax1.text(0.98, 0.02, textstr, transform=ax1.transAxes, fontsize=9,
             verticalalignment='bottom', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # =========================================================================
    # PLOT 2: Residual vs. Depth
    # =========================================================================
    ax2 = plt.subplot(4, 5, idx*5 + 2)
    
    if has_depth:
        ax2.hexbin(Depth, residuals, gridsize=50, cmap='YlOrRd', bins='log', mincnt=1)
        ax2.axhline(0, color='black', linestyle='--', linewidth=2)
        ax2.axhline(rmse, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
        ax2.axhline(-rmse, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
        
        ax2.set_xlabel('Depth (m)', fontsize=10)
        ax2.set_ylabel('Residual (μmol/kg)', fontsize=10)
        ax2.set_title('Residual vs. Depth', fontsize=11, fontweight='bold')
        ax2.invert_xaxis()  # Deep ocean on right
        ax2.grid(True, alpha=0.3, linestyle=':')
        
        # Check for depth-dependent bias
        shallow_mask = Depth < 500
        deep_mask = Depth > 2000
        if np.sum(shallow_mask) > 100 and np.sum(deep_mask) > 100:
            shallow_bias = np.mean(residuals[shallow_mask])
            deep_bias = np.mean(residuals[deep_mask])
            bias_diff = abs(shallow_bias - deep_bias)
            
            if bias_diff > 5:
                ax2.text(0.5, 0.95, f'⚠️ Depth bias: {bias_diff:.1f}', 
                        transform=ax2.transAxes, fontsize=9, ha='center', va='top',
                        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    else:
        ax2.text(0.5, 0.5, 'Depth data\nnot available', 
                transform=ax2.transAxes, ha='center', va='center', fontsize=12)
        ax2.set_xticks([])
        ax2.set_yticks([])
    
    # =========================================================================
    # PLOT 3: Residual vs. Latitude
    # =========================================================================
    ax3 = plt.subplot(4, 5, idx*5 + 3)
    
    ax3.hexbin(Lat_basin, residuals, gridsize=50, cmap='YlOrRd', bins='log', mincnt=1)
    ax3.axhline(0, color='black', linestyle='--', linewidth=2)
    ax3.axhline(rmse, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
    ax3.axhline(-rmse, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
    
    ax3.set_xlabel('Latitude (°N)', fontsize=10)
    ax3.set_ylabel('Residual (μmol/kg)', fontsize=10)
    ax3.set_title('Residual vs. Latitude', fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3, linestyle=':')
    
    # Check for latitudinal bias
    if basin_name != 'Southern Ocean':
        tropical_mask = (Lat_basin > -30) & (Lat_basin < 30)
        high_lat_mask = (np.abs(Lat_basin) > 45)
        
        if np.sum(tropical_mask) > 100 and np.sum(high_lat_mask) > 100:
            trop_bias = np.mean(residuals[tropical_mask])
            high_bias = np.mean(residuals[high_lat_mask])
            lat_bias_diff = abs(trop_bias - high_bias)
            
            if lat_bias_diff > 5:
                ax3.text(0.5, 0.95, f'⚠️ Latitude bias: {lat_bias_diff:.1f}', 
                        transform=ax3.transAxes, fontsize=9, ha='center', va='top',
                        bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    # =========================================================================
    # PLOT 4: Residual Distribution
    # =========================================================================
    ax4 = plt.subplot(4, 5, idx*5 + 4)
    
    ax4.hist(residuals, bins=60, edgecolor='black', alpha=0.7, color='skyblue', density=True)
    
    # Overlay normal distribution
    mu, sigma = np.mean(residuals), np.std(residuals)
    x = np.linspace(residuals.min(), residuals.max(), 100)
    ax4.plot(x, stats.norm.pdf(x, mu, sigma), 'r-', linewidth=2, label='Normal fit')
    
    ax4.axvline(0, color='black', linestyle='--', linewidth=2, label='Zero')
    ax4.set_xlabel('Residual (μmol/kg)', fontsize=10)
    ax4.set_ylabel('Density', fontsize=10)
    ax4.set_title('Residual Distribution', fontsize=11, fontweight='bold')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3, linestyle=':')
    
    # Normality test
    _, p_value = stats.normaltest(residuals)
    if p_value < 0.01:
        ax4.text(0.02, 0.98, '⚠️ Non-normal\n(p<0.01)', 
                transform=ax4.transAxes, fontsize=9, va='top',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    else:
        ax4.text(0.02, 0.98, '✓ Normal', 
                transform=ax4.transAxes, fontsize=9, va='top',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # =========================================================================
    # PLOT 5: Q-Q Plot
    # =========================================================================
    ax5 = plt.subplot(4, 5, idx*5 + 5)
    
    stats.probplot(residuals, dist="norm", plot=ax5)
    ax5.set_title('Q-Q Plot', fontsize=11, fontweight='bold')
    ax5.grid(True, alpha=0.3, linestyle=':')
    
    # Add reference line
    ax5.get_lines()[0].set_markerfacecolor('skyblue')
    ax5.get_lines()[0].set_markersize(3)
    ax5.get_lines()[0].set_alpha(0.6)

plt.tight_layout()
plt.savefig('Residual_Analysis_Comprehensive.png', dpi=300, bbox_inches='tight')
print("\n✓ Residual analysis plots saved to 'Residual_Analysis_Comprehensive.png'")

# =============================================================================
# SPATIAL RESIDUAL MAP (if longitude available)
# =============================================================================
if has_lon:
    print("\nCreating spatial residual map...")
    
    fig_map, axes_map = plt.subplots(2, 2, figsize=(18, 10))
    axes_map = axes_map.flatten()
    
    for idx, (basin_name, b_id) in enumerate(basin_map.items()):
        # Apply masks
        if basin_name == 'Southern Ocean':
            mask = (Lat < -30) & (Region != 4)
        else:
            mask = (np.isin(Region, b_id)) & (Lat > -30)
        
        mask &= (S_all > 30) & (Y_all > 1800) & (AOU_all > 0)
        S, T, A, Y = S_all[mask], T_all[mask], AOU_all[mask], Y_all[mask]
        Lat_basin = Lat[mask]
        Lon_basin = Lon_all[mask]
        
        # Fit and calculate residuals
        popt, _ = curve_fit(stable_model, (S, T, A), Y, p0=[0.1, -1.0, 5.0, 10.0, 1500])
        Y_pred = stable_model((S, T, A), *popt)
        residuals = Y - Y_pred
        
        # Plot spatial distribution
        ax = axes_map[idx]
        scatter = ax.scatter(Lon_basin, Lat_basin, c=residuals, s=2, alpha=0.6,
                           cmap='RdBu_r', vmin=-50, vmax=50)
        
        ax.set_xlabel('Longitude (°E)', fontsize=11)
        ax.set_ylabel('Latitude (°N)', fontsize=11)
        ax.set_title(f'{basin_name} - Spatial Residuals', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle=':')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Residual (μmol/kg)', fontsize=10)
        
        # Add RMSE annotation
        rmse = np.sqrt(np.mean(residuals**2))
        ax.text(0.02, 0.98, f'RMSE = {rmse:.2f}', transform=ax.transAxes,
               fontsize=10, va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('Spatial_Residual_Map.png', dpi=300, bbox_inches='tight')
    print("✓ Spatial residual map saved to 'Spatial_Residual_Map.png'")

# =============================================================================
# INTERPRETATION GUIDE
# =============================================================================
print("\n" + "="*80)
print("RESIDUAL ANALYSIS INTERPRETATION GUIDE")
print("="*80)

print("""
What to Look For:

1. RESIDUAL vs. MEASURED TCO2:
   ✓ Random scatter around zero = good
   ❌ Systematic pattern (e.g., curve, funnel) = model bias
   ❌ Heteroscedasticity (variance changes) = consider transformation

2. RESIDUAL vs. DEPTH:
   ✓ Random scatter = model captures vertical processes well
   ❌ Surface bias (e.g., positive at 0-100m) = missing biology/gas exchange
   ❌ Deep bias = missing deep remineralization signals

3. RESIDUAL vs. LATITUDE:
   ✓ Random scatter = model captures meridional gradients
   ❌ Tropical bias = missing equatorial processes
   ❌ High-latitude bias = missing polar water mass effects

4. DISTRIBUTION:
   ✓ Normal distribution = error assumptions valid
   ❌ Skewed/heavy tails = outliers or missing nonlinearity
   ❌ Bimodal = two different regimes not captured

5. Q-Q PLOT:
   ✓ Points on diagonal = normality
   ❌ Deviation at tails = extreme values poorly modeled
   ❌ S-curve = systematic skewness

6. SPATIAL MAP:
   ✓ Random spatial distribution = no geographic bias
   ❌ Regional clusters of high/low residuals = missing regional processes
   ❌ Coastal vs. open ocean patterns = missing shelf dynamics

NEXT STEPS based on patterns:
- Depth bias → Add depth-dependent term or stratify model by depth
- Latitude bias → Add latitude interaction terms
- Non-normal residuals → Investigate outliers, consider robust fitting
- Spatial patterns → Add regional water mass classifications
- Heteroscedasticity → Consider weighted regression or transformation
""")

print("="*80)
print("✓ Residual analysis complete. Review plots carefully before publication.")
print("="*80)