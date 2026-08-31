import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# 1. LOAD DATA
mat = sio.loadmat('data/raw/GLODAPv2.2023_Merged_Master_File.mat', squeeze_me=True)
S_all, T_all, AOU_all, Y_all = mat['G2salinity'], mat['G2temperature'], mat['G2aou'], mat['G2tco2']
Region, Lat = mat['G2region'], mat['G2latitude']

# 2. THE STABLE QUADRATIC MODEL
def stable_model(X, alpha, beta, gamma, delta, epsilon):
    S, T, A = X
    A_norm = A / 100.0
    core = (alpha * S) + (beta * T) + (gamma * A_norm) + delta
    return (core**2) + epsilon

# 3. INITIALIZATION
basin_map = {'Atlantic': [1], 'Pacific': [8, 2], 'Indian': [16, 3], 'Southern Ocean': 'Lat'}
results_table = []

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
axes = axes.flatten()

print(f"{'Basin':<15} | {'R2':<8} | {'Alpha':<8} | {'Beta':<8} | {'Gamma':<8} | {'Epsilon':<8}")
print("-" * 75)

# 4. PROCESSING LOOP
for i, (name, b_id) in enumerate(basin_map.items()):
    # Apply Basin Mask
    if name == 'Southern Ocean':
        mask = (Lat < -30) & (Region != 4)
    else:
        mask = (np.isin(Region, b_id)) & (Lat > -30)
    
    # Quality Control Mask
    mask &= (S_all > 30) & (Y_all > 1800) & (AOU_all > 0)
    S, T, A, Y = S_all[mask], T_all[mask], AOU_all[mask], Y_all[mask]
    
    # PERFORM FIT
    # p0: alpha, beta, gamma, delta, epsilon
    popt, pcov = curve_fit(stable_model, (S, T, A), Y, p0=[0.1, -1.0, 5.0, 10.0, 1500])
    
    # Calculate Accuracy
    Y_pred = stable_model((S, T, A), *popt)
    r2 = r2_score(Y, Y_pred)
    
    # Store Results
    results_table.append({
        'Basin': name,
        'R2': round(r2, 4),
        'Alpha (S)': round(popt[0], 4),
        'Beta (T)': round(popt[1], 4),
        'Gamma (AOU)': round(popt[2], 4),
        'Delta': round(popt[3], 4),
        'Epsilon': round(popt[4], 2)
    })

    # 5. VISUALIZATION
    ax = axes[i]
    hb = ax.hexbin(Y, Y_pred, gridsize=60, cmap='YlGnBu', bins='log', mincnt=1)
    ax.plot([1900, 2500], [1900, 2500], 'r--', lw=2.5, label='1:1 Line')
    
    ax.set_title(f"{name}\n$R^2 = {r2:.4f}$", fontsize=15, fontweight='bold')
    ax.set_xlabel("Measured $TCO_2$ ($\mu mol/kg$)", fontsize=12)
    ax.set_ylabel("Predicted $TCO_2$ ($\mu mol/kg$)", fontsize=12)
    ax.set_xlim([1950, 2400])
    ax.set_ylim([1950, 2400])
    ax.grid(True, linestyle=':', alpha=0.6)

    print(f"{name:<15} | {r2:<8.4f} | {popt[0]:<8.4f} | {popt[1]:<8.4f} | {popt[2]:<8.4f} | {popt[4]:<8.2f}")

# 6. EXPORT & FINALIZE
plt.tight_layout()
plt.savefig('Final_Carbon_Validation_Plots.png', dpi=300)
plt.show()

# Save coefficients to CSV
df = pd.DataFrame(results_table)
df.to_csv('ocean_carbon_coefficients.csv', index=False)
print("\nSuccess! Coefficients saved to 'ocean_carbon_coefficients.csv'")