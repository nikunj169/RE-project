import numpy as np
import scipy.io as sio
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import pandas as pd

# 1. LOAD DATA
print("Loading GLODAP Master File...")
mat = sio.loadmat('data/raw/GLODAPv2.2023_Merged_Master_File.mat', squeeze_me=True)

# Extract core variables
S_all = mat['G2salinity']
T_all = mat['G2temperature']
AOU_all = mat['G2aou']
Y_all = mat['G2tco2'] # Target: Total Carbon (or change to your target)
Region = mat['G2region'] # The column with the IDs
Lat = mat['G2latitude']  # Needed for Southern Ocean separation

# 2. DEFINE THE UNIVERSAL EQUATION
def ocean_model(X, alpha, Ac, beta, gamma, delta, epsilon):
    S, T, A = X
    # Ac is the singularity. We add a tiny buffer (1e-6) to avoid division by zero
    denom = Ac - A 
    # Structural Form: ( alpha*(S/(Ac-A)) + beta*T - (gamma*A + delta) )^2 + epsilon
    term1 = alpha * (S / (denom + 1e-6))
    return (term1 + beta * T - (gamma * A + delta))**2 + epsilon

# 3. DEFINE BASINS (Using Your Bitwise IDs)
# We use a dictionary to map IDs to Names. 
# We also allow for 'Standard' IDs (2,3) just in case they appear mixed.
basin_defs = {
    'Atlantic': [1],        # Bitwise 1
    'Pacific':  [8, 2],     # Bitwise 8, or Standard 2
    'Indian':   [16, 3],    # Bitwise 16, or Standard 3
    'Arctic':   [4],        # Bitwise 4
}

# 4. CONSTRAINTS (Physical Realism)
# [alpha, Ac, beta, gamma, delta, epsilon]
# We force Alpha/Beta > 0.01 to ensure S and T actually contribute (no "parameter vanishing")
lower_bounds = [0.01, 2.50, 0.01, 0.01, -200, -10000]
upper_bounds = [50.0, 4.00, 50.0, 50.0,  500,  10000]
p0 = [1.0, 3.05, 1.0, 1.0, 6.47, 10.15]

results_list = []

print("\n" + "="*80)
print(f"{'BASIN':<12} | {'Ac (Singular)':<14} | {'Alpha (S)':<10} | {'Beta (T)':<10} | {'R2 Score':<8} | {'Status'}")
print("="*80)

# 5. EXECUTION LOOP
# We loop through our defined basins, plus a custom "Southern Ocean" check
basins_to_process = list(basin_defs.keys()) + ['Southern Ocean']

for name in basins_to_process:
    
    # --- STEP A: FILTERING ---
    if name == 'Southern Ocean':
        # Southern Ocean is defined geographically (South of 30S or 35S)
        # We exclude the Arctic (ID 4) just to be safe
        mask = (Lat < -30) & (Region != 4) & (S_all > 0) & (Y_all > -900) & (np.abs(AOU_all - 3.0) > 0.1)
    else:
        # For other oceans, we filter by their specific IDs AND exclude the Southern Ocean (Lat > -30)
        # This prevents "double counting" the southern parts of Atlantic/Pacific
        target_ids = basin_defs[name]
        mask = (np.isin(Region, target_ids)) & (Lat > -30) & (S_all > 0) & (Y_all > -900) & (np.abs(AOU_all - 3.0) > 0.1)

    # --- STEP B: CHECK DATA COUNT ---
    count = np.sum(mask)
    if count < 1000:
        print(f"{name:<12} | {'-':<14} | {'-':<10} | {'-':<10} | {'-':<8} | Skipped (<1k pts)")
        continue

    # --- STEP C: FITTING ---
    S, T, AOU, Y = S_all[mask], T_all[mask], AOU_all[mask], Y_all[mask]
    
    try:
        # Use 'trf' method for robust bounds handling
        popt, _ = curve_fit(ocean_model, (S, T, AOU), Y, p0=p0, 
                            bounds=(lower_bounds, upper_bounds), 
                            method='trf', maxfev=50000)
        
        # Calculate Accuracy
        Y_pred = ocean_model((S, T, AOU), *popt)
        r2 = r2_score(Y, Y_pred)
        
        # Print Row
        print(f"{name:<12} | {popt[1]:.5f}        | {popt[0]:.3f}      | {popt[2]:.3f}      | {r2:.4f}   | Success")
        
        # Save for summary
        results_list.append({
            'Basin': name, 'Ac': popt[1], 'Alpha': popt[0], 'Beta': popt[2], 
            'Gamma': popt[3], 'Delta': popt[4], 'Epsilon': popt[5], 'R2': r2
        })
        
    except Exception as e:
        print(f"{name:<12} | {'ERROR':<14} | {'-':<10} | {'-':<10} | {'-':<8} | Failed: {str(e)[:20]}")

# 6. SUMMARY STATISTICS
if results_list:
    df = pd.DataFrame(results_list)
    print("\n" + "="*80)
    print("FINAL CONSTANT ANALYSIS")
    print("="*80)
    print(df[['Basin', 'Ac', 'Alpha', 'Beta', 'Epsilon', 'R2']])
    print("\nGlobal Ac Stability (Std Dev):", df['Ac'].std())