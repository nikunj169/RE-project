import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

# 1. LOAD DATA
mat = sio.loadmat('data/raw/GLODAPv2.2023_Merged_Master_File.mat', squeeze_me=True)
S_all, T_all, AOU_all, Y_all = mat['G2salinity'], mat['G2temperature'], mat['G2aou'], mat['G2tco2']
Region = mat['G2region']

def stable_model(X, alpha, beta, gamma, delta, epsilon):
    S, T, A = X
    A_norm = A / 100.0
    core = (alpha * S) + (beta * T) + (gamma * A_norm) + delta
    return (core**2) + epsilon

# 2. YOUR SUCCESSFUL PACIFIC CONSTANTS
# Using the values you just generated
pac_params = [1.1962, -0.2874, 1.821, -22.2795, 1828.43]

# 3. APPLY TO ARCTIC DATA (Region 4)
mask = (Region == 4) & (S_all > 15) & (Y_all > 1500) & (AOU_all > 0)
S_arc, T_arc, A_arc, Y_arc = S_all[mask], T_all[mask], AOU_all[mask], Y_all[mask]

# Predict Arctic TCO2 using Pacific Physics
Y_pred = stable_model((S_arc, T_arc, A_arc), *pac_params)
r2_cross = r2_score(Y_arc, Y_pred)

# 4. PLOT THE "FAILURE"
plt.figure(figsize=(8, 6))
plt.hexbin(Y_arc, Y_pred, gridsize=50, cmap='plasma', bins='log', mincnt=1)
plt.plot([1900, 2450], [1900, 2450], 'r--', lw=2, label='1:1 Line (Perfect Fit)')
plt.title(f"Cross-Validation: Pacific Model on Arctic Data\n$R^2 = {r2_cross:.4f}$", fontsize=12)
plt.xlabel("Actual Arctic $TCO_2$")
plt.ylabel("Predicted $TCO_2$ (Using Pacific Constants)")
plt.show()

print(f"Cross-Basin R2: {r2_cross}")