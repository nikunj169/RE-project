import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error

# 1. Load the clean file
data = np.load('processed_ocean_data.npz')
S, T, AOU, Y = data['S'], data['T'], data['AOU'], data['Y']

# 2. Define your equation
def equation_to_fit(X, alpha, Ac, beta, gamma, delta, epsilon):
    S, T, A = X
    # We use the calibrated Ac as a variable to see if it shifts
    denom = Ac - A
    term1 = alpha * (S / denom)
    return (term1 + beta * T - (gamma * A + delta))**2 + epsilon

# 3. Fit
# p0 = [alpha, Ac, beta, gamma, delta, epsilon]
p0 = [1.0, 3.054312, 1.0, 1.0, 6.475268, 10.15233]

print("Starting calibration...")
try:
    popt, _ = curve_fit(equation_to_fit, (S, T, AOU), Y, p0=p0, maxfev=20000)
    
    # 4. Accuracy Calculation
    Y_pred = equation_to_fit((S, T, AOU), *popt)
    r2 = r2_score(Y, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y, Y_pred))
    
    print("\n--- FINAL RESULTS ---")
    print(f"Optimal Ac (Singularity): {popt[1]:.6f}")
    print(f"Optimal epsilon (Offset): {popt[5]:.6f}")
    print(f"Final Accuracy (R^2): {r2:.4f}")
    print(f"Final Error (RMSE): {rmse:.4f}")
    
except np.linalg.LinAlgError:
    print("Error: The data is too noisy or initial guess is too far off.")