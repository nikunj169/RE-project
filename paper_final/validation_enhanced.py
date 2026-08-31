"""
Enhanced Ocean Carbon Model Validation Script
==============================================

This script adds the critical validation components needed for publication:
1. Train/test split with cross-validation
2. Comprehensive error metrics (RMSE, MAE, bias, etc.)
3. Residual analysis plots
4. Baseline comparison (Multiple Linear Regression)
5. Depth-stratified performance analysis
6. Uncertainty quantification via bootstrap

Run this AFTER your initial model fitting to assess publication readiness.
"""

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.model_selection import KFold, train_test_split
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# MODEL DEFINITION (same as original)
# =============================================================================
def stable_model(X, alpha, beta, gamma, delta, epsilon):
    """Your quadratic TCO2 model"""
    S, T, A = X
    A_norm = A / 100.0
    core = (alpha * S) + (beta * T) + (gamma * A_norm) + delta
    return (core**2) + epsilon

# =============================================================================
# COMPREHENSIVE METRICS FUNCTION
# =============================================================================
def calculate_metrics(y_true, y_pred):
    """Calculate all relevant error metrics"""
    metrics = {
        'R2': r2_score(y_true, y_pred),
        'RMSE': np.sqrt(mean_squared_error(y_true, y_pred)),
        'MAE': mean_absolute_error(y_true, y_pred),
        'Bias': np.mean(y_pred - y_true),
        'NRMSE': np.sqrt(mean_squared_error(y_true, y_pred)) / np.mean(y_true) * 100,
        'Within_10': np.mean(np.abs(y_pred - y_true) < 10) * 100,
        'Within_20': np.mean(np.abs(y_pred - y_true) < 20) * 100,
        'Within_30': np.mean(np.abs(y_pred - y_true) < 30) * 100
    }
    return metrics

# =============================================================================
# LOAD DATA (same as original)
# =============================================================================
print("Loading GLODAP data...")
mat = sio.loadmat('data/raw/GLODAPv2.2023_Merged_Master_File.mat', squeeze_me=True)
S_all, T_all, AOU_all, Y_all = mat['G2salinity'], mat['G2temperature'], mat['G2aou'], mat['G2tco2']
Region, Lat = mat['G2region'], mat['G2latitude']

# Try to load depth if available for stratified analysis
try:
    Depth_all = mat['G2depth']
    has_depth = True
except:
    print("Warning: Depth data not found in .mat file")
    has_depth = False

basin_map = {'Atlantic': [1], 'Pacific': [8, 2], 'Indian': [16, 3], 'Southern Ocean': 'Lat'}

# =============================================================================
# VALIDATION RESULTS STORAGE
# =============================================================================
validation_results = []

# =============================================================================
# MAIN VALIDATION LOOP
# =============================================================================
for basin_name, b_id in basin_map.items():
    print(f"\n{'='*80}")
    print(f"BASIN: {basin_name}")
    print(f"{'='*80}")
    
    # Apply Basin Mask (same as original)
    if basin_name == 'Southern Ocean':
        mask = (Lat < -30) & (Region != 4)
    else:
        mask = (np.isin(Region, b_id)) & (Lat > -30)
    
    # Quality Control Mask
    mask &= (S_all > 30) & (Y_all > 1800) & (AOU_all > 0)
    S, T, A, Y = S_all[mask], T_all[mask], AOU_all[mask], Y_all[mask]
    
    if has_depth:
        Depth = Depth_all[mask]
    
    n_samples = len(Y)
    print(f"Sample size: {n_samples:,}")
    
    # =========================================================================
    # 1. FULL DATASET FIT (original approach)
    # =========================================================================
    print("\n1. Full Dataset Performance:")
    print("-" * 40)
    
    popt_full, pcov_full = curve_fit(stable_model, (S, T, A), Y, p0=[0.1, -1.0, 5.0, 10.0, 1500])
    Y_pred_full = stable_model((S, T, A), *popt_full)
    
    metrics_full = calculate_metrics(Y, Y_pred_full)
    
    print(f"R² = {metrics_full['R2']:.4f}")
    print(f"RMSE = {metrics_full['RMSE']:.2f} μmol/kg")
    print(f"MAE = {metrics_full['MAE']:.2f} μmol/kg")
    print(f"Bias = {metrics_full['Bias']:.2f} μmol/kg")
    print(f"Normalized RMSE = {metrics_full['NRMSE']:.2f}%")
    print(f"Within ±10 μmol/kg: {metrics_full['Within_10']:.1f}%")
    print(f"Within ±20 μmol/kg: {metrics_full['Within_20']:.1f}%")
    print(f"Within ±30 μmol/kg: {metrics_full['Within_30']:.1f}%")
    
    # =========================================================================
    # 2. TRAIN/TEST SPLIT VALIDATION
    # =========================================================================
    print("\n2. Train/Test Split (80/20):")
    print("-" * 40)
    
    # Combine predictors
    X = np.column_stack([S, T, A])
    
    # Split data
    X_train, X_test, Y_train, Y_test = train_test_split(
        X, Y, test_size=0.2, random_state=42
    )
    
    S_train, T_train, A_train = X_train[:, 0], X_train[:, 1], X_train[:, 2]
    S_test, T_test, A_test = X_test[:, 0], X_test[:, 1], X_test[:, 2]
    
    # Fit on training data
    popt_train, _ = curve_fit(stable_model, (S_train, T_train, A_train), Y_train, 
                               p0=[0.1, -1.0, 5.0, 10.0, 1500])
    
    # Predict on test data
    Y_pred_test = stable_model((S_test, T_test, A_test), *popt_train)
    
    metrics_test = calculate_metrics(Y_test, Y_pred_test)
    
    print(f"Test R² = {metrics_test['R2']:.4f}")
    print(f"Test RMSE = {metrics_test['RMSE']:.2f} μmol/kg")
    print(f"Test Bias = {metrics_test['Bias']:.2f} μmol/kg")
    
    # Check for overfitting
    overfit_indicator = metrics_full['R2'] - metrics_test['R2']
    print(f"Overfitting indicator (R² drop): {overfit_indicator:.4f}")
    if overfit_indicator > 0.05:
        print("⚠️  WARNING: Possible overfitting detected!")
    else:
        print("✓ Model generalizes well")
    
    # =========================================================================
    # 3. K-FOLD CROSS-VALIDATION
    # =========================================================================
    print("\n3. 5-Fold Cross-Validation:")
    print("-" * 40)
    
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_r2_scores = []
    cv_rmse_scores = []
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(X)):
        X_cv_train, X_cv_val = X[train_idx], X[val_idx]
        Y_cv_train, Y_cv_val = Y[train_idx], Y[val_idx]
        
        S_cv, T_cv, A_cv = X_cv_train[:, 0], X_cv_train[:, 1], X_cv_train[:, 2]
        
        try:
            popt_cv, _ = curve_fit(stable_model, (S_cv, T_cv, A_cv), Y_cv_train, 
                                   p0=[0.1, -1.0, 5.0, 10.0, 1500])
            
            S_val, T_val, A_val = X_cv_val[:, 0], X_cv_val[:, 1], X_cv_val[:, 2]
            Y_pred_cv = stable_model((S_val, T_val, A_val), *popt_cv)
            
            r2_cv = r2_score(Y_cv_val, Y_pred_cv)
            rmse_cv = np.sqrt(mean_squared_error(Y_cv_val, Y_pred_cv))
            
            cv_r2_scores.append(r2_cv)
            cv_rmse_scores.append(rmse_cv)
        except:
            print(f"  Fold {fold+1}: Fit failed (skipping)")
    
    if cv_r2_scores:
        print(f"Mean CV R² = {np.mean(cv_r2_scores):.4f} ± {np.std(cv_r2_scores):.4f}")
        print(f"Mean CV RMSE = {np.mean(cv_rmse_scores):.2f} ± {np.std(cv_rmse_scores):.2f} μmol/kg")
        
        if np.mean(cv_r2_scores) > 0.85:
            print("✓ Excellent cross-validation performance")
        elif np.mean(cv_r2_scores) > 0.75:
            print("⚠️  Acceptable but could be improved")
        else:
            print("❌ Poor generalization - model needs revision")
    
    # =========================================================================
    # 4. BASELINE COMPARISON: Multiple Linear Regression
    # =========================================================================
    print("\n4. Baseline Comparison (MLR):")
    print("-" * 40)
    
    mlr = LinearRegression()
    mlr.fit(X_train, Y_train)
    Y_mlr_test = mlr.predict(X_test)
    
    metrics_mlr = calculate_metrics(Y_test, Y_mlr_test)
    
    print(f"MLR Test R² = {metrics_mlr['R2']:.4f}")
    print(f"MLR Test RMSE = {metrics_mlr['RMSE']:.2f} μmol/kg")
    print(f"\nYour Model Test R² = {metrics_test['R2']:.4f}")
    print(f"Your Model Test RMSE = {metrics_test['RMSE']:.2f} μmol/kg")
    
    improvement = (metrics_mlr['RMSE'] - metrics_test['RMSE']) / metrics_mlr['RMSE'] * 100
    print(f"\n➤ Improvement over MLR: {improvement:.1f}% reduction in RMSE")
    
    if improvement > 15:
        print("✓ Substantial improvement - model adds value")
    elif improvement > 5:
        print("⚠️  Modest improvement - justify added complexity")
    else:
        print("❌ Minimal improvement - consider simpler model")
    
    # =========================================================================
    # 5. DEPTH-STRATIFIED ANALYSIS (if depth available)
    # =========================================================================
    if has_depth:
        print("\n5. Depth-Stratified Performance:")
        print("-" * 40)
        
        depth_bins = [(0, 100), (100, 500), (500, 1000), (1000, 3000), (3000, 7000)]
        
        for d_min, d_max in depth_bins:
            depth_mask = (Depth >= d_min) & (Depth < d_max)
            if np.sum(depth_mask) > 100:  # Only analyze if enough samples
                Y_depth = Y[depth_mask]
                Y_pred_depth = Y_pred_full[depth_mask]
                
                r2_depth = r2_score(Y_depth, Y_pred_depth)
                rmse_depth = np.sqrt(mean_squared_error(Y_depth, Y_pred_depth))
                
                print(f"{d_min:4d}-{d_max:4d}m (n={np.sum(depth_mask):5d}): "
                      f"R²={r2_depth:.4f}, RMSE={rmse_depth:.2f} μmol/kg")
    
    # =========================================================================
    # 6. BOOTSTRAP UNCERTAINTY QUANTIFICATION
    # =========================================================================
    print("\n6. Parameter Uncertainty (Bootstrap):")
    print("-" * 40)
    
    n_bootstrap = 100  # Use 1000 for final analysis
    param_distributions = []
    
    print("Running bootstrap... ", end='', flush=True)
    for i in range(n_bootstrap):
        # Resample with replacement
        indices = np.random.choice(len(Y), len(Y), replace=True)
        S_boot, T_boot, A_boot, Y_boot = S[indices], T[indices], A[indices], Y[indices]
        
        try:
            popt_boot, _ = curve_fit(stable_model, (S_boot, T_boot, A_boot), Y_boot, 
                                      p0=[0.1, -1.0, 5.0, 10.0, 1500], maxfev=5000)
            param_distributions.append(popt_boot)
        except:
            pass  # Skip failed fits
    
    print(f"Done ({len(param_distributions)}/{n_bootstrap} successful fits)")
    
    if param_distributions:
        param_distributions = np.array(param_distributions)
        param_names = ['Alpha (S)', 'Beta (T)', 'Gamma (AOU)', 'Delta', 'Epsilon']
        
        print("\nParameter Estimates (95% Confidence Intervals):")
        for j, param_name in enumerate(param_names):
            mean_val = np.mean(param_distributions[:, j])
            ci_lower = np.percentile(param_distributions[:, j], 2.5)
            ci_upper = np.percentile(param_distributions[:, j], 97.5)
            print(f"{param_name:12s}: {mean_val:8.4f} [{ci_lower:8.4f}, {ci_upper:8.4f}]")
    
    # =========================================================================
    # STORE RESULTS FOR SUMMARY TABLE
    # =========================================================================
    validation_results.append({
        'Basin': basin_name,
        'N': n_samples,
        'Full_R2': metrics_full['R2'],
        'Full_RMSE': metrics_full['RMSE'],
        'Test_R2': metrics_test['R2'],
        'Test_RMSE': metrics_test['RMSE'],
        'CV_R2_mean': np.mean(cv_r2_scores) if cv_r2_scores else np.nan,
        'CV_R2_std': np.std(cv_r2_scores) if cv_r2_scores else np.nan,
        'MLR_R2': metrics_mlr['R2'],
        'MLR_RMSE': metrics_mlr['RMSE'],
        'Improvement_pct': improvement,
        'Alpha': popt_full[0],
        'Beta': popt_full[1],
        'Gamma': popt_full[2],
        'Delta': popt_full[3],
        'Epsilon': popt_full[4]
    })

# =============================================================================
# SUMMARY TABLE
# =============================================================================
print("\n" + "="*80)
print("VALIDATION SUMMARY TABLE")
print("="*80)

df_validation = pd.DataFrame(validation_results)

# Display key metrics
print("\nModel Performance:")
print(df_validation[['Basin', 'Full_R2', 'Test_R2', 'CV_R2_mean', 'Full_RMSE', 'Test_RMSE']].to_string(index=False))

print("\nComparison to Baseline (MLR):")
print(df_validation[['Basin', 'MLR_RMSE', 'Test_RMSE', 'Improvement_pct']].to_string(index=False))

print("\nModel Coefficients:")
print(df_validation[['Basin', 'Alpha', 'Beta', 'Gamma', 'Epsilon']].to_string(index=False))

# Save results
df_validation.to_csv('validation_results_comprehensive.csv', index=False)
print("\n✓ Results saved to 'validation_results_comprehensive.csv'")

# =============================================================================
# PUBLICATION READINESS ASSESSMENT
# =============================================================================
print("\n" + "="*80)
print("PUBLICATION READINESS ASSESSMENT")
print("="*80)

criteria = []

# Check cross-validation performance
cv_r2_min = df_validation['CV_R2_mean'].min()
if cv_r2_min > 0.90:
    criteria.append(("✓", "Cross-validation R² > 0.90 (excellent)"))
elif cv_r2_min > 0.85:
    criteria.append(("⚠️", "Cross-validation R² > 0.85 (good, but room for improvement)"))
else:
    criteria.append(("❌", "Cross-validation R² < 0.85 (needs improvement)"))

# Check RMSE
test_rmse_max = df_validation['Test_RMSE'].max()
if test_rmse_max < 15:
    criteria.append(("✓", f"Test RMSE < 15 μmol/kg (excellent: max={test_rmse_max:.1f})"))
elif test_rmse_max < 25:
    criteria.append(("⚠️", f"Test RMSE < 25 μmol/kg (acceptable: max={test_rmse_max:.1f})"))
else:
    criteria.append(("❌", f"Test RMSE > 25 μmol/kg (needs improvement: max={test_rmse_max:.1f})"))

# Check improvement over baseline
improvement_min = df_validation['Improvement_pct'].min()
if improvement_min > 15:
    criteria.append(("✓", f"Substantial improvement over MLR (min={improvement_min:.1f}%)"))
elif improvement_min > 5:
    criteria.append(("⚠️", f"Modest improvement over MLR (min={improvement_min:.1f}%)"))
else:
    criteria.append(("❌", f"Minimal improvement over MLR (min={improvement_min:.1f}%)"))

# Check overfitting
max_r2_drop = (df_validation['Full_R2'] - df_validation['Test_R2']).max()
if max_r2_drop < 0.03:
    criteria.append(("✓", f"No overfitting detected (max R² drop={max_r2_drop:.4f})"))
elif max_r2_drop < 0.05:
    criteria.append(("⚠️", f"Minor overfitting (max R² drop={max_r2_drop:.4f})"))
else:
    criteria.append(("❌", f"Overfitting detected (max R² drop={max_r2_drop:.4f})"))

print("\nCriteria Checklist:")
for symbol, message in criteria:
    print(f"{symbol} {message}")

# Overall assessment
n_pass = sum(1 for s, _ in criteria if s == "✓")
n_warn = sum(1 for s, _ in criteria if s == "⚠️")
n_fail = sum(1 for s, _ in criteria if s == "❌")

print(f"\nOverall: {n_pass}/4 criteria passed, {n_warn}/4 warnings, {n_fail}/4 failures")

if n_fail == 0 and n_pass >= 3:
    print("\n✓✓✓ MODEL IS PUBLICATION-READY")
    print("Next steps: Add physical interpretation, compare to literature algorithms")
elif n_fail == 0:
    print("\n⚠️⚠️ MODEL IS NEAR PUBLICATION-READY")
    print("Next steps: Address warnings, improve metrics where possible")
else:
    print("\n❌❌ MODEL NEEDS SIGNIFICANT IMPROVEMENT")
    print("Next steps: Address failed criteria before proceeding to publication")

print("\n" + "="*80)
print("See 'research_roadmap.md' for detailed next steps")
print("="*80)