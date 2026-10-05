"""
Southern Ocean abyssal failure analysis.

Reproduces and investigates the global-vs-local model failure below 3000m
in the Southern Ocean.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
from sklearn.metrics import r2_score, mean_squared_error
import logging

from ..models.rank1_quadratic import rank1_quad_func

logger = logging.getLogger(__name__)


def analyze_abyssal_failure(S, T, A, Y, depth, year,
                             train_end=2015, val_end=2018,
                             depth_min=3000):
    """
    Analyze the Southern Ocean abyssal model failure.

    Compares:
    1. Global basin model applied to abyssal test data
    2. Locally fitted abyssal model
    3. Predicts worse than mean?

    Returns
    -------
    dict with detailed failure analysis.
    """
    # Masks
    abyss = depth >= depth_min
    train = abyss & (year < train_end)
    test = abyss & (year >= train_end) & (year < val_end)

    n_train = train.sum()
    n_test = test.sum()

    if n_train < 50 or n_test < 20:
        return {"error": "Insufficient abyssal data", "n_train": n_train, "n_test": n_test}

    # Fit global model (all basin training data)
    all_train = year < train_end
    p0 = [1.0, -0.3, 2.5, 20.0, 2000.0]

    popt_global, _ = curve_fit(rank1_quad_func,
                                (S[all_train], T[all_train], A[all_train]),
                                Y[all_train], p0=p0, maxfev=30000)

    # Fit local model (abyssal training data only)
    popt_local, _ = curve_fit(rank1_quad_func,
                               (S[train], T[train], A[train]),
                               Y[train], p0=p0, maxfev=30000)

    # Predictions
    X_test = (S[test], T[test], A[test])
    Y_test = Y[test]

    Yp_global = rank1_quad_func(X_test, *popt_global)
    Yp_local = rank1_quad_func(X_test, *popt_local)

    # Mean baseline
    Y_mean = np.full_like(Y_test, np.mean(Y[train]))

    # Metrics
    r2_global = r2_score(Y_test, Yp_global)
    rmse_global = np.sqrt(mean_squared_error(Y_test, Yp_global))
    r2_local = r2_score(Y_test, Yp_local)
    rmse_local = np.sqrt(mean_squared_error(Y_test, Yp_local))
    rmse_mean = np.sqrt(mean_squared_error(Y_test, Y_mean))

    # Residual analysis
    res_global = Y_test - Yp_global
    r_res_T, _ = stats.pearsonr(T[test], res_global)
    r_res_S, _ = stats.pearsonr(S[test], res_global)
    r_res_A, _ = stats.pearsonr(A[test], res_global)

    return {
        "n_train_abyssal": int(n_train),
        "n_test_abyssal": int(n_test),
        "global_r2": r2_global,
        "global_rmse": rmse_global,
        "local_r2": r2_local,
        "local_rmse": rmse_local,
        "mean_rmse": rmse_mean,
        "mean_bias_global": float(np.mean(res_global)),
        "sd_residual_global": float(np.std(res_global)),
        "corr_residual_T": r_res_T,
        "corr_residual_S": r_res_S,
        "corr_residual_A": r_res_A,
        "test_S_range": (float(S[test].min()), float(S[test].max())),
        "test_T_range": (float(T[test].min()), float(T[test].max())),
        "test_AOU_range": (float(A[test].min()), float(A[test].max())),
        "test_TCO2_range": (float(Y_test.min()), float(Y_test.max())),
        "global_coeffs": {k: v for k, v in zip(
            ["alpha", "beta", "gamma", "delta", "epsilon"], popt_global)},
        "local_coeffs": {k: v for k, v in zip(
            ["alpha", "beta", "gamma", "delta", "epsilon"], popt_local)},
        "failure_is_coefficient_transfer": r2_local > 0.5,
    }
