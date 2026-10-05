"""
Metric computation functions.

All metrics are computed from raw arrays — no model assumptions.
"""

import numpy as np
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error


def compute_metrics(y_true, y_pred):
    """
    Compute a comprehensive set of error metrics.

    Parameters
    ----------
    y_true, y_pred : array-like

    Returns
    -------
    dict
        Metric name -> value.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    residuals = y_true - y_pred
    abs_errors = np.abs(residuals)

    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    bias = np.mean(residuals)
    mean_true = np.mean(y_true)

    return {
        "n": len(y_true),
        "RMSE": rmse,
        "MAE": mae,
        "R2": r2,
        "Bias": bias,
        "NRMSE_pct": (rmse / mean_true * 100) if mean_true != 0 else np.nan,
        "Mean_AE": np.mean(abs_errors),
        "Median_AE": np.median(abs_errors),
        "P90_AE": np.percentile(abs_errors, 90),
        "Within_5": np.mean(abs_errors < 5) * 100,
        "Within_10": np.mean(abs_errors < 10) * 100,
        "Within_20": np.mean(abs_errors < 20) * 100,
        "Within_30": np.mean(abs_errors < 30) * 100,
    }


def compute_paired_metrics(y_true, y_pred_a, y_pred_b, label_a="A", label_b="B"):
    """
    Compute paired comparison metrics between two models.

    Positive values indicate model A is better (lower error).

    Returns
    -------
    dict with paired error statistics.
    """
    ae_a = np.abs(y_true - y_pred_a)
    ae_b = np.abs(y_true - y_pred_b)
    d = ae_b - ae_a  # positive = A is better

    return {
        "n": len(d),
        f"RMSE_{label_a}": np.sqrt(mean_squared_error(y_true, y_pred_a)),
        f"RMSE_{label_b}": np.sqrt(mean_squared_error(y_true, y_pred_b)),
        "paired_AE_mean_diff": np.mean(d),
        "paired_AE_median_diff": np.median(d),
        "paired_AE_std_diff": np.std(d, ddof=1),
        f"wins_{label_a}": int(np.sum(d > 0)),
        f"wins_{label_b}": int(np.sum(d < 0)),
        "ties": int(np.sum(d == 0)),
    }
