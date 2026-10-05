"""
Depth-stratified performance analysis.
"""

import numpy as np
import pandas as pd
from .metrics import compute_metrics


DEPTH_BINS = [0, 100, 500, 1000, 3000, 7000]
DEPTH_LABELS = ["0-100m", "100-500m", "500-1000m", "1000-3000m", "3000-7000m"]


def depth_stratified_metrics(y_true, y_pred, depth, bins=None, labels=None):
    """
    Compute metrics per depth bin.

    Parameters
    ----------
    y_true, y_pred, depth : np.ndarray
    bins : list, optional
    labels : list, optional

    Returns
    -------
    pd.DataFrame
        One row per depth bin with metrics and sample size.
    """
    if bins is None:
        bins = DEPTH_BINS
    if labels is None:
        labels = DEPTH_LABELS

    results = []
    for i in range(len(bins) - 1):
        mask = (depth >= bins[i]) & (depth < bins[i + 1])
        if mask.sum() < 10:
            continue
        m = compute_metrics(y_true[mask], y_pred[mask])
        m["depth_bin"] = labels[i]
        m["depth_min"] = bins[i]
        m["depth_max"] = bins[i + 1]
        results.append(m)

    return pd.DataFrame(results)


def depth_stratified_comparison(y_true, y_pred_a, y_pred_b, depth,
                                 label_a="ModelA", label_b="ModelB",
                                 bins=None, labels=None):
    """
    Compare two models in each depth bin.
    """
    if bins is None:
        bins = DEPTH_BINS
    if labels is None:
        labels = DEPTH_LABELS

    results = []
    for i in range(len(bins) - 1):
        mask = (depth >= bins[i]) & (depth < bins[i + 1])
        if mask.sum() < 10:
            continue

        m_a = compute_metrics(y_true[mask], y_pred_a[mask])
        m_b = compute_metrics(y_true[mask], y_pred_b[mask])

        improvement = (m_b["RMSE"] - m_a["RMSE"]) / m_b["RMSE"] * 100

        results.append({
            "depth_bin": labels[i],
            "n": int(mask.sum()),
            f"RMSE_{label_a}": m_a["RMSE"],
            f"RMSE_{label_b}": m_b["RMSE"],
            f"R2_{label_a}": m_a["R2"],
            f"R2_{label_b}": m_b["R2"],
            "Improvement_pct": improvement,
            f"Bias_{label_a}": m_a["Bias"],
            f"Bias_{label_b}": m_b["Bias"],
        })

    return pd.DataFrame(results)
