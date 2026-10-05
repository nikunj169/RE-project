"""
Spatial and cruise-blocked validation.

Implements GroupKFold-based validation using cruise IDs or spatial blocks.
"""

import numpy as np
import pandas as pd
import logging
from sklearn.model_selection import GroupKFold

from ..analysis.metrics import compute_metrics
from ..basins import assign_basin, get_basin_mask
from ..splits import assign_spatial_blocks, assign_cruise_blocks

logger = logging.getLogger(__name__)


def spatial_blocked_cv(model_class, df, n_splits=5, block_type="spatial",
                       model_kwargs=None):
    """
    Perform GroupKFold cross-validation with spatial or cruise blocks.

    Parameters
    ----------
    model_class : class
        Model class with fit() and predict().
    df : pd.DataFrame
        Full dataset (all basins or single basin).
    n_splits : int
        Number of CV folds.
    block_type : str
        'spatial' or 'cruise'.
    model_kwargs : dict, optional

    Returns
    -------
    dict with per-fold and aggregate metrics.
    """
    model_kwargs = model_kwargs or {}

    S = df["salinity"].values
    T = df["temperature"].values
    A = df["aou"].values
    Y = df["tco2"].values

    if block_type == "cruise" and "cruise" in df.columns:
        groups = df["cruise"].values
    else:
        groups = assign_spatial_blocks(df).values

    gkf = GroupKFold(n_splits=n_splits)
    fold_results = []

    for fold, (train_idx, test_idx) in enumerate(gkf.split(
            np.zeros(len(Y)), Y, groups)):
        try:
            model = model_class(**model_kwargs)
            model.fit(S[train_idx], T[train_idx], A[train_idx], Y[train_idx])
            y_pred = model.predict(S[test_idx], T[test_idx], A[test_idx])
            metrics = compute_metrics(Y[test_idx], y_pred)
            metrics["fold"] = fold
            fold_results.append(metrics)
        except Exception as e:
            logger.warning(f"Fold {fold} failed: {e}")

    if not fold_results:
        return None

    df_folds = pd.DataFrame(fold_results)

    # Aggregate
    agg = {
        "n_folds": len(fold_results),
        "RMSE_mean": df_folds["RMSE"].mean(),
        "RMSE_std": df_folds["RMSE"].std(),
        "MAE_mean": df_folds["MAE"].mean(),
        "MAE_std": df_folds["MAE"].std(),
        "R2_mean": df_folds["R2"].mean(),
        "R2_std": df_folds["R2"].std(),
        "Bias_mean": df_folds["Bias"].mean(),
        "block_type": block_type,
        "n_splits": n_splits,
    }

    return {"per_fold": df_folds, "aggregate": agg}


def compare_temporal_vs_spatial(model_class, df, n_splits=5,
                                 model_kwargs=None):
    """
    Compare temporal validation with spatial/cruise-blocked validation.

    Returns both validation results for comparison.
    """
    from ..splits import split_temporal

    model_kwargs = model_kwargs or {}

    # Temporal validation
    splits = split_temporal(df)
    train_mask = splits["train"]
    val_mask = splits["validation"]

    S = df["salinity"].values
    T = df["temperature"].values
    A = df["aou"].values
    Y = df["tco2"].values

    temporal_result = None
    try:
        model = model_class(**model_kwargs)
        model.fit(S[train_mask], T[train_mask], A[train_mask], Y[train_mask])
        y_pred = model.predict(S[val_mask], T[val_mask], A[val_mask])
        temporal_result = compute_metrics(Y[val_mask], y_pred)
        temporal_result["validation_type"] = "temporal"
    except Exception as e:
        logger.warning(f"Temporal validation failed: {e}")

    # Spatial blocked CV
    spatial_result = spatial_blocked_cv(
        model_class, df, n_splits=n_splits,
        block_type="spatial", model_kwargs=model_kwargs)

    # Cruise blocked CV
    cruise_result = spatial_blocked_cv(
        model_class, df, n_splits=n_splits,
        block_type="cruise", model_kwargs=model_kwargs)

    return {
        "temporal": temporal_result,
        "spatial_blocked": spatial_result,
        "cruise_blocked": cruise_result,
    }
