"""
Dependence-aware bootstrap for uncertainty quantification.

Implements:
- Observation-level bootstrap (iid assumption)
- Cluster bootstrap (resample by group)
- Spatial block bootstrap (resample by lat/lon blocks)
"""

import numpy as np
from scipy.optimize import curve_fit
import logging

logger = logging.getLogger(__name__)


def observation_bootstrap(model_class, S, T, A, Y, S_test, T_test, A_test,
                           n_bootstrap=1000, confidence=0.95, random_state=42,
                           model_kwargs=None):
    """
    Standard observation-level bootstrap.

    Resamples training observations with replacement, refits model,
    generates predictions, and adds residual noise.

    Parameters
    ----------
    model_class : class
        Model class with fit() and predict() methods.
    S, T, A, Y : np.ndarray
        Training data.
    S_test, T_test, A_test : np.ndarray
        Test predictors.
    n_bootstrap : int
    confidence : float
    random_state : int
    model_kwargs : dict, optional
        Additional kwargs for model constructor.

    Returns
    -------
    dict with prediction arrays and PI bounds.
    """
    rng = np.random.RandomState(random_state)
    model_kwargs = model_kwargs or {}

    n_train = len(Y)
    n_test = len(S_test)

    # Fit full model for residual sigma
    try:
        full_model = model_class(**model_kwargs)
        full_model.fit(S, T, A, Y)
        train_pred = full_model.predict(S, T, A)
        sigma_resid = np.std(Y - train_pred)
    except Exception as e:
        logger.warning(f"Full model fit failed: {e}")
        sigma_resid = np.std(Y) * 0.1  # fallback

    boot_preds = np.full((n_bootstrap, n_test), np.nan)
    n_success = 0

    for b in range(n_bootstrap):
        idx = rng.choice(n_train, size=n_train, replace=True)
        try:
            model = model_class(**model_kwargs)
            model.fit(S[idx], T[idx], A[idx], Y[idx])
            pred = model.predict(S_test, T_test, A_test)
            # Add residual noise
            noise = rng.normal(0, sigma_resid, n_test)
            boot_preds[b] = pred + noise
            n_success += 1
        except Exception:
            pass

    # Compute PI bounds from successful fits
    valid = boot_preds[~np.any(np.isnan(boot_preds), axis=1)]
    if len(valid) == 0:
        logger.error("All bootstrap fits failed!")
        return None

    alpha = 1 - confidence
    ci_lo = np.percentile(valid, 100 * alpha / 2, axis=0)
    ci_hi = np.percentile(valid, 100 * (1 - alpha / 2), axis=0)

    return {
        "n_bootstrap": n_bootstrap,
        "n_success": n_success,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "mean_pred": np.mean(valid, axis=0),
        "sigma_resid": sigma_resid,
    }


def cluster_bootstrap(model_class, S, T, A, Y, S_test, T_test, A_test,
                       groups, n_bootstrap=1000, confidence=0.95,
                       random_state=42, model_kwargs=None):
    """
    Cluster bootstrap: resample by group (cruise or spatial block).

    Parameters
    ----------
    groups : np.ndarray
        Group labels for each training observation.

    Other parameters same as observation_bootstrap.
    """
    rng = np.random.RandomState(random_state)
    model_kwargs = model_kwargs or {}

    n_test = len(S_test)
    groups = np.asarray(groups)
    unique_groups = np.unique(groups)
    n_groups = len(unique_groups)

    # Fit full model for residual sigma
    try:
        full_model = model_class(**model_kwargs)
        full_model.fit(S, T, A, Y)
        sigma_resid = np.std(Y - full_model.predict(S, T, A))
    except Exception:
        sigma_resid = np.std(Y) * 0.1

    boot_preds = np.full((n_bootstrap, n_test), np.nan)
    n_success = 0

    for b in range(n_bootstrap):
        sampled_groups = rng.choice(unique_groups, size=n_groups, replace=True)
        mask = np.isin(groups, sampled_groups)
        if mask.sum() < 10:
            continue
        try:
            model = model_class(**model_kwargs)
            model.fit(S[mask], T[mask], A[mask], Y[mask])
            pred = model.predict(S_test, T_test, A_test)
            noise = rng.normal(0, sigma_resid, n_test)
            boot_preds[b] = pred + noise
            n_success += 1
        except Exception:
            pass

    valid = boot_preds[~np.any(np.isnan(boot_preds), axis=1)]
    if len(valid) == 0:
        return None

    alpha = 1 - confidence
    ci_lo = np.percentile(valid, 100 * alpha / 2, axis=0)
    ci_hi = np.percentile(valid, 100 * (1 - alpha / 2), axis=0)

    return {
        "n_bootstrap": n_bootstrap,
        "n_success": n_success,
        "n_groups": n_groups,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "mean_pred": np.mean(valid, axis=0),
        "sigma_resid": sigma_resid,
    }
