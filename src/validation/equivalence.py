"""
Equivalence testing: TOST and dependence-aware alternatives.

Implements Two One-Sided Tests (TOST) on paired absolute-error differences,
with support for multiple equivalence margins and dependence-aware bootstrap.
"""

import numpy as np
from scipy import stats
import logging

logger = logging.getLogger(__name__)


def tost_paired_ae(y_true, y_pred_a, y_pred_b, delta, alpha=0.05):
    """
    TOST on per-sample paired absolute-error differences.

    d_i = |y_true_i - y_pred_b_i| - |y_true_i - y_pred_a_i|
    Positive d_i means model A is more accurate.

    Parameters
    ----------
    y_true, y_pred_a, y_pred_b : np.ndarray
    delta : float
        Equivalence margin (μmol/kg).
    alpha : float
        Significance level.

    Returns
    -------
    dict with TOST results.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred_a = np.asarray(y_pred_a, dtype=float)
    y_pred_b = np.asarray(y_pred_b, dtype=float)

    ae_a = np.abs(y_true - y_pred_a)
    ae_b = np.abs(y_true - y_pred_b)
    d = ae_b - ae_a  # positive = A better

    n = len(d)
    d_bar = np.mean(d)
    d_sd = np.std(d, ddof=1)
    se = d_sd / np.sqrt(n)

    # TOST: two one-sided tests
    # Test 1: H0: d_bar <= -delta vs H1: d_bar > -delta
    # Reject if d_bar is significantly above -delta
    t1 = (d_bar - (-delta)) / se
    p1 = 1.0 - stats.t.cdf(t1, df=n - 1)  # right-tail p-value

    # Test 2: H0: d_bar >= +delta vs H1: d_bar < +delta
    # Reject if d_bar is significantly below +delta
    t2 = (d_bar - delta) / se
    p2 = stats.t.cdf(t2, df=n - 1)  # left-tail p-value

    # TOST p-value: max of individual one-sided p-values
    # Small tost_p means both tests reject → equivalence
    tost_p = max(p1, p2)
    equiv = tost_p < alpha

    # Confidence interval
    tcrit = stats.t.ppf(1 - alpha, df=n - 1)
    ci_lo = d_bar - tcrit * se
    ci_hi = d_bar + tcrit * se
    ci_within = (ci_lo > -delta) and (ci_hi < delta)

    return {
        "n": n,
        "delta": delta,
        "alpha": alpha,
        "mean_AE_diff": d_bar,
        "SE": se,
        "tost_p": tost_p,
        "equiv_by_test": equiv,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "ci_within_bounds": ci_within,
    }


def tost_multiple_margins(y_true, y_pred_a, y_pred_b,
                          margins=(1.0, 2.0, 5.0, 10.0), alpha=0.05):
    """
    Run TOST at multiple equivalence margins.

    Returns list of dicts, one per margin.
    """
    results = []
    for delta in margins:
        result = tost_paired_ae(y_true, y_pred_a, y_pred_b, delta, alpha)
        results.append(result)
    return results


def bootstrap_paired_ae_diff(y_true, y_pred_a, y_pred_b,
                              n_bootstrap=1000, confidence=0.95,
                              groups=None, random_state=42):
    """
    Bootstrap confidence interval for paired AE difference.

    If groups are provided, uses cluster bootstrap (resample by group).

    Parameters
    ----------
    y_true, y_pred_a, y_pred_b : np.ndarray
    n_bootstrap : int
    confidence : float
    groups : np.ndarray, optional
        Group labels for cluster bootstrap.
    random_state : int

    Returns
    -------
    dict with bootstrap CI for mean paired AE difference.
    """
    rng = np.random.RandomState(random_state)
    y_true = np.asarray(y_true, dtype=float)
    y_pred_a = np.asarray(y_pred_a, dtype=float)
    y_pred_b = np.asarray(y_pred_b, dtype=float)

    ae_a = np.abs(y_true - y_pred_a)
    ae_b = np.abs(y_true - y_pred_b)
    d = ae_b - ae_a  # positive = A better

    boot_means = []
    alpha = 1 - confidence

    if groups is not None:
        groups = np.asarray(groups)
        unique_groups = np.unique(groups)
        n_groups = len(unique_groups)

        for _ in range(n_bootstrap):
            # Resample groups with replacement
            sampled_groups = rng.choice(unique_groups, size=n_groups, replace=True)
            mask = np.isin(groups, sampled_groups)
            if mask.sum() > 0:
                boot_means.append(np.mean(d[mask]))
    else:
        n = len(d)
        for _ in range(n_bootstrap):
            idx = rng.choice(n, size=n, replace=True)
            boot_means.append(np.mean(d[idx]))

    boot_means = np.array(boot_means)
    ci_lo = np.percentile(boot_means, 100 * alpha / 2)
    ci_hi = np.percentile(boot_means, 100 * (1 - alpha / 2))

    return {
        "n_bootstrap": len(boot_means),
        "confidence": confidence,
        "mean_diff": np.mean(d),
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "boot_std": np.std(boot_means),
        "uses_cluster": groups is not None,
    }
