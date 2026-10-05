"""
Derivative / sensitivity analysis for the rank-1 quadratic model.

Computes and evaluates partial derivatives of TCO2 with respect to S, T, AOU
across the observed predictor domain.
"""

import numpy as np
from scipy import stats


def compute_derivatives(S, T, A, alpha, beta, gamma, delta, aou_divisor=100.0):
    """
    Compute partial derivatives of the rank-1 quadratic model.

    For TCO2 = (alpha*S + beta*T + gamma*(A/div) + delta)^2 + epsilon:

        dTCO2/dS = 2*alpha*z
        dTCO2/dT = 2*beta*z
        dTCO2/dA = 2*gamma*(1/div)*z

    where z = alpha*S + beta*T + gamma*(A/div) + delta.

    Parameters
    ----------
    S, T, A : np.ndarray
    alpha, beta, gamma, delta : float
    aou_divisor : float

    Returns
    -------
    dict with arrays: core (z), dS, dT, dA
    """
    S = np.asarray(S, dtype=float)
    T = np.asarray(T, dtype=float)
    A = np.asarray(A, dtype=float)

    core = alpha * S + beta * T + gamma * (A / aou_divisor) + delta
    dS = 2 * alpha * core
    dT = 2 * beta * core
    dA = 2 * gamma * (1.0 / aou_divisor) * core

    return {"core": core, "dS": dS, "dT": dT, "dA": dA}


def derivative_summary(S, T, A, alpha, beta, gamma, delta, aou_divisor=100.0):
    """
    Compute summary statistics for the derivatives across the observed domain.

    Returns
    -------
    dict with summary statistics.
    """
    derivs = compute_derivatives(S, T, A, alpha, beta, gamma, delta, aou_divisor)

    core = derivs["core"]
    dS = derivs["dS"]
    dT = derivs["dT"]
    dA = derivs["dA"]

    # Fraction with expected physical signs
    frac_dS_pos = np.mean(dS > 0) * 100
    frac_dT_neg = np.mean(dT < 0) * 100
    frac_dA_pos = np.mean(dA > 0) * 100

    return {
        "core_median": np.median(core),
        "core_min": np.min(core),
        "core_max": np.max(core),
        "frac_dS_positive_pct": frac_dS_pos,
        "frac_dT_negative_pct": frac_dT_neg,
        "frac_dA_positive_pct": frac_dA_pos,
        "dS_median": np.median(dS),
        "dT_median": np.median(dT),
        "dA_median": np.median(dA),
        "dS_q25": np.percentile(dS, 25),
        "dS_q75": np.percentile(dS, 75),
        "dT_q25": np.percentile(dT, 25),
        "dT_q75": np.percentile(dT, 75),
        "dA_q25": np.percentile(dA, 25),
        "dA_q75": np.percentile(dA, 75),
        "core_positive_fraction": np.mean(core > 0) * 100,
    }
