"""
Model 2: Rank-1 Quadratic
TCO2 = (alpha*S + beta*T + gamma*(AOU/100) + delta)^2 + epsilon

This is the old discovered structure. It imposes a rank-1 quadratic
dependence — a one-dimensional nonlinear curve in 3D predictor space.
"""

import numpy as np
from scipy.optimize import curve_fit
from .base import BaseModel


def rank1_quad_func(X, alpha, beta, gamma, delta, epsilon):
    """The rank-1 quadratic model function."""
    S, T, A = X
    core = alpha * S + beta * T + gamma * (A / 100.0) + delta
    return core ** 2 + epsilon


class Rank1Quadratic(BaseModel):
    name = "Rank-1 Quadratic"
    parameter_count = 5
    effective_rank = 1
    complexity = 8
    category = "quadratic"

    DEFAULT_P0 = [1.0, -0.3, 2.5, 20.0, 2000.0]
    MAXFEV = 30000

    def __init__(self, p0=None, maxfev=None):
        super().__init__()
        self.p0 = p0 or self.DEFAULT_P0
        self.maxfev = maxfev or self.MAXFEV
        self.pcov_ = None

    def fit(self, S, T, A, Y):
        X = (np.asarray(S, dtype=float),
             np.asarray(T, dtype=float),
             np.asarray(A, dtype=float))
        Y = np.asarray(Y, dtype=float)

        popt, pcov = curve_fit(rank1_quad_func, X, Y,
                               p0=self.p0, maxfev=self.maxfev)
        self.params_ = {
            "alpha": popt[0],
            "beta": popt[1],
            "gamma": popt[2],
            "delta": popt[3],
            "epsilon": popt[4],
        }
        self.pcov_ = pcov
        self.is_fitted = True
        return self

    def predict(self, S, T, A):
        X = (np.asarray(S, dtype=float),
             np.asarray(T, dtype=float),
             np.asarray(A, dtype=float))
        p = self.params_
        return rank1_quad_func(X, p["alpha"], p["beta"],
                               p["gamma"], p["delta"], p["epsilon"])

    def get_core(self, S, T, A):
        """
        Compute the core variable z = alpha*S + beta*T + gamma*(AOU/100) + delta.
        """
        p = self.params_
        return (p["alpha"] * np.asarray(S) +
                p["beta"] * np.asarray(T) +
                p["gamma"] * (np.asarray(A) / 100.0) +
                p["delta"])

    def get_hessian_eigenvalue(self):
        """
        For f(x) = (v^T x + delta)^2 + epsilon, the Hessian is H = 2vv^T.
        The sole non-zero eigenvalue is 2||v||^2.
        """
        p = self.params_
        v = np.array([p["alpha"], p["beta"], p["gamma"] / 100.0])
        return 2.0 * np.dot(v, v)

    def get_coefficient_std(self):
        """Return standard errors of fitted coefficients from covariance matrix."""
        if self.pcov_ is None:
            return None
        return np.sqrt(np.diag(self.pcov_))
