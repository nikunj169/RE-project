"""
Model 5: Cubic Polynomial Benchmark
Full third-order polynomial with all cross-terms.

This includes all terms up to degree 3 in S, T, A (where A = AOU).
"""

import numpy as np
from sklearn.linear_model import LinearRegression
from .base import BaseModel


def cubic_features(S, T, A):
    """
    Construct features for a full cubic polynomial.

    Degree 1: S, T, A                          (3 terms)
    Degree 2: S^2, T^2, A^2, ST, SA, TA       (6 terms)
    Degree 3: S^3, T^3, A^3, S^2T, S^2A, T^2S, T^2A, A^2S, A^2T, STA (10 terms)
    Total: 19 terms + intercept = 20 parameters
    """
    return np.column_stack([
        # Degree 1
        S, T, A,
        # Degree 2
        S**2, T**2, A**2, S*T, S*A, T*A,
        # Degree 3
        S**3, T**3, A**3,
        S**2 * T, S**2 * A,
        T**2 * S, T**2 * A,
        A**2 * S, A**2 * T,
        S * T * A,
    ])


class CubicPolynomial(BaseModel):
    name = "Cubic Polynomial"
    parameter_count = 20  # 19 features + intercept
    effective_rank = None
    complexity = 30
    category = "polynomial"

    def __init__(self):
        super().__init__()
        self.model_ = None

    def fit(self, S, T, A, Y):
        X = cubic_features(
            np.asarray(S, dtype=float),
            np.asarray(T, dtype=float),
            np.asarray(A, dtype=float))
        Y = np.asarray(Y, dtype=float)

        self.model_ = LinearRegression().fit(X, Y)
        self.params_ = {"intercept": self.model_.intercept_}
        for i, val in enumerate(self.model_.coef_):
            self.params_[f"coef_{i}"] = val
        self.is_fitted = True
        return self

    def predict(self, S, T, A):
        X = cubic_features(
            np.asarray(S, dtype=float),
            np.asarray(T, dtype=float),
            np.asarray(A, dtype=float))
        return self.model_.predict(X)
