"""
Model 4: Full Second-Order Polynomial
TCO2 = b0 + b1*S + b2*T + b3*A + b4*S^2 + b5*T^2 + b6*A^2 + b7*ST + b8*SA + b9*TA

where A = AOU (raw, not scaled — matching the old benchmark).
10 free parameters (9 features + intercept).
"""

import numpy as np
from sklearn.linear_model import LinearRegression
from .base import BaseModel


def full_quad_features(S, T, A):
    """
    Construct the 9 features for the full second-order polynomial.

    Returns array with columns: S, T, A, S^2, T^2, A^2, S*T, S*A, T*A
    """
    return np.column_stack([
        S, T, A,
        S**2, T**2, A**2,
        S * T, S * A, T * A,
    ])


class FullQuadratic(BaseModel):
    name = "Full Second-Order Polynomial"
    parameter_count = 10
    effective_rank = 3
    complexity = 20
    category = "polynomial"

    def __init__(self):
        super().__init__()
        self.model_ = None

    def fit(self, S, T, A, Y):
        X = full_quad_features(
            np.asarray(S), np.asarray(T), np.asarray(A))
        Y = np.asarray(Y)
        self.model_ = LinearRegression().fit(X, Y)
        coef_names = ["S", "T", "A", "S2", "T2", "A2", "ST", "SA", "TA"]
        self.params_ = {"intercept": self.model_.intercept_}
        for name, val in zip(coef_names, self.model_.coef_):
            self.params_[f"coef_{name}"] = val
        self.is_fitted = True
        return self

    def predict(self, S, T, A):
        X = full_quad_features(
            np.asarray(S), np.asarray(T), np.asarray(A))
        return self.model_.predict(X)
