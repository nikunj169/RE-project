"""
Model 1: Linear Regression
TCO2 = b0 + b1*S + b2*T + b3*A

where A = AOU/100 (consistently scaled).
"""

import numpy as np
from sklearn.linear_model import LinearRegression
from .base import BaseModel


class LinearModel(BaseModel):
    name = "Linear Regression"
    parameter_count = 4
    effective_rank = 0
    complexity = 4
    category = "linear"

    AOU_DIVISOR = 100.0

    def __init__(self):
        super().__init__()
        self.model_ = None

    def _features(self, S, T, A):
        return np.column_stack([S, T, A / self.AOU_DIVISOR])

    def fit(self, S, T, A, Y):
        X = self._features(S, T, A)
        self.model_ = LinearRegression().fit(X, Y)
        self.params_ = {
            "intercept": self.model_.intercept_,
            "coef_S": self.model_.coef_[0],
            "coef_T": self.model_.coef_[1],
            "coef_A": self.model_.coef_[2],
        }
        self.is_fitted = True
        return self

    def predict(self, S, T, A):
        X = self._features(S, T, A)
        return self.model_.predict(X)
