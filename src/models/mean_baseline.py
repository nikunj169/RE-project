"""
Model 0: Baseline Mean
Predict the training-set mean TCO2 for all observations.
"""

import numpy as np
from .base import BaseModel


class MeanBaseline(BaseModel):
    name = "Baseline Mean"
    parameter_count = 1
    effective_rank = 0
    complexity = 1
    category = "baseline"

    def __init__(self):
        super().__init__()
        self.mean_ = None

    def fit(self, S, T, A, Y):
        self.mean_ = np.mean(Y)
        self.params_ = {"mean": self.mean_}
        self.is_fitted = True
        return self

    def predict(self, S, T, A):
        return np.full_like(S, self.mean_, dtype=float)
