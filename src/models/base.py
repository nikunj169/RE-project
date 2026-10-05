"""
Abstract base class for all model families.
"""

from abc import ABC, abstractmethod
import numpy as np


class BaseModel(ABC):
    """Abstract base class for TCO2 prediction models."""

    name: str = "BaseModel"
    parameter_count: int = 0
    effective_rank: int = 0
    complexity: int = 0
    category: str = "base"

    def __init__(self):
        self.is_fitted = False
        self.params_ = None

    @abstractmethod
    def fit(self, S, T, A, Y):
        """
        Fit the model to training data.

        Parameters
        ----------
        S, T, A : np.ndarray
            Salinity, temperature, AOU (raw, not scaled).
        Y : np.ndarray
            Target TCO2.
        """
        pass

    @abstractmethod
    def predict(self, S, T, A):
        """
        Predict TCO2 from predictors.

        Parameters
        ----------
        S, T, A : np.ndarray
            Salinity, temperature, AOU (raw, not scaled).

        Returns
        -------
        np.ndarray
            Predicted TCO2.
        """
        pass

    def get_params_dict(self):
        """Return fitted parameters as a dict."""
        return self.params_ if self.params_ else {}

    def get_description(self):
        """Return a human-readable description."""
        return {
            "name": self.name,
            "parameter_count": self.parameter_count,
            "effective_rank": self.effective_rank,
            "complexity": self.complexity,
            "category": self.category,
        }
