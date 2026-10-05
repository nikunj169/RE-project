"""
Model 3: Rank-2 Quadratic
TCO2 = z1^2 + z2^2 + epsilon

where:
    z1 = a1*S + b1*T + c1*(AOU/100) + d1
    z2 = a2*S + b2*T + c2*(AOU/100) + d2

This is parameterized as a sum of two squared linear forms.
Total parameters: 2*(3+1) + 1 = 9.

Identifiability note: The decomposition is not unique (rotation ambiguity),
but the resulting model function is unique.
"""

import numpy as np
from scipy.optimize import minimize
from .base import BaseModel


def rank2_quad_predict(S, T, A, params):
    """
    Predict using rank-2 quadratic model.

    params = [a1, b1, c1, d1, a2, b2, c2, d2, epsilon]
    """
    a1, b1, c1, d1, a2, b2, c2, d2, eps = params
    z1 = a1 * S + b1 * T + c1 * (A / 100.0) + d1
    z2 = a2 * S + b2 * T + c2 * (A / 100.0) + d2
    return z1**2 + z2**2 + eps


class Rank2Quadratic(BaseModel):
    name = "Rank-2 Quadratic"
    parameter_count = 9
    effective_rank = 2
    complexity = 16
    category = "quadratic"

    def __init__(self, n_restarts=5, random_state=42):
        super().__init__()
        self.n_restarts = n_restarts
        self.random_state = random_state

    def fit(self, S, T, A, Y):
        S = np.asarray(S, dtype=float)
        T = np.asarray(T, dtype=float)
        A = np.asarray(A, dtype=float)
        Y = np.asarray(Y, dtype=float)

        def objective(params):
            pred = rank2_quad_predict(S, T, A, params)
            return np.mean((Y - pred) ** 2)

        rng = np.random.RandomState(self.random_state)
        best_result = None
        best_loss = np.inf

        # Initialize from rank-1 solution as one starting point
        from .rank1_quadratic import Rank1Quadratic
        r1 = Rank1Quadratic()
        r1.fit(S, T, A, Y)
        p = r1.params_
        p0_r1 = [p["alpha"], p["beta"], p["gamma"], p["delta"],
                  0.1, 0.01, 0.05, 0.0, p["epsilon"]]

        all_starts = [p0_r1]
        for _ in range(self.n_restarts - 1):
            noise = rng.randn(9) * 0.5
            all_starts.append(np.array(p0_r1) + noise)

        for p0 in all_starts:
            try:
                result = minimize(objective, p0, method="Nelder-Mead",
                                  options={"maxiter": 50000, "xatol": 1e-8,
                                           "fatol": 1e-8})
                if result.fun < best_loss:
                    best_loss = result.fun
                    best_result = result
            except Exception:
                continue

        if best_result is None:
            raise RuntimeError("Rank-2 quadratic fitting failed for all restarts")

        p_opt = best_result.x
        self.params_ = {
            "a1": p_opt[0], "b1": p_opt[1], "c1": p_opt[2], "d1": p_opt[3],
            "a2": p_opt[4], "b2": p_opt[5], "c2": p_opt[6], "d2": p_opt[7],
            "epsilon": p_opt[8],
        }
        self.is_fitted = True
        return self

    def predict(self, S, T, A):
        p = self.params_
        params = [p["a1"], p["b1"], p["c1"], p["d1"],
                  p["a2"], p["b2"], p["c2"], p["d2"], p["epsilon"]]
        return rank2_quad_predict(
            np.asarray(S), np.asarray(T), np.asarray(A), params)
