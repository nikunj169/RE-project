"""
Model 6+: PySR Symbolic Regression Candidates

Wraps PySR to produce Pareto-optimal candidates and records the
complete candidate library for analysis.
"""

import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)

try:
    from pysr import PySRRegressor
    PYSR_AVAILABLE = True
except ImportError:
    PYSR_AVAILABLE = False
    logger.warning("PySR not installed. Install with: pip install pysr && python -m pysr install")


class SymbolicRegressionSearch:
    """
    Run PySR symbolic regression and extract the Pareto-optimal candidate library.
    """

    def __init__(self, search_config, random_state=42):
        """
        Parameters
        ----------
        search_config : dict
            Search space configuration (from pysr_config.yaml).
        random_state : int
            Random seed.
        """
        self.search_config = search_config
        self.random_state = random_state
        self.model_ = None
        self.candidates_df_ = None

    def fit(self, S, T, A, Y, subsample_max=None):
        """
        Run PySR symbolic regression.

        Parameters
        ----------
        S, T, A, Y : np.ndarray
            Training data (raw, not scaled).
        subsample_max : int, optional
            Subsample training data if too large.

        Returns
        -------
        self
        """
        if not PYSR_AVAILABLE:
            raise ImportError("PySR not installed")

        S = np.asarray(S, dtype=float)
        T = np.asarray(T, dtype=float)
        A = np.asarray(A, dtype=float)
        Y = np.asarray(Y, dtype=float)

        # Subsample if needed
        n = len(Y)
        if subsample_max and n > subsample_max:
            rng = np.random.RandomState(self.random_state)
            idx = rng.choice(n, subsample_max, replace=False)
            S, T, A, Y = S[idx], T[idx], A[idx], Y[idx]
            logger.info(f"Subsampled to {subsample_max} training observations")

        X = np.column_stack([S, T, A / 100.0])

        cfg = self.search_config
        common = cfg.get("common", {})

        self.model_ = PySRRegressor(
            binary_operators=cfg.get("binary_operators", ["+", "-", "*", "/"]),
            unary_operators=cfg.get("unary_operators", ["exp", "log"]),
            maxsize=cfg.get("maxsize", 30),
            populations=cfg.get("populations", 30),
            population_size=cfg.get("population_size", 50),
            niterations=cfg.get("niterations", 100),
            parsimony=cfg.get("parsimony", 0.003),
            constraints=cfg.get("constraints", {}),
            model_selection=common.get("model_selection", "best"),
            batching=common.get("batching", True),
            batch_size=common.get("batch_size", 500),
            turbo=common.get("turbo", False),
            precision=common.get("precision", 64),
            fast_cycle=common.get("fast_cycle", True),
            random_state=self.random_state,
            deterministic=True,
            parallelism="serial",
            temp_equation_file=True,
            verbosity=0,
        )

        logger.info(f"Running PySR with seed={self.random_state}, "
                     f"operators={self.model_.binary_operators + self.model_.unary_operators}")

        self.model_.fit(X, Y, variable_names=["sal", "tmp", "aou100"])

        # Extract full candidate library
        self._extract_candidates()

        return self

    def _extract_candidates(self):
        """Extract all Pareto-optimal candidates into a DataFrame."""
        eqs = self.model_.equations_
        rows = []
        for idx, row in eqs.iterrows():
            entry = {
                "complexity": int(row["complexity"]),
                "loss": float(row["loss"]),
                "score": float(row["score"]) if "score" in row.index else np.nan,
                "equation": str(row["equation"]),
            }
            if "sympy_format" in row.index:
                entry["sympy_format"] = str(row["sympy_format"])
            rows.append(entry)
        self.candidates_df_ = pd.DataFrame(rows)

    def get_candidates(self):
        """Return the full candidate library as a DataFrame."""
        return self.candidates_df_

    def predict_from_equation(self, equation_idx, S, T, A):
        """Predict using a specific equation from the hall of fame."""
        X = np.column_stack([
            np.asarray(S, dtype=float),
            np.asarray(T, dtype=float),
            np.asarray(A, dtype=float) / 100.0,
        ])
        return self.model_.predict(X, index=equation_idx)

    def get_best_equation(self):
        """Return the best equation (complexity-accuracy tradeoff)."""
        return str(self.model_.get_best())
