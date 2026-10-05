"""Frozen symbolic candidate adapter for leakage-safe validation."""

from __future__ import annotations

import re

import numpy as np
import sympy as sp


class FrozenSymbolicModel:
    """Evaluate a preselected PySR expression without refitting it."""

    parameter_count = None
    effective_rank = None
    complexity = None
    category = "symbolic"

    def __init__(self, equation, sympy_format=None, complexity=None):
        self.equation = str(equation)
        self.sympy_format = str(sympy_format or equation)
        self.complexity = complexity
        self.is_fitted = False
        self._compiled = None

    def fit(self, S, T, A, Y):
        """Validate the frozen expression and leave its constants unchanged."""
        if self._compiled is None:
            self._compiled = self._compile()
        self.is_fitted = True
        return self

    def _compile(self):
        sal, tmp, aou100 = sp.symbols("sal tmp aou100")
        expression = self.sympy_format
        expression = re.sub(r"\bcbrt\(", "real_cbrt(", expression)
        expression = re.sub(r"\bcube\(", "cube(", expression)
        expression = sp.sympify(
            expression,
            locals={
                "sal": sal,
                "tmp": tmp,
                "aou100": aou100,
                "real_cbrt": lambda x: sp.real_root(x, 3),
                "cube": lambda x: x**3,
            },
        )
        return sp.lambdify(
            (sal, tmp, aou100), expression, modules=["numpy"]
        )

    def predict(self, S, T, A):
        if self._compiled is None:
            self._compiled = self._compile()
        prediction = self._compiled(
            np.asarray(S, dtype=float),
            np.asarray(T, dtype=float),
            np.asarray(A, dtype=float) / 100.0,
        )
        return np.asarray(prediction, dtype=float)
