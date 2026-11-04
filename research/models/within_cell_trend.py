"""Within-cell fixed-effects trend estimator with cruise-cluster bootstrap.

Model (per basin x density layer):
    y_i = alpha_cell(i) + f(theta, S, AOU) + beta * (t_i - 2000) + e_i
estimated by OLS after within-cell demeaning (Frisch-Waugh-Lovell). Several response
columns can be fitted jointly so that bootstrap replicates are shared (joint CIs).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

T_REF = 2000.0


def add_cells(df: pd.DataFrame, cell_deg: float) -> pd.Series:
    """Cell id from a regular lat/lon grid (longitude wrapped to [-180, 180))."""
    lon = np.where(df["longitude"] > 180, df["longitude"] - 360, df["longitude"])
    return (np.floor(df["latitude"].to_numpy() / cell_deg).astype(int).astype(str) + "_"
            + np.floor(lon / cell_deg).astype(int).astype(str))


def select_cells(df: pd.DataFrame, min_span: float) -> pd.DataFrame:
    """Keep cells whose observations span at least ``min_span`` years."""
    span = df.groupby("cell")["decimal_year"].agg(lambda s: s.max() - s.min())
    return df[df["cell"].isin(span[span >= min_span].index)]


def conditioning_terms(df: pd.DataFrame, form: str, variables: tuple[str, ...]) -> np.ndarray:
    """Standardised conditioning covariates: 'quadratic' (incl. interactions) or 'linear'."""
    if not variables:
        return np.empty((len(df), 0))
    Z = df[list(variables)].to_numpy(float)
    Z = (Z - Z.mean(0)) / Z.std(0)
    cols = [Z[:, j] for j in range(Z.shape[1])]
    if form == "quadratic":
        k = Z.shape[1]
        cols += [Z[:, a] * Z[:, b] for a in range(k) for b in range(a, k)]
    elif form != "linear":
        raise ValueError(form)
    return np.column_stack(cols)


def _demean(M: np.ndarray, codes: np.ndarray, n_groups: int, w: np.ndarray | None = None) -> np.ndarray:
    """Subtract (weighted) group means column-wise. ``w`` = multiplicity weights."""
    w = np.ones(len(codes)) if w is None else w
    sums = np.zeros((n_groups, M.shape[1]))
    np.add.at(sums, codes, M * w[:, None])
    cnt = np.bincount(codes, weights=w, minlength=n_groups)
    means = sums / np.maximum(cnt, 1e-12)[:, None]
    return M - means[codes]


@dataclass
class TrendFit:
    beta: dict[str, float]
    n: int
    n_cells: int
    n_cruises: int
    median_span: float
    se_cr1: dict[str, float]
    boot: dict[str, np.ndarray] = field(default_factory=dict)

    def ci(self, key: str, level: float = 0.95) -> tuple[float, float]:
        a = (1 - level) / 2
        b = self.boot[key]
        b = b[np.isfinite(b)]
        return float(np.quantile(b, a)), float(np.quantile(b, 1 - a))


def fit_trend(df: pd.DataFrame, responses: list[str], *, form: str = "quadratic",
              variables: tuple[str, ...] = ("theta", "salinity", "aou"),
              extra: list[str] | None = None, n_boot: int = 0, seed: int = 20261104,
              ratio: tuple[str, str] | None = None) -> TrendFit:
    """Fit beta (time coefficient) for each response column.

    ``ratio=(num, den)`` additionally stores bootstrap distributions of num/den and num-den.
    Bootstrap resamples cruises with replacement and keeps duplicate draws (multiplicity
    weights), unlike the repository's np.isin implementation.
    """
    cells, cell_codes = np.unique(df["cell"].to_numpy(), return_inverse=True)
    cruises, cruise_codes = np.unique(df["cruise"].to_numpy(), return_inverse=True)
    Xc = conditioning_terms(df, form, variables)
    extras = [df[c].to_numpy(float)[:, None] for c in (extra or [])]
    t = (df["decimal_year"].to_numpy(float) - T_REF)[:, None]
    X = np.hstack([Xc, *extras, t]) if extras else np.hstack([Xc, t])
    Y = df[responses].to_numpy(float)

    def solve(w: np.ndarray | None) -> np.ndarray:
        Xd = _demean(X, cell_codes, len(cells), w)
        Yd = _demean(Y, cell_codes, len(cells), w)
        if w is None:
            coef, *_ = np.linalg.lstsq(Xd, Yd, rcond=None)
        else:
            sw = np.sqrt(w)[:, None]
            coef, *_ = np.linalg.lstsq(Xd * sw, Yd * sw, rcond=None)
        return coef

    coef = solve(None)
    beta = {r: float(coef[-1, j]) for j, r in enumerate(responses)}

    # CR1 cluster-robust SE (cruise) for the time coefficient
    Xd = _demean(X, cell_codes, len(cells))
    Yd = _demean(Y, cell_codes, len(cells))
    bread = np.linalg.pinv(Xd.T @ Xd)
    se = {}
    G, n, k = len(cruises), len(df), X.shape[1] + len(cells)
    for j, r in enumerate(responses):
        u = Yd[:, j] - Xd @ coef[:, j]
        S = np.zeros((X.shape[1], X.shape[1]))
        sc = np.zeros((G, X.shape[1]))
        np.add.at(sc, cruise_codes, Xd * u[:, None])
        S = sc.T @ sc
        V = bread @ S @ bread * (G / (G - 1)) * ((n - 1) / max(n - k, 1))
        se[r] = float(np.sqrt(V[-1, -1]))

    span = df.groupby("cell")["decimal_year"].agg(lambda s: s.max() - s.min())
    fit = TrendFit(beta=beta, n=len(df), n_cells=len(cells), n_cruises=G,
                   median_span=float(span.median()), se_cr1=se)

    if n_boot:
        rng = np.random.default_rng(seed)
        draws = {r: np.full(n_boot, np.nan) for r in responses}
        for b in range(n_boot):
            counts = np.bincount(rng.integers(0, G, G), minlength=G).astype(float)
            w = counts[cruise_codes]
            keep = w > 0
            if keep.sum() < X.shape[1] + 10:
                continue
            # cells absent from the replicate contribute nothing (zero weight)
            cb = solve_subset(X[keep], Y[keep], cell_codes[keep], w[keep])
            for j, r in enumerate(responses):
                draws[r][b] = cb[-1, j]
        fit.boot = draws
        if ratio:
            num, den = ratio
            fit.boot[f"{num}/{den}"] = draws[num] / draws[den]
            fit.boot[f"{num}-{den}"] = draws[num] - draws[den]
    return fit


def solve_subset(X: np.ndarray, Y: np.ndarray, codes: np.ndarray, w: np.ndarray) -> np.ndarray:
    uniq, c = np.unique(codes, return_inverse=True)
    Xd = _demean(X, c, len(uniq), w)
    Yd = _demean(Y, c, len(uniq), w)
    sw = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(Xd * sw, Yd * sw, rcond=None)
    return coef
