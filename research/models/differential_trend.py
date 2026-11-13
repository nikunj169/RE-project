"""Drift-immune differential trend estimator (Round 2, E06).

For one basin, a target density layer is stacked with a deep reference layer set, restricted to
cruises that sampled both. The model

    y = alpha_{cell x layer} + gamma_cruise + sum_layer f_layer(theta, S, AOU) + dbeta * (t - 2000) * 1[target] + e

is fitted with two-way fixed effects (cell x layer, cruise) by alternating projections. The cruise
effects absorb any additive per-cruise offset, including offsets that drift with time, so ``dbeta``
is the trend of the target layer *relative to* the reference layers, free of such offsets.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from research.models.within_cell_trend import T_REF, conditioning_terms


def _group_demean(M: np.ndarray, codes: np.ndarray, n: int, w: np.ndarray) -> np.ndarray:
    cnt = np.bincount(codes, weights=w, minlength=n)
    cnt = np.maximum(cnt, 1e-12)
    out = np.empty_like(M)
    for j in range(M.shape[1]):
        means = np.bincount(codes, weights=M[:, j] * w, minlength=n) / cnt
        out[:, j] = M[:, j] - means[codes]
    return out


def twoway_demean(M: np.ndarray, a: np.ndarray, b: np.ndarray, w: np.ndarray,
                  tol: float = 1e-9, max_iter: int = 500) -> np.ndarray:
    """Project out two sets of fixed effects (codes ``a`` and ``b``) by alternating projections."""
    na, nb = a.max() + 1, b.max() + 1
    scale = np.maximum(np.abs(M).max(0), 1e-12)
    cur = M.copy()
    for _ in range(max_iter):
        nxt = _group_demean(_group_demean(cur, a, na, w), b, nb, w)
        if np.max(np.abs(nxt - cur).max(0) / scale) < tol:
            return nxt
        cur = nxt
    return cur


def exact_twoway(X: np.ndarray, Y: np.ndarray, a: np.ndarray, b: np.ndarray, w: np.ndarray,
                 return_resid: bool = False):
    """Exact weighted OLS of Y on X with fixed effects for codes ``a`` (many) and ``b`` (few).

    ``a`` effects are swept out by within-group demeaning; ``b`` effects enter through the Schur
    complement of the normal equations, assembled from weighted group sums (no n x G matrix).
    Returns the X coefficients (k x m) and, optionally, residuals (n x m).
    """
    na, nb = a.max() + 1, b.max() + 1
    Z = np.hstack([X, Y])
    wa = np.maximum(np.bincount(a, weights=w, minlength=na), 1e-12)
    wb = np.bincount(b, weights=w, minlength=nb)
    A = np.column_stack([np.bincount(a, weights=Z[:, j] * w, minlength=na) for j in range(Z.shape[1])])
    Zt = Z - (A / wa[:, None])[a]
    N = sp.csr_matrix((w, (a, b)), shape=(na, nb))
    ZZ = (Zt * w[:, None]).T @ Zt
    DZ = np.column_stack([np.bincount(b, weights=Zt[:, j] * w, minlength=nb) for j in range(Z.shape[1])])
    DD = np.diag(wb) - (N.T @ sp.diags(1 / wa) @ N).toarray()
    DDp = np.linalg.pinv(DD, rcond=1e-10, hermitian=True)
    S = ZZ - DZ.T @ DDp @ DZ
    k = X.shape[1]
    coef = np.linalg.lstsq(S[:k, :k], S[:k, k:], rcond=None)[0]
    if not return_resid:
        return coef
    gamma = DDp @ (DZ[:, k:] - DZ[:, :k] @ coef)
    Dt_g = gamma[b] - ((N @ gamma) / wa[:, None])[a]
    return coef, Zt[:, k:] - Zt[:, :k] @ coef - Dt_g


def design(df: pd.DataFrame, target_mask: np.ndarray, variables: tuple[str, ...],
           form: str = "quadratic") -> np.ndarray:
    """Layer-specific conditioning terms (standardised within each part) + t*1[target] last."""
    parts = []
    for m in (target_mask, ~target_mask):
        Z = np.zeros((len(df), 0))
        if variables:
            sub = conditioning_terms(df[m], form, variables)
            Z = np.zeros((len(df), sub.shape[1]))
            Z[m] = sub
        parts.append(Z)
    t = (df["decimal_year"].to_numpy(float) - T_REF) * target_mask
    return np.hstack(parts + [t[:, None]])


def _codes(s: pd.Series) -> np.ndarray:
    return np.unique(s.to_numpy(), return_inverse=True)[1]


def fit_differential(df: pd.DataFrame, target_mask: np.ndarray, responses: list[str], *,
                     variables: tuple[str, ...] = ("theta", "salinity", "aou"), form: str = "quadratic",
                     n_boot: int = 0, seed: int = 20261114) -> dict:
    """Return point estimates and bootstrap draws of dbeta for each response (+ DIC - Cant)."""
    X = design(df, target_mask, variables, form)
    Y = df[responses].to_numpy(float)
    a = _codes(df["cell"].astype(str) + "|" + np.where(target_mask, "T", "R"))
    cru = _codes(df["cruise"])
    G = cru.max() + 1

    def solve(w: np.ndarray, keep: np.ndarray) -> np.ndarray:
        aa, bb = np.unique(a[keep], return_inverse=True)[1], np.unique(cru[keep], return_inverse=True)[1]
        return exact_twoway(X[keep], Y[keep], aa, bb, w[keep])[-1]

    ones = np.ones(len(df))
    est = solve(ones, ones > 0)
    out = {"point": dict(zip(responses, est.tolist())), "n": len(df), "n_cruises": int(G),
           "n_target": int(target_mask.sum())}
    if n_boot:
        rng = np.random.default_rng(seed)
        draws = np.full((n_boot, len(responses)), np.nan)
        for i in range(n_boot):
            w = np.bincount(rng.integers(0, G, G), minlength=G).astype(float)[cru]
            keep = w > 0
            if target_mask[keep].sum() < 50 or (~target_mask[keep]).sum() < 50:
                continue
            draws[i] = solve(w, keep)
        out["boot"] = {r: draws[:, j] for j, r in enumerate(responses)}
    return out
