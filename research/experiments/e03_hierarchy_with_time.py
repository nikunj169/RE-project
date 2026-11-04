"""E03 (H3): the repository's model hierarchy with and without a time/Cant covariate,
under temporal, spatial-block and cruise-block validation.

Independent re-implementation (does not import repository code):
  mean, linear, rank1, rank2, full_quadratic, cubic   (S, T in situ, AOU/100 as in repo)
  each optionally + kappa*(t-2000)  [covariate 'year'] or + kappa*Cant_TRACE [covariate 'canth']

Outputs research/results/e03_*.csv.

Usage: python -m research.experiments.e03_hierarchy_with_time [--qc relaxed|strict] [--covariates none year]
"""

from __future__ import annotations

import argparse
from itertools import combinations_with_replacement

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from sklearn.model_selection import GroupKFold

from research.data_io import ROOT, load_config

BASINS = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]


# ---------------------------------------------------------------- features
def raw_X(d: pd.DataFrame) -> np.ndarray:
    return np.column_stack([d.salinity, d.temperature, d.aou / 100.0])


def poly(Z: np.ndarray, degree: int) -> np.ndarray:
    cols = [np.ones(len(Z))]
    for deg in range(1, degree + 1):
        for comb in combinations_with_replacement(range(Z.shape[1]), deg):
            cols.append(np.prod(Z[:, comb], axis=1))
    return np.column_stack(cols)


class Model:
    """Fit/predict wrapper; covariate c (n,) enters linearly if given."""

    def __init__(self, kind: str):
        self.kind = kind

    def fit(self, X, y, c=None):
        self.mu, self.sd = X.mean(0), X.std(0)
        Z = (X - self.mu) / self.sd
        self.use_c = c is not None
        k = self.kind
        if k == "mean":
            A = np.ones((len(y), 1))
        elif k in ("linear", "full_quadratic", "cubic"):
            A = poly(Z, {"linear": 1, "full_quadratic": 2, "cubic": 3}[k])
        if k in ("mean", "linear", "full_quadratic", "cubic"):
            if self.use_c:
                A = np.column_stack([A, c])
            self.coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            return self
        # rank-1 / rank-2: nonlinear least squares on standardised predictors
        self.rank = 1 if k == "rank1" else 2
        lin = np.linalg.lstsq(np.column_stack([np.ones(len(y)), Z] + ([c] if self.use_c else [])), y, rcond=None)[0]
        # initial rank-1: (v.z + d)^2 + e with v along linear gradient direction
        g = lin[1:4]
        s = np.sqrt(np.abs(g).sum()) + 1e-6
        p0 = list(g / (2 * s * 10)) + [10.0, y.mean() - 100.0]
        if self.rank == 2:
            p0 = p0[:4] + [0.01, 0.01, 0.01, 0.0] + [p0[4]]
        if self.use_c:
            p0 = p0 + [lin[-1]]
        res = least_squares(lambda p: self._rank_pred(Z, p, c) - y, np.array(p0, float),
                            method="trf", x_scale="jac", max_nfev=2000)
        self.p = res.x
        return self

    def _rank_pred(self, Z, p, c):
        if self.rank == 1:
            out = (Z @ p[0:3] + p[3]) ** 2 + p[4]
            j = 5
        else:
            out = (Z @ p[0:3] + p[3]) ** 2 + (Z @ p[4:7] + p[7]) ** 2 + p[8]
            j = 9
        if c is not None:
            out = out + p[j] * c
        return out

    def predict(self, X, c=None):
        Z = (X - self.mu) / self.sd
        k = self.kind
        if k in ("rank1", "rank2"):
            return self._rank_pred(Z, self.p, c if self.use_c else None)
        if k == "mean":
            A = np.ones((len(X), 1))
        else:
            A = poly(Z, {"linear": 1, "full_quadratic": 2, "cubic": 3}[k])
        if self.use_c:
            A = np.column_stack([A, c])
        return A @ self.coef


KINDS = ["mean", "linear", "rank1", "rank2", "full_quadratic", "cubic"]


def covariate(d: pd.DataFrame, name: str):
    if name == "none":
        return None
    if name == "year":
        return d.decimal_year.to_numpy() - 2000.0
    if name == "canth":
        return d.canth.to_numpy()
    raise ValueError(name)


def metrics(y, yhat):
    r = y - yhat
    return dict(RMSE=float(np.sqrt(np.mean(r ** 2))), Bias=float(np.mean(r)), MAE=float(np.mean(np.abs(r))))


def cruise_paired_ci(y, a, b, cruise, n_boot, rng):
    """Correct cruise-cluster bootstrap of mean(|err_b| - |err_a|) (positive -> a better)."""
    d = np.abs(y - b) - np.abs(y - a)
    cr, codes = np.unique(cruise, return_inverse=True)
    s = np.bincount(codes, weights=d, minlength=len(cr))
    n = np.bincount(codes, minlength=len(cr)).astype(float)
    draws = []
    for _ in range(n_boot):
        k = np.bincount(rng.integers(0, len(cr), len(cr)), minlength=len(cr))
        draws.append((k * s).sum() / (k * n).sum())
    return float(d.mean()), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", default="relaxed", choices=["relaxed", "strict"])
    ap.add_argument("--covariates", nargs="+", default=["none", "year"])
    ap.add_argument("--boot", type=int, default=1000)
    args = ap.parse_args()
    cfg = load_config()
    df = pd.read_csv(ROOT / cfg["data"]["path"])
    df = df[df.basin.isin(BASINS)]
    if args.qc == "strict":
        df = df[df.qc_strict]
    if "canth" in args.covariates:
        tr = pd.read_csv(ROOT / "research" / "data" / "trace_canth.csv.gz").set_index("row_id")
        df = df.join(tr[["canth"]], how="inner").dropna(subset=["canth"])
    rng = np.random.default_rng(cfg["bootstrap"]["seed"])
    temporal_rows, blocked_rows, paired_rows = [], [], []

    for basin in BASINS:
        d = df[df.basin == basin]
        tr_, va, ex = d[d.year < 2015], d[(d.year >= 2015) & (d.year < 2018)], d[d.year >= 2018]
        for cov in args.covariates:
            preds = {}
            for kind in KINDS:
                m = Model(kind).fit(raw_X(tr_), tr_.tco2.to_numpy(), covariate(tr_, cov))
                for split, part in (("train", tr_), ("validation", va), ("external", ex)):
                    yhat = m.predict(raw_X(part), covariate(part, cov))
                    preds[(kind, split)] = yhat
                    temporal_rows.append(dict(basin=basin, covariate=cov, model=kind, split=split,
                                              n=len(part), **metrics(part.tco2.to_numpy(), yhat)))
            for split, part in (("validation", va), ("external", ex)):
                for a, b in (("rank1", "full_quadratic"), ("rank1", "rank2"), ("rank1", "cubic")):
                    md, lo, hi = cruise_paired_ci(part.tco2.to_numpy(), preds[(a, split)], preds[(b, split)],
                                                  part.cruise.to_numpy(), args.boot, rng)
                    paired_rows.append(dict(basin=basin, covariate=cov, split=split, a=a, b=b,
                                            n_cruises=part.cruise.nunique(), mean_AE_diff=md,
                                            ci_lo=lo, ci_hi=hi, equivalent_pm4=bool(lo > -4 and hi < 4)))
            # blocked CV on pre-2018 data
            pre = d[d.year < 2018].reset_index(drop=True)
            lon = np.where(pre.longitude > 180, pre.longitude - 360, pre.longitude)
            block_ids = {"spatial": (np.floor(pre.latitude / 5).astype(int).astype(str) + "_"
                                     + np.floor(lon / 10).astype(int).astype(str)).to_numpy(),
                         "cruise": pre.cruise.astype(int).astype(str).to_numpy()}
            for btype, grp in block_ids.items():
                for fold, (i_tr, i_te) in enumerate(GroupKFold(5).split(pre, groups=grp)):
                    a_, b_ = pre.iloc[i_tr], pre.iloc[i_te]
                    for kind in KINDS:
                        m = Model(kind).fit(raw_X(a_), a_.tco2.to_numpy(), covariate(a_, cov))
                        yhat = m.predict(raw_X(b_), covariate(b_, cov))
                        blocked_rows.append(dict(basin=basin, covariate=cov, block=btype, fold=fold,
                                                 model=kind, n=len(b_), **metrics(b_.tco2.to_numpy(), yhat)))
            print(f"{basin} {cov} done", flush=True)

    out = ROOT / "research" / "results"
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{args.qc}_{'-'.join(args.covariates)}"
    pd.DataFrame(temporal_rows).to_csv(out / f"e03_temporal_{tag}.csv", index=False)
    b = pd.DataFrame(blocked_rows)
    b.to_csv(out / f"e03_blocked_folds_{tag}.csv", index=False)
    b.groupby(["basin", "covariate", "block", "model"]).agg(RMSE_mean=("RMSE", "mean"), RMSE_sd=("RMSE", "std"),
                                                             Bias_mean=("Bias", "mean")).reset_index() \
        .to_csv(out / f"e03_blocked_summary_{tag}.csv", index=False)
    pd.DataFrame(paired_rows).to_csv(out / f"e03_paired_{tag}.csv", index=False)


if __name__ == "__main__":
    main()
