"""E01: synthetic recovery test of the within-cell trend estimator (methodology section 9).

For representative basin x layer groups, the observed DIC is replaced by
    y = fitted(alpha_cell + f(theta,S,AOU)) + beta_true*(t-2000) + cruise offset + resampled residual
and the estimator (incl. cruise-cluster bootstrap) is re-applied.

Pass criteria: |mean bias| < 0.02 umol/kg/yr; 95% bootstrap CI coverage in [90, 99] %.

Usage: python -m research.experiments.e01_synthetic_recovery [--reps 100] [--boot 200]
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.models.within_cell_trend import (_demean, conditioning_terms, fit_trend)

TARGETS = [("Pacific", "(26.0, 26.5]"), ("Atlantic", "(27.8, 28.0]"),
           ("Southern Ocean", "(28.0, 28.1]"), ("Indian", "(27.0, 27.5]")]


def null_components(g: pd.DataFrame):
    """Fitted cell+conditioning part and residuals of the no-trend model on observed DIC."""
    cells, codes = np.unique(g["cell"], return_inverse=True)
    X = conditioning_terms(g, "quadratic", ("theta", "salinity", "aou"))
    y = g["tco2"].to_numpy(float)
    Xd = _demean(X, codes, len(cells))
    yd = _demean(y[:, None], codes, len(cells))[:, 0]
    coef, *_ = np.linalg.lstsq(Xd, yd, rcond=None)
    resid = yd - Xd @ coef
    return y - resid, resid


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--boot", type=int, default=200)
    args = ap.parse_args()
    cfg = load_config()
    df = load_data(cfg)
    sel = {(b, l): g for b, l, g in groups(df, cfg) if (b, l) in TARGETS}
    rng = np.random.default_rng(cfg["bootstrap"]["seed"])
    rows = []
    for (basin, layer), g in sel.items():
        g = g.copy()
        fitted, resid = null_components(g)
        t = g["decimal_year"].to_numpy() - 2000
        cruises, ccodes = np.unique(g["cruise"], return_inverse=True)
        for beta_true in (0.0, 0.1, 0.5, 1.0):
            for sigma_c in (0.0, 2.0, 4.0):
                est, cover = [], []
                for r in range(args.reps):
                    off = rng.normal(0, sigma_c, len(cruises))[ccodes]
                    g["y_syn"] = fitted + beta_true * t + off + rng.permutation(resid)
                    fit = fit_trend(g, ["y_syn"], n_boot=args.boot, seed=int(rng.integers(1 << 31)))
                    lo, hi = fit.ci("y_syn")
                    est.append(fit.beta["y_syn"])
                    cover.append(lo <= beta_true <= hi)
                rows.append(dict(basin=basin, layer=layer, beta_true=beta_true, sigma_cruise=sigma_c,
                                 mean_est=np.mean(est), bias=np.mean(est) - beta_true,
                                 sd_est=np.std(est), coverage_pct=100 * np.mean(cover), reps=args.reps))
                print(rows[-1], flush=True)
    out = pd.DataFrame(rows)
    out_dir = ROOT / "research" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / "e01_synthetic_recovery.csv", index=False)
    verdict = {"max_abs_bias": float(out.bias.abs().max()),
               "coverage_min": float(out.coverage_pct.min()), "coverage_max": float(out.coverage_pct.max()),
               "pass_bias": bool(out.bias.abs().max() < 0.02),
               "pass_coverage": bool(out.coverage_pct.between(90, 99).all())}
    (out_dir / "e01_verdict.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
