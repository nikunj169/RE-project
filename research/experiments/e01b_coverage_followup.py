"""E01b (deviation, logged): coverage of percentile-bootstrap vs CR1 vs conservative-union
intervals in the E01 settings that failed the coverage criterion.

Usage: python -m research.experiments.e01b_coverage_followup [--reps 200] [--boot 200]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.experiments.e01_synthetic_recovery import null_components
from research.models.within_cell_trend import fit_trend

CASES = [("Atlantic", "(27.8, 28.0]", 0.0, 4.0), ("Atlantic", "(27.8, 28.0]", 0.1, 4.0),
         ("Pacific", "(26.0, 26.5]", 0.5, 4.0)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--boot", type=int, default=200)
    args = ap.parse_args()
    cfg = load_config()
    df = load_data(cfg)
    sel = {(b, l): g for b, l, g in groups(df, cfg)}
    rng = np.random.default_rng(777)
    rows = []
    for basin, layer, beta_true, sig in CASES:
        g = sel[(basin, layer)].copy()
        fitted, resid = null_components(g)
        t = g.decimal_year.to_numpy() - 2000
        cr, cc = np.unique(g.cruise, return_inverse=True)
        cov = {"boot": [], "cr1": [], "union": []}
        for _ in range(args.reps):
            g["y"] = fitted + beta_true * t + rng.normal(0, sig, len(cr))[cc] + rng.permutation(resid)
            f = fit_trend(g, ["y"], n_boot=args.boot, seed=int(rng.integers(1 << 31)))
            lo_b, hi_b = f.ci("y")
            lo_c, hi_c = f.beta["y"] - 1.96 * f.se_cr1["y"], f.beta["y"] + 1.96 * f.se_cr1["y"]
            cov["boot"].append(lo_b <= beta_true <= hi_b)
            cov["cr1"].append(lo_c <= beta_true <= hi_c)
            cov["union"].append(min(lo_b, lo_c) <= beta_true <= max(hi_b, hi_c))
        rows.append(dict(basin=basin, layer=layer, beta_true=beta_true, sigma_cruise=sig, reps=args.reps,
                         **{f"coverage_{k}": 100 * np.mean(v) for k, v in cov.items()}))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e01b_coverage_followup.csv", index=False)


if __name__ == "__main__":
    main()
