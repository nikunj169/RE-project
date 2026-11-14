"""E11: proportional (scale) form of the drift-immune excess (Round 2, post hoc — see EXPERIMENT_LOG).

From the E06 differential fits, estimate per basin (and globally) the least-squares slope through
the origin of the layer excess on the TRACE differential rate:

    Ddelta_l = m * dbeta_Cant,l   ->   m = sum(Ddelta * dbeta_Cant) / sum(dbeta_Cant^2)

m is the fractional amount by which DIC-based accumulation (relative to the deep reference)
exceeds the TRACE prediction. Joint CIs use one cruise resample per replicate shared by all
layer frames (cruise ids global), so layers that share cruises are resampled consistently.

Usage: python -m research.experiments.e11_scale_factor [--boot 500] [--cant default|atm] [--ref primary|upperdeep|abyssal]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from research.data_io import ROOT, load_config
from research.experiments.e06_differential_trends import basin_frames
from research.models.differential_trend import _codes, design, exact_twoway

SEED = 20261114
CANT_FILES = {"default": "trace_canth.csv.gz", "atm": "trace_canth_atm.csv.gz"}
REFS = {"primary": (27.8, np.inf), "upperdeep": (27.8, 28.0), "abyssal": (28.0, np.inf)}  # E14 (GAP_AUDIT_ROUND3 §4)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=500)
    ap.add_argument("--cant", default="default", choices=list(CANT_FILES))
    ap.add_argument("--ref", default="primary", choices=list(REFS))
    args = ap.parse_args()
    sfx = ("" if args.cant == "default" else f"_{args.cant}") + ("" if args.ref == "primary" else f"_ref_{args.ref}")
    cfg = load_config()
    frames = []
    for basin, layer, st, mask in basin_frames(cfg, CANT_FILES[args.cant], ref=REFS[args.ref]):
        X = design(st, mask, ("theta", "salinity", "aou"))
        Y = st[["tco2", "canth"]].to_numpy(float)
        a = _codes(st["cell"].astype(str) + "|" + np.where(mask, "T", "R"))
        frames.append(dict(basin=basin, layer=layer, X=X, Y=Y, a=a, cruise=st["cruise"].to_numpy()))
    all_cru = np.unique(np.concatenate([f["cruise"] for f in frames]))
    for f in frames:
        f["ci"] = np.searchsorted(all_cru, f["cruise"])

    def est(wc):
        res = []
        for f in frames:
            w = np.ones(len(f["ci"])) if wc is None else wc[f["ci"]]
            k = w > 0
            aa = np.unique(f["a"][k], return_inverse=True)[1]
            bb = np.unique(f["ci"][k], return_inverse=True)[1]
            res.append(exact_twoway(f["X"][k], f["Y"][k], aa, bb, w[k])[-1])
        return np.array(res)  # layers x (dic, cant)

    def scale(r, idx):
        dd, dc = r[idx, 0] - r[idx, 1], r[idx, 1]
        return float(np.sum(dd * dc) / np.sum(dc * dc))

    point = est(None)
    rng = np.random.default_rng(SEED)
    G = len(all_cru)
    boots = [est(np.bincount(rng.integers(0, G, G), minlength=G).astype(float)) for _ in range(args.boot)]
    sets = {b: [i for i, f in enumerate(frames) if f["basin"] == b] for b in sorted({f["basin"] for f in frames})}
    sets["all"] = list(range(len(frames)))
    sets["all, gamma<=27.0"] = [i for i, f in enumerate(frames) if float(f["layer"].strip("(]").split(",")[1]) <= 27.0]
    rows = []
    for name, idx in sets.items():
        bs = np.array([scale(b, idx) for b in boots])
        rows.append(dict(set=name, n_layers=len(idx), m=scale(point, idx),
                         m_lo=np.quantile(bs, 0.025), m_hi=np.quantile(bs, 0.975)))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / f"e11_scale_factor{sfx}.csv", index=False)


if __name__ == "__main__":
    main()
