"""E09: pooled estimates with a global cruise bootstrap, and piecewise-linear time (GAP_AUDIT_ROUND2).

(a) One cruise resample per replicate shared by all basin x layer groups, so pooled statistics
    across groups have valid joint CIs. Pooled over "upper" (layers with upper edge <= 27.0) and
    "intermediate" (27.0 < upper edge <= 27.8) groups, sample-size weighted:
        delta_pool = sum n_g delta_g / sum n_g,   R_pool = sum n_g bDIC_g / sum n_g bCant_g
(b) Piecewise-linear time with a knot at 2000: slopes 1972-2000 and 2000-2021 for DIC and Cant.

Usage: python -m research.experiments.e09_pooled_and_timeshape [--boot 500]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.models.within_cell_trend import T_REF, conditioning_terms, solve_subset

SEED = 20261114
CANT_FILES = {"default": "trace_canth.csv.gz", "atm": "trace_canth_atm.csv.gz"}


def upper(layer: str) -> float:
    return float(layer.strip("(]").split(",")[1])


def prep(g: pd.DataFrame, piecewise: bool):
    t = g.decimal_year.to_numpy(float) - T_REF
    Xc = conditioning_terms(g, "quadratic", ("theta", "salinity", "aou"))
    T = np.column_stack([np.minimum(t, 0), np.maximum(t, 0)]) if piecewise else t[:, None]
    X = np.hstack([Xc, T])
    Y = g[["tco2", "canth"]].to_numpy(float)
    cell = np.unique(g.cell.to_numpy(), return_inverse=True)[1]
    return X, Y, cell


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=500)
    ap.add_argument("--cant", default="default", choices=list(CANT_FILES))
    args = ap.parse_args()
    sfx = "" if args.cant == "default" else f"_{args.cant}"
    cfg = load_config()
    df = load_data(cfg)
    tr = pd.read_csv(ROOT / "research" / "data" / CANT_FILES[args.cant]).set_index("row_id")
    df = df.join(tr[["canth"]], how="inner")
    df = df[np.isfinite(df.canth)]
    all_cruises = np.unique(df.cruise.to_numpy())
    G = len(all_cruises)

    grp = []
    for basin, layer, g in groups(df, cfg):
        cidx = np.searchsorted(all_cruises, g.cruise.to_numpy())
        grp.append(dict(basin=basin, layer=layer, n=len(g), cidx=cidx,
                        lin=prep(g, False), pw=prep(g, True)))

    def estimate(w_cruise: np.ndarray | None):
        out = []
        for d in grp:
            w = np.ones(d["n"]) if w_cruise is None else w_cruise[d["cidx"]]
            keep = w > 0
            res = {}
            for key in ("lin", "pw"):
                X, Y, cell = d[key]
                if keep.sum() < X.shape[1] + 10:
                    res[key] = np.full((2 if key == "pw" else 1, 2), np.nan)
                    continue
                coef = solve_subset(X[keep], Y[keep], cell[keep], w[keep])
                res[key] = coef[-1:] if key == "lin" else coef[-2:]
            out.append(res)
        return out

    rng = np.random.default_rng(SEED)
    point = estimate(None)
    boots = [estimate(np.bincount(rng.integers(0, G, G), minlength=G).astype(float)) for _ in range(args.boot)]

    # per-group piecewise slopes: pre = coef(t_pre) ... model a*min(t,0)+b*max(t,0): pre slope a, post slope b
    rows = []
    for i, d in enumerate(grp):
        pw = point[i]["pw"]
        pb = np.array([b[i]["pw"] for b in boots])  # boot x 2 (pre, post) x 2 (dic, cant)
        dpre, dpost = pb[:, 0, 0] - pb[:, 0, 1], pb[:, 1, 0] - pb[:, 1, 1]
        q = lambda x: np.nanquantile(x, [0.025, 0.975])
        rows.append(dict(basin=d["basin"], layer=d["layer"], n=d["n"],
                         dic_pre=pw[0, 0], cant_pre=pw[0, 1], delta_pre=pw[0, 0] - pw[0, 1],
                         delta_pre_lo=q(dpre)[0], delta_pre_hi=q(dpre)[1],
                         dic_post=pw[1, 0], cant_post=pw[1, 1], delta_post=pw[1, 0] - pw[1, 1],
                         delta_post_lo=q(dpost)[0], delta_post_hi=q(dpost)[1],
                         ddelta_post_minus_pre=(pw[1, 0] - pw[1, 1]) - (pw[0, 0] - pw[0, 1]),
                         ddelta_lo=q(dpost - dpre)[0], ddelta_hi=q(dpost - dpre)[1]))
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / f"e09_piecewise_time{sfx}.csv", index=False)

    # pooled
    pooled = []
    sets = {"upper (<=27.0)": lambda u: u <= 27.0, "intermediate (27.0-27.8]": lambda u: 27.0 < u <= 27.8,
            "deep (>27.8)": lambda u: u > 27.8}
    basins = sorted({d["basin"] for d in grp})
    for set_name, rule in sets.items():
        for bsel in ["all"] + basins:
            idx = [i for i, d in enumerate(grp) if rule(upper(d["layer"])) and bsel in ("all", d["basin"])]
            if not idx:
                continue
            n = np.array([grp[i]["n"] for i in idx], float)

            def pool(est):
                bd = np.array([est[i]["lin"][0, 0] for i in idx])
                bc = np.array([est[i]["lin"][0, 1] for i in idx])
                ok = np.isfinite(bd) & np.isfinite(bc)
                return (np.sum(n[ok] * (bd - bc)[ok]) / n[ok].sum(),
                        np.sum(n[ok] * bd[ok]) / np.sum(n[ok] * bc[ok]))
            p = pool(point)
            bs = np.array([pool(b) for b in boots])
            pooled.append(dict(layer_set=set_name, basin=bsel, n_groups=len(idx), n=int(n.sum()),
                               delta=p[0], delta_lo=np.nanquantile(bs[:, 0], 0.025),
                               delta_hi=np.nanquantile(bs[:, 0], 0.975), R=p[1],
                               R_lo=np.nanquantile(bs[:, 1], 0.025), R_hi=np.nanquantile(bs[:, 1], 0.975)))
            print(pooled[-1], flush=True)
    pd.DataFrame(pooled).to_csv(ROOT / "research" / "results" / f"e09_pooled{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
