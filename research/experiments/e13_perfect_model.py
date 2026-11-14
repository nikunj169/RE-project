"""E13/E15 step 2: perfect-model tests in an ESM sampled like GLODAP (GAP_AUDIT_ROUND3 §1-2, §8).

Model data (from e13_sample_model) replace the observed theta, S, in-situ T, AOU and DIC at every
E02 sample; basin, layer, cell, cruise and dates stay observed. Responses:
  tco2         model DIC
  canth_true   model dissic - dissicnat (exact anthropogenic carbon)
  canth_trace  TRACE (default settings) run on model S, in-situ T at the samples
Tests (one global cruise resample per replicate, shared by all groups/frames):
  (a) estimator:   tco2 vs canth_true        H13a: pooled upper R in [0.95, 1.05], |delta| <= 0.03
  (b) TRACE:       canth_true vs canth_trace (no directional prediction)
  (c) model analogue of the observed comparison: tco2 vs canth_trace
Single-layer pooled form as E09; differential (cruise FE vs gamma > 27.8) scale factor as E11.

--model selects the sampled model (E13: mri_esm2_0; E15: noresm2_lm); --max-year restricts the samples
(E15 like-for-like: 2014) before the inclusion rules are applied.

Usage: python -m research.experiments.e13_perfect_model [--boot 500] [--trace-variant default|atm]
       [--model mri_esm2_0|noresm2_lm] [--max-year 2014]
"""

from __future__ import annotations

import argparse
import contextlib

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.experiments.e02a_run_trace import run_trace
from research.experiments.e06_differential_trends import basin_frames, upper
from research.experiments.e12_trace_unadjusted import unadjusted_trajectory
from research.models.differential_trend import _codes, design, exact_twoway
from research.models.within_cell_trend import T_REF, conditioning_terms, solve_subset

SEED = 20261114
DATA = ROOT / "research" / "data"
RES = ROOT / "research" / "results"
SFX = ""  # set from --model, --max-year, --trace-variant
TRACE_CACHE = {"mri_esm2_0": "trace_canth_model_mri"}  # E13 cache name; other models use trace_canth_model_<name>
RESP = ["tco2", "canth_true", "canth_trace"]
PAIRS = {"a_dic_vs_true": ("tco2", "canth_true"), "b_true_vs_trace": ("canth_true", "canth_trace"),
         "c_dic_vs_trace": ("tco2", "canth_trace")}
SETS = {"upper (<=27.0)": lambda u: u <= 27.0, "intermediate (27.0-27.8]": lambda u: 27.0 < u <= 27.8,
        "deep (>27.8)": lambda u: u > 27.8}


def model_frame(cfg, trace_variant: str = "default", model: str = "mri_esm2_0",
                max_year: int | None = None) -> pd.DataFrame:
    df = load_data(cfg)
    m = pd.read_csv(DATA / f"model_{model}.csv.gz").set_index("row_id")
    m = m[np.isfinite(m[["salinity", "temperature"]]).all(axis=1)]  # NorESM2-LM: no samples after 2014
    stem = TRACE_CACHE.get(model, f"trace_canth_model_{model}")
    tf = DATA / (f"{stem}.csv.gz" if trace_variant == "default" else f"{stem}_{trace_variant}.csv.gz")
    if not tf.exists():
        lon = np.mod(df.loc[m.index, "longitude"].to_numpy(float) + 180, 360) - 180
        g = df.loc[m.index]
        with (unadjusted_trajectory() if trace_variant == "atm" else contextlib.nullcontext()):
            out = run_trace(lon, g.latitude.to_numpy(float), g.depth.to_numpy(float), m.salinity.to_numpy(float),
                            m.temperature.to_numpy(float), g.decimal_year.to_numpy(float), cfg)
        pd.DataFrame({"row_id": m.index, "canth": np.asarray(out["canth"], float)}).to_csv(tf, index=False)
    tr = pd.read_csv(tf).set_index("row_id")
    obs_tr = pd.read_csv(DATA / "trace_canth.csv.gz").set_index("row_id")
    df = df.loc[m.index].copy()
    df["tco2_obs"], df["canth_trace_obs"] = df.tco2, obs_tr.loc[m.index, "canth"]
    for c in ("theta", "salinity", "temperature", "aou", "tco2"):
        df[c] = m[c]
    df["canth_true"], df["canth_trace"] = m.canth_true, tr.loc[m.index, "canth"]
    df["canth"] = df.canth_true  # basin_frames filters on canth
    ok = np.isfinite(df[["theta", "salinity", "aou"] + RESP]).all(axis=1)
    if max_year is not None:
        ok &= df.year <= max_year
    print(f"model rows {len(df)}, used {int(ok.sum())}", flush=True)
    return df[ok]


def boot_weights(G, n_boot):
    rng = np.random.default_rng(SEED)
    return [np.bincount(rng.integers(0, G, G), minlength=G).astype(float) for _ in range(n_boot)]


def q(x):
    return np.nanquantile(x, [0.025, 0.975])


def single_layer(df, cfg, n_boot):
    cru = np.unique(df.cruise.to_numpy())
    grp = []
    for basin, layer, g in groups(df, cfg):
        X = np.hstack([conditioning_terms(g, "quadratic", ("theta", "salinity", "aou")),
                       (g.decimal_year.to_numpy(float) - T_REF)[:, None]])
        grp.append(dict(basin=basin, layer=layer, n=len(g), X=X, Y=g[RESP].to_numpy(float),
                        cell=np.unique(g.cell.to_numpy(), return_inverse=True)[1],
                        ci=np.searchsorted(cru, g.cruise.to_numpy())))

    def est(wc):
        out = np.full((len(grp), len(RESP)), np.nan)
        for i, d in enumerate(grp):
            w = np.ones(d["n"]) if wc is None else wc[d["ci"]]
            k = w > 0
            if k.sum() >= d["X"].shape[1] + 10:
                out[i] = solve_subset(d["X"][k], d["Y"][k], d["cell"][k], w[k])[-1]
        return out

    point, boots = est(None), [est(w) for w in boot_weights(len(cru), n_boot)]
    B = np.array(boots)
    j = {r: i for i, r in enumerate(RESP)}
    rows = []
    for i, d in enumerate(grp):
        rec = dict(basin=d["basin"], layer=d["layer"], n=d["n"], **{f"beta_{r}": point[i, j[r]] for r in RESP})
        for name, (y, x) in PAIRS.items():
            dl = B[:, i, j[y]] - B[:, i, j[x]]
            rec[f"delta_{name}"] = point[i, j[y]] - point[i, j[x]]
            rec[f"delta_{name}_lo"], rec[f"delta_{name}_hi"] = q(dl)
        rows.append(rec)
    pd.DataFrame(rows).to_csv(RES / f"e13_model_layer_trends{SFX}.csv", index=False)

    pooled = []
    for set_name, rule in SETS.items():
        for bsel in ["all"] + sorted({d["basin"] for d in grp}):
            idx = [i for i, d in enumerate(grp) if rule(upper(d["layer"])) and bsel in ("all", d["basin"])]
            if not idx:
                continue
            n = np.array([grp[i]["n"] for i in idx], float)
            for name, (y, x) in PAIRS.items():
                def pool(e):
                    by, bx = e[idx, j[y]], e[idx, j[x]]
                    ok = np.isfinite(by) & np.isfinite(bx)
                    return np.sum(n[ok] * (by - bx)[ok]) / n[ok].sum(), np.sum(n[ok] * by[ok]) / np.sum(n[ok] * bx[ok])
                p, bs = pool(point), np.array([pool(b) for b in boots])
                pooled.append(dict(layer_set=set_name, basin=bsel, test=name, n_groups=len(idx), n=int(n.sum()),
                                   delta=p[0], delta_lo=q(bs[:, 0])[0], delta_hi=q(bs[:, 0])[1],
                                   R=p[1], R_lo=q(bs[:, 1])[0], R_hi=q(bs[:, 1])[1]))
                if bsel == "all":
                    print(pooled[-1], flush=True)
    pd.DataFrame(pooled).to_csv(RES / f"e13_model_pooled{SFX}.csv", index=False)


def differential(df, cfg, n_boot):
    frames = []
    for basin, layer, st, mask in basin_frames(cfg, df=df):
        frames.append(dict(basin=basin, layer=layer, X=design(st, mask, ("theta", "salinity", "aou")),
                           Y=st[RESP].to_numpy(float), a=_codes(st["cell"].astype(str) + "|" + np.where(mask, "T", "R")),
                           cruise=st["cruise"].to_numpy()))
    cru = np.unique(np.concatenate([f["cruise"] for f in frames]))
    for f in frames:
        f["ci"] = np.searchsorted(cru, f["cruise"])

    def est(wc):
        res = []
        for f in frames:
            w = np.ones(len(f["ci"])) if wc is None else wc[f["ci"]]
            k = w > 0
            aa = np.unique(f["a"][k], return_inverse=True)[1]
            bb = np.unique(f["ci"][k], return_inverse=True)[1]
            res.append(exact_twoway(f["X"][k], f["Y"][k], aa, bb, w[k])[-1])
        return np.array(res)

    point, boots = est(None), [est(w) for w in boot_weights(len(cru), n_boot)]
    B = np.array(boots)
    j = {r: i for i, r in enumerate(RESP)}
    rows = []
    for i, f in enumerate(frames):
        rec = dict(basin=f["basin"], layer=f["layer"], **{f"dbeta_{r}": point[i, j[r]] for r in RESP})
        for name, (y, x) in PAIRS.items():
            rec[f"ddelta_{name}"] = point[i, j[y]] - point[i, j[x]]
            rec[f"ddelta_{name}_lo"], rec[f"ddelta_{name}_hi"] = q(B[:, i, j[y]] - B[:, i, j[x]])
        rows.append(rec)
    pd.DataFrame(rows).to_csv(RES / f"e13_model_differential{SFX}.csv", index=False)

    sets = {b: [i for i, f in enumerate(frames) if f["basin"] == b] for b in sorted({f["basin"] for f in frames})}
    sets["all"] = list(range(len(frames)))
    out = []
    for name, (y, x) in PAIRS.items():
        def scale(r, idx):
            dd, dc = r[idx, j[y]] - r[idx, j[x]], r[idx, j[x]]
            return float(np.sum(dd * dc) / np.sum(dc * dc))
        for s, idx in sets.items():
            bs = np.array([scale(b, idx) for b in boots])
            out.append(dict(test=name, set=s, n_layers=len(idx), m=scale(point, idx), m_lo=q(bs)[0], m_hi=q(bs)[1]))
            print(out[-1], flush=True)
    pd.DataFrame(out).to_csv(RES / f"e13_model_scale_factor{SFX}.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=500)
    ap.add_argument("--trace-variant", default="default", choices=["default", "atm"])  # atm: E13d, post hoc
    ap.add_argument("--model", default="mri_esm2_0", choices=["mri_esm2_0", "noresm2_lm"])
    ap.add_argument("--max-year", type=int, default=None)
    args = ap.parse_args()
    global SFX
    SFX = ("" if args.model == "mri_esm2_0" else f"_{args.model}") + \
          ("" if args.max_year is None else f"_to{args.max_year}") + \
          ("" if args.trace_variant == "default" else f"_{args.trace_variant}")
    cfg = load_config()
    df = model_frame(cfg, args.trace_variant, args.model, args.max_year)
    single_layer(df, cfg, args.boot)
    differential(df, cfg, args.boot)


if __name__ == "__main__":
    main()
