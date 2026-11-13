"""E07: TRACE structural sensitivity of the DIC-vs-Cant comparison (pre-registered, GAP_AUDIT_ROUND2).

TRACE variants, each evaluated at every strict-QC sample:
  dg1.0 / dg1.6   : TTD shape Delta/Gamma (default 1.30); age scale factors re-estimated by TRACE
  diseq0.9 / 1.1  : scale on the atmospheric CO2 anomaly (air-sea disequilibrium)
  sf0.8 / sf1.2   : default mean-age scale factors x 0.8 / x 1.2 (younger / older water)
  theta           : potential instead of in-situ temperature as input
Then the E02 estimator is refitted jointly for DIC and each Cant variant (shared bootstrap).

Usage: python -m research.experiments.e07_trace_sensitivity [--boot 300] [--skip-trace]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.models.within_cell_trend import fit_trend

VARIANTS = {
    "dg1.0": {"delta_over_gamma": 1.0},
    "dg1.6": {"delta_over_gamma": 1.6},
    "diseq0.9": {"canth_diseq": 0.9},
    "diseq1.1": {"canth_diseq": 1.1},
    "sf0.8": {"sf_mult": 0.8},
    "sf1.2": {"sf_mult": 1.2},
    "theta": {"temp_col": "theta"},
}
OUT = ROOT / "research" / "data" / "trace_variants.csv.gz"


def trace_call(g: pd.DataFrame, lon: np.ndarray, cfg: dict, temp_col="temperature", **kw) -> pd.DataFrame:
    from tracepy import trace

    return trace(
        output_coordinates=np.column_stack([lon, g.latitude, g.depth]),
        dates=g.decimal_year.to_numpy(float),
        predictor_measurements=np.column_stack([g.salinity, g[temp_col]]),
        predictor_types=np.array([1, 2]),
        atm_co2_trajectory=cfg["trace"]["atm_co2_trajectory"],
        canth_diseq=kw.get("canth_diseq", cfg["trace"]["canth_diseq"]),
        eos=cfg["trace"]["eos"],
        delta_over_gamma=kw.get("delta_over_gamma"),
        scale_factors=kw.get("scale_factors"),
        output_format="pandas",
        verbose_tf=False,
    )


def run_variants(cfg: dict, chunk: int = 25000) -> None:
    df = load_data(cfg).reset_index().rename(columns={"index": "row_id"})
    lon = np.where(df.longitude > 180, df.longitude - 360, df.longitude)
    parts, t0 = [], time.time()
    for s in range(0, len(df), chunk):
        sl = slice(s, s + chunk)
        g = df.iloc[sl]
        base = trace_call(g, lon[sl], cfg)
        res = {"row_id": g.row_id.to_numpy(), "canth_default": base.canth.to_numpy(float),
               "sf_default": base.scale_factors.to_numpy(float), "mean_age": base.mean_age.to_numpy(float)}
        for name, kw in VARIANTS.items():
            kw = dict(kw)
            if "sf_mult" in kw:
                sf = res["sf_default"] * kw.pop("sf_mult")
                ok = np.isfinite(sf)
                out = np.full(len(g), np.nan)
                if ok.any():
                    o = trace_call(g[ok], lon[sl][ok], cfg, scale_factors=sf[ok], **kw)
                    out[ok] = o.canth.to_numpy(float)
                res[f"canth_{name}"] = out
            else:
                res[f"canth_{name}"] = trace_call(g, lon[sl], cfg, **kw).canth.to_numpy(float)
        parts.append(pd.DataFrame(res))
        print(f"{min(s + chunk, len(df))}/{len(df)} {time.time() - t0:.0f}s", flush=True)
    pd.concat(parts).to_csv(OUT, index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=300)
    ap.add_argument("--skip-trace", action="store_true")
    args = ap.parse_args()
    cfg = load_config()
    if not args.skip_trace or not OUT.exists():
        run_variants(cfg)
    tv = pd.read_csv(OUT).set_index("row_id")
    df = load_data(cfg).join(tv, how="inner")
    cols = ["canth_default"] + [f"canth_{v}" for v in VARIANTS]
    df = df[np.isfinite(df[cols]).all(axis=1)]
    rows = []
    for basin, layer, g in groups(df, cfg):
        fit = fit_trend(g, ["tco2"] + cols, n_boot=args.boot, seed=cfg["bootstrap"]["seed"])
        for c in cols:
            d = fit.boot["tco2"] - fit.boot[c]
            r = fit.boot["tco2"] / fit.boot[c]
            rows.append(dict(basin=basin, layer=layer, variant=c.removeprefix("canth_"), n=fit.n,
                             beta_dic=fit.beta["tco2"], beta_cant=fit.beta[c],
                             delta=fit.beta["tco2"] - fit.beta[c],
                             delta_lo=np.nanquantile(d, 0.025), delta_hi=np.nanquantile(d, 0.975),
                             R=fit.beta["tco2"] / fit.beta[c],
                             R_lo=np.nanquantile(r, 0.025), R_hi=np.nanquantile(r, 0.975),
                             mean_age_median=float(g.mean_age.median())))
        print(basin, layer, {c.removeprefix('canth_'): round(fit.beta['tco2'] - fit.beta[c], 3) for c in cols}, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e07_trace_sensitivity.csv", index=False)


if __name__ == "__main__":
    main()
