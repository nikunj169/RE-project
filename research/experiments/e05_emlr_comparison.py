"""E05: compare observed within-cell DIC trends and TRACE Cant rates with the independent
eMLR(C*) Cant change of Gruber et al. (2019), 1994-2007, per basin x density layer.

Design and decision rule are fixed in research_docs/EXPERIMENT_LOG.md (E05).

Usage: python -m research.experiments.e05_emlr_comparison [--boot 300]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import xarray as xr

from research.data_io import ROOT, groups, load_config, load_data
from research.experiments.e02a_run_trace import run_trace
from research.models.within_cell_trend import fit_trend

NC = ROOT / "data" / "external" / "gruber2019" / "dcant_emlr_cstar_gruber_94-07_vs1.nc"
YEARS = 13.0


def emlr_at(ds: xr.Dataset, var: str, lat, lon, depth) -> np.ndarray:
    lon360 = np.mod(lon, 360)
    da = ds[var]
    v = da.interp(DEPTH=xr.DataArray(depth, dims="p"),
                  LATITUDE=xr.DataArray(lat, dims="p"),
                  LONGITUDE=xr.DataArray(lon360, dims="p"), method="nearest")
    # depth: linear between levels at the nearest horizontal point
    vlin = da.sel(LATITUDE=xr.DataArray(lat, dims="p"), LONGITUDE=xr.DataArray(lon360, dims="p"),
                  method="nearest").interp(DEPTH=xr.DataArray(depth, dims="p"))
    out = np.asarray(vlin, float)
    fill = ~np.isfinite(out)
    out[fill] = np.asarray(v, float)[fill]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=300)
    args = ap.parse_args()
    cfg = load_config()
    ds = xr.open_dataset(NC)
    df = load_data(cfg)
    tr = pd.read_csv(ROOT / "research" / "data" / "trace_canth.csv.gz").set_index("row_id")
    df = df.join(tr[["canth"]], how="inner").dropna(subset=["canth"])
    df = df[(df.year >= 1989) & (df.year <= 2012)]
    lon = np.where(df.longitude > 180, df.longitude - 360, df.longitude)
    df = df.assign(lon180=lon)

    rows = []
    for basin, layer, g in groups(df, cfg):
        fit = fit_trend(g, ["tco2", "canth"], n_boot=args.boot, seed=cfg["bootstrap"]["seed"],
                        ratio=("tco2", "canth"))
        lat, lo, dep = g.latitude.to_numpy(), g.lon180.to_numpy(), g.depth.to_numpy()
        ens = np.array([np.nanmean(emlr_at(ds, f"DCANT_{k:02d}", lat, lo, dep)) / YEARS for k in range(1, 15)])
        e01 = ens[0]
        n_valid = int(np.isfinite(emlr_at(ds, "DCANT_01", lat, lo, dep)).sum())
        sig = np.nanmean(emlr_at(ds, "SIGMA_DCANT_01", lat, lo, dep)) / YEARS
        a = run_trace(lo, lat, dep, g.salinity, g.temperature, np.full(len(g), 1994.0), cfg)
        b = run_trace(lo, lat, dep, g.salinity, g.temperature, np.full(len(g), 2007.0), cfg)
        trace_rate = float(np.nanmean((np.asarray(b["canth"], float) - np.asarray(a["canth"], float)) / YEARS))
        bd, bc = fit.beta["tco2"], fit.beta["canth"]
        delta = bd - bc
        if delta > 0:
            verdict = ("TRACE_low" if e01 >= bd - 0.5 * delta else
                       "nonCant_or_drift" if e01 <= trace_rate + 0.5 * delta else "inconclusive")
        else:
            verdict = "n/a (delta<=0)"
        rows.append(dict(basin=basin, layer=layer, n=fit.n, n_emlr_valid=n_valid,
                         beta_dic=bd, beta_dic_lo=fit.ci("tco2")[0], beta_dic_hi=fit.ci("tco2")[1],
                         beta_cant_trace=bc, delta=delta, delta_lo=fit.ci("tco2-canth")[0],
                         delta_hi=fit.ci("tco2-canth")[1],
                         trace_rate_94_07=trace_rate, emlr_rate=e01, emlr_sigma=sig,
                         emlr_ens_min=float(ens.min()), emlr_ens_max=float(ens.max()), verdict=verdict))
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "research" / "results" / "e05_emlr_comparison.csv", index=False)
    print(out.verdict.value_counts())


if __name__ == "__main__":
    main()
