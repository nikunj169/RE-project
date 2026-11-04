"""E02a: evaluate TRACE anthropogenic carbon at every strict-QC sample (requires tracepy).

Outputs research/data/trace_canth.csv.gz with columns: row_id, canth, u_canth,
canth_plus, canth_minus (Cant at t +/- 0.5 yr with identical inputs, for the direct
dCant/dt like-for-like check).

Usage: python -m research.experiments.e02a_run_trace [--chunk 20000]
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from research.data_io import ROOT, load_config, load_data


def run_trace(lon, lat, depth, sal, temp, year, cfg):
    from tracepy import trace  # imported lazily: optional dependency

    out = trace(
        output_coordinates=np.column_stack([lon, lat, depth]),
        dates=np.asarray(year, float),
        predictor_measurements=np.column_stack([sal, temp]),
        predictor_types=np.array([1, 2]),
        atm_co2_trajectory=cfg["trace"]["atm_co2_trajectory"],
        canth_diseq=cfg["trace"]["canth_diseq"],
        eos=cfg["trace"]["eos"],
        output_format="pandas",
        verbose_tf=False,
    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", type=int, default=20000)
    args = ap.parse_args()
    cfg = load_config()
    df = load_data(cfg).reset_index().rename(columns={"index": "row_id"})
    lon = np.where(df.longitude > 180, df.longitude - 360, df.longitude)
    parts = []
    t0 = time.time()
    for s in range(0, len(df), args.chunk):
        sl = slice(s, s + args.chunk)
        g = df.iloc[sl]
        res = {}
        for tag, dt in (("canth", 0.0), ("canth_plus", 0.5), ("canth_minus", -0.5)):
            o = run_trace(lon[sl], g.latitude, g.depth, g.salinity, g.temperature,
                          g.decimal_year + dt, cfg)
            res[tag] = np.asarray(o["canth"], float)
            if tag == "canth":
                res["u_canth"] = np.asarray(o["u_canth"], float)
        parts.append(pd.DataFrame({"row_id": g.row_id.to_numpy(), **res}))
        print(f"{min(s + args.chunk, len(df))}/{len(df)} rows, {time.time() - t0:.0f}s", flush=True)
    out = pd.concat(parts, ignore_index=True)
    out.to_csv(ROOT / "research" / "data" / "trace_canth.csv.gz", index=False)
    print(out.describe().T)


if __name__ == "__main__":
    main()
