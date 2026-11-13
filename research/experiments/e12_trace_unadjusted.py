"""E12 (post hoc, see GAP_AUDIT_ROUND2 §4): TRACE with the unadjusted atmospheric CO2 trajectory.

TRACE hard-codes data/CO2TrajectoriesAdjusted.txt (a damped surface-ocean boundary condition). Here
numpy.loadtxt is wrapped for the duration of each TRACE call so that file name resolves to the
unadjusted CO2Trajectories.txt; nothing in the installed package is modified.

Output: research/data/trace_canth_atm.csv.gz (row_id, canth). Refits are run with
  python -m research.experiments.e09_pooled_and_timeshape --cant atm
  python -m research.experiments.e11_scale_factor --cant atm
  python -m research.experiments.e12_trace_unadjusted --fit
"""

from __future__ import annotations

import argparse
import contextlib
import time

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.experiments.e02a_run_trace import run_trace
from research.models.within_cell_trend import fit_trend

OUT = ROOT / "research" / "data" / "trace_canth_atm.csv.gz"


@contextlib.contextmanager
def unadjusted_trajectory():
    orig = np.loadtxt

    def patched(fname, *a, **k):
        if str(fname).endswith("CO2TrajectoriesAdjusted.txt"):
            fname = str(fname).replace("CO2TrajectoriesAdjusted.txt", "CO2Trajectories.txt")
        return orig(fname, *a, **k)

    np.loadtxt = patched
    try:
        yield
    finally:
        np.loadtxt = orig


def run(cfg: dict, chunk: int = 25000) -> None:
    df = load_data(cfg).reset_index().rename(columns={"index": "row_id"})
    lon = np.where(df.longitude > 180, df.longitude - 360, df.longitude)
    parts, t0 = [], time.time()
    with unadjusted_trajectory():
        for s in range(0, len(df), chunk):
            sl = slice(s, s + chunk)
            g = df.iloc[sl]
            o = run_trace(lon[sl], g.latitude, g.depth, g.salinity, g.temperature, g.decimal_year, cfg)
            parts.append(pd.DataFrame({"row_id": g.row_id.to_numpy(), "canth": np.asarray(o["canth"], float)}))
            print(f"{min(s + chunk, len(df))}/{len(df)} {time.time() - t0:.0f}s", flush=True)
    pd.concat(parts).to_csv(OUT, index=False)


def fit(cfg: dict, n_boot: int) -> None:
    df = load_data(cfg)
    a = pd.read_csv(ROOT / "research" / "data" / "trace_canth.csv.gz").set_index("row_id")[["canth"]]
    b = pd.read_csv(OUT).set_index("row_id")[["canth"]].rename(columns={"canth": "canth_atm"})
    df = df.join(a, how="inner").join(b, how="inner")
    df = df[np.isfinite(df.canth) & np.isfinite(df.canth_atm)]
    rows = []
    for basin, layer, g in groups(df, cfg):
        f = fit_trend(g, ["tco2", "canth", "canth_atm"], n_boot=n_boot, seed=cfg["bootstrap"]["seed"])
        d = f.boot["tco2"] - f.boot["canth_atm"]
        rows.append(dict(basin=basin, layer=layer, n=f.n, beta_dic=f.beta["tco2"], beta_cant=f.beta["canth"],
                         beta_cant_atm=f.beta["canth_atm"], delta_atm=f.beta["tco2"] - f.beta["canth_atm"],
                         delta_atm_lo=np.nanquantile(d, 0.025), delta_atm_hi=np.nanquantile(d, 0.975),
                         R_atm=f.beta["tco2"] / f.beta["canth_atm"], rate_ratio_atm_default=f.beta["canth_atm"] / f.beta["canth"]))
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e12_trace_unadjusted.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", action="store_true")
    ap.add_argument("--boot", type=int, default=300)
    args = ap.parse_args()
    cfg = load_config()
    if not args.fit or not OUT.exists():
        run(cfg)
    if args.fit:
        fit(cfg, args.boot)


if __name__ == "__main__":
    main()
