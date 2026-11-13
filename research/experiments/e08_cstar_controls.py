"""E08: pump-corrected response (C*) and control variables (pre-registered, GAP_AUDIT_ROUND2).

C* = DIC - rCP*PO4 - 0.5*(TAlk + rNP*PO4), rCP = 117, rNP = 16 (Gruber et al., 1996), on samples
with strict-QC DIC and flag-2/secondary-QC-1 PO4 and TAlk. On these rows:
  fit A (conditioning theta, S):       C*, TRACE Cant, TAlk, AOU, PO4
  fit B (conditioning theta, S, AOU):  DIC, TRACE Cant   (Round-1 estimator on the same rows)
Both fits share the cruise bootstrap draws (same seed, same cruise set).

Usage: python -m research.experiments.e08_cstar_controls [--boot 300]
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.models.within_cell_trend import fit_trend

R_CP, R_NP = 117.0, 16.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=300)
    args = ap.parse_args()
    cfg = load_config()
    df = load_data(cfg)
    tr = pd.read_csv(ROOT / "research" / "data" / "trace_canth.csv.gz").set_index("row_id")
    df = df.join(tr[["canth"]], how="inner")
    ok = (np.isfinite(df.canth) & (df.phosphate_f == 2) & (df.phosphate_qc == 1)
          & (df.talk_f == 2) & (df.talk_qc == 1) & df.phosphate.notna() & df.talk.notna())
    df = df[ok].copy()
    df["cstar"] = df.tco2 - R_CP * df.phosphate - 0.5 * (df.talk + R_NP * df.phosphate)
    print(f"rows with C* inputs: {len(df)}")
    seed = cfg["bootstrap"]["seed"]
    rows = []
    q = lambda x: tuple(np.nanquantile(x, [0.025, 0.975]))
    for basin, layer, g in groups(df, cfg):
        a = fit_trend(g, ["cstar", "canth", "talk", "aou", "phosphate"], variables=("theta", "salinity"),
                      n_boot=args.boot, seed=seed)
        b = fit_trend(g, ["tco2", "canth"], n_boot=args.boot, seed=seed)
        dc = a.boot["cstar"] - a.boot["canth"]
        dd = b.boot["tco2"] - b.boot["canth"]
        rec = dict(basin=basin, layer=layer, n=a.n, n_cruises=a.n_cruises,
                   beta_cstar=a.beta["cstar"], beta_cant_ts=a.beta["canth"],
                   delta_cstar=a.beta["cstar"] - a.beta["canth"], delta_cstar_lo=q(dc)[0], delta_cstar_hi=q(dc)[1],
                   beta_dic_aou=b.beta["tco2"], beta_cant_tsa=b.beta["canth"],
                   delta_aou=b.beta["tco2"] - b.beta["canth"], delta_aou_lo=q(dd)[0], delta_aou_hi=q(dd)[1],
                   beta_talk=a.beta["talk"], beta_talk_lo=q(a.boot["talk"])[0], beta_talk_hi=q(a.boot["talk"])[1],
                   beta_aou=a.beta["aou"], beta_aou_lo=q(a.boot["aou"])[0], beta_aou_hi=q(a.boot["aou"])[1],
                   beta_po4=a.beta["phosphate"])
        rows.append(rec)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e08_cstar_controls.csv", index=False)


if __name__ == "__main__":
    main()
