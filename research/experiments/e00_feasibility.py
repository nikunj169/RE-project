"""E00 (exploratory, pre-registration): data coverage for a within-location trend design.

Reports, per basin x neutral-density layer, how many 4x4-degree cells have strict-QC
observations whose time span is >= 10 years, and a first-look within-cell DIC trend
conditional on (theta, S, AOU). Results are exploratory and are NOT used for inference.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GAMMA_EDGES = [0, 26.0, 26.5, 27.0, 27.5, 27.8, 28.0, 28.1, 28.2, 29.5]


def main() -> None:
    df = pd.read_csv(ROOT / "research" / "data" / "glodap_research.csv.gz")
    df = df[df.qc_strict & df.gamma.notna() & df.basin.isin(
        ["Atlantic", "Indian", "Pacific", "Southern Ocean"])].copy()
    df["layer"] = pd.cut(df.gamma, GAMMA_EDGES)
    lon = np.where(df.longitude > 180, df.longitude - 360, df.longitude)
    df["cell"] = (np.floor(df.latitude / 4).astype(int).astype(str) + "_"
                  + np.floor(lon / 4).astype(int).astype(str))
    rows = []
    for (b, lay), g in df.groupby(["basin", "layer"], observed=True):
        span = g.groupby("cell").decimal_year.agg(lambda s: s.max() - s.min())
        good = span[span >= 10].index
        gg = g[g.cell.isin(good)]
        rec = dict(basin=b, layer=str(lay), n=len(g), cells=g.cell.nunique(),
                   cells_span10=len(good), n_in_span10=len(gg), cruises=gg.cruise.nunique())
        if len(gg) > 500:
            # within-cell demeaning (cell fixed effects) of DIC, predictors and time
            cols = ["tco2", "theta", "salinity", "aou", "decimal_year"]
            dm = gg[cols] - gg.groupby("cell")[cols].transform("mean")
            X = np.column_stack([dm.theta, dm.salinity, dm.aou, dm.decimal_year])
            beta, *_ = np.linalg.lstsq(X, dm.tco2.to_numpy(), rcond=None)
            rec["dDIC_dt_within_cell"] = round(beta[-1], 3)
        rows.append(rec)
    out = pd.DataFrame(rows)
    print(out.to_string(index=False))
    out.to_csv(ROOT / "research" / "data" / "e00_feasibility.csv", index=False)


if __name__ == "__main__":
    main()
