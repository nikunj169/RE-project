"""E04 (S8): leave-one-cruise-out influence on beta_DIC, beta_Cant and delta per basin x layer.

Reports, per group, the largest absolute change in each estimate when any single cruise is
removed, the identity of that cruise, and the range of delta across all LOCO fits.

Usage: python -m research.experiments.e04_loco
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.models.within_cell_trend import fit_trend


def main() -> None:
    cfg = load_config()
    df = load_data(cfg)
    tr = pd.read_csv(ROOT / "research" / "data" / "trace_canth.csv.gz").set_index("row_id")
    df = df.join(tr[["canth"]], how="inner").dropna(subset=["canth"])
    expo = df.drop_duplicates("cruise").set_index("cruise")["expocode"].to_dict()
    rows = []
    for basin, layer, g in groups(df, cfg):
        full = fit_trend(g, ["tco2", "canth"])
        d0 = full.beta["tco2"] - full.beta["canth"]
        res = []
        for c in g.cruise.unique():
            h = g[g.cruise != c]
            if h.cell.nunique() < cfg["min_cells"]:
                continue
            f = fit_trend(h, ["tco2", "canth"])
            res.append((c, f.beta["tco2"], f.beta["canth"], f.beta["tco2"] - f.beta["canth"]))
        r = pd.DataFrame(res, columns=["cruise", "b_dic", "b_cant", "delta"])
        i = (r.b_dic - full.beta["tco2"]).abs().idxmax()
        rows.append(dict(basin=basin, layer=layer, n_cruises=len(r), beta_dic=full.beta["tco2"], delta=d0,
                         max_abs_change_beta_dic=float((r.b_dic - full.beta["tco2"]).abs().max()),
                         most_influential_cruise=expo.get(r.cruise[i], str(r.cruise[i])),
                         max_abs_change_delta=float((r.delta - d0).abs().max()),
                         delta_loco_min=float(r.delta.min()), delta_loco_max=float(r.delta.max())))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(ROOT / "research" / "results" / "e04_loco.csv", index=False)


if __name__ == "__main__":
    main()
