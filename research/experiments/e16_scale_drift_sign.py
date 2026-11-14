"""E16b (GAP_AUDIT_ROUND4 §1): can a uniform DIC scale drift explain both excess estimates?

Under DIC_obs = DIC_true * (1 + kappa * (t - 2000)), a layer's single-layer excess is kappa * DIC_layer
and the differential excess (target vs reference, cruise FE) is kappa * (DIC_target - DIC_ref).
For every E06 frame this reports the kappa each observed excess would require, using the
layer-mean observed DIC of the same rows. Opposite signs mean the drift cannot explain both.

Usage: python -m research.experiments.e16_scale_drift_sign
"""

from __future__ import annotations

import pandas as pd

from research.data_io import ROOT, load_config
from research.experiments.e06_differential_trends import basin_frames


def main() -> None:
    cfg = load_config()
    e06 = pd.read_csv(ROOT / "research" / "results" / "e06_differential_trends.csv").set_index(["basin", "layer"])
    rows = []
    for basin, layer, st, mask in basin_frames(cfg):
        r = e06.loc[(basin, layer)]
        d_t, d_r = st.tco2[mask].mean(), st.tco2[~mask].mean()
        rows.append(dict(basin=basin, layer=layer, dic_target=d_t, dic_ref=d_r,
                         delta_single=r.delta_target_same_rows, ddelta=r.ddelta,
                         kappa_single=r.delta_target_same_rows / d_t, kappa_diff=r.ddelta / (d_t - d_r)))
    out = pd.DataFrame(rows)
    out["opposite_sign"] = out.kappa_single * out.kappa_diff < 0
    out.to_csv(ROOT / "research" / "results" / "e16_scale_drift_sign.csv", index=False)
    print(out.round(6).to_string())
    up = out[out.layer.map(lambda l: float(l.strip("(]").split(",")[1]) <= 27.0)]
    print("upper layers with positive ddelta:", int((up.ddelta > 0).sum()), "/", len(up),
          "; of these with opposite-sign kappa:", int((up.opposite_sign & (up.ddelta > 0)).sum()))


if __name__ == "__main__":
    main()
