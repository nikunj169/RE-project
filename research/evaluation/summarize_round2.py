"""Round-2 summaries: E06 pre-registered verdicts, pooled eMLR/TRACE ratios (E10), C* robustness (E08).

Usage: python -m research.evaluation.summarize_round2   -> research/results/summary_round2_*.csv
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from research.data_io import ROOT

RES = ROOT / "research" / "results"


def upper(layer: str) -> float:
    return float(layer.strip("(]").split(",")[1])


def main() -> None:
    pd.set_option("display.width", 250)
    e02 = pd.read_csv(RES / "e02_layer_trends_primary.csv").set_index(["basin", "layer"])
    e06 = pd.read_csv(RES / "e06_differential_trends.csv").set_index(["basin", "layer"])
    # E06 rule (GAP_AUDIT_ROUND2 §2): excess not explained by cruise offsets if Ddelta CI > 0;
    # attributed to offsets/deep signal if CI includes 0 and |Ddelta| < 0.5 * |delta_E02|.
    v = e06[["ddelta", "ddelta_lo", "ddelta_hi"]].join(e02[["delta", "delta_ci_lo", "delta_ci_hi"]])
    v["verdict"] = np.select(
        [v.ddelta_lo > 0, (v.ddelta_lo <= 0) & (v.ddelta.abs() < 0.5 * v.delta.abs())],
        ["excess not explained by cruise offsets", "attributable to offset/deep signal"], "inconclusive")
    v.reset_index().to_csv(RES / "summary_round2_e06_verdicts.csv", index=False)
    print(v.round(3).to_string(), "\n", v.groupby(v.reset_index().layer.map(upper).le(27.0).values).verdict.value_counts())

    e10 = pd.read_csv(RES / "e10_emlr_products.csv")
    rows = []
    for name, rule in {"upper (<=27.0)": lambda u: u <= 27.0, "intermediate (27.0-27.8]": lambda u: 27.0 < u <= 27.8,
                       "deep (>27.8)": lambda u: u > 27.8}.items():
        for b in ["all"] + sorted(e10.basin.unique()):
            s = e10[e10.layer.map(upper).map(rule) & ((e10.basin == b) | (b == "all"))]
            w = s.n.to_numpy(float)
            rec = dict(layer_set=name, basin=b, n_groups=len(s))
            for key, src in (("g19_94_07", "g19_94_07"), ("m23_94_04", "m23_94_04"), ("m23_04_14", "m23_04_14")):
                rec[f"ratio_{key}"] = np.sum(w * s[f"emlr_{src}"]) / np.sum(w * s[f"trace_{key}"])
            rec["ratio_m23_94_04_ensmin"] = np.sum(w * s.emlr_m23_94_04_min) / np.sum(w * s.trace_m23_94_04)
            rec["ratio_m23_94_04_ensmax"] = np.sum(w * s.emlr_m23_94_04_max) / np.sum(w * s.trace_m23_94_04)
            rec["ratio_g19_ensmin"] = np.sum(w * s.emlr_g19_min) / np.sum(w * s.trace_g19_94_07)
            rec["ratio_g19_ensmax"] = np.sum(w * s.emlr_g19_max) / np.sum(w * s.trace_g19_94_07)
            rec["R_obs_pooled"] = np.sum(w * s.beta_dic) / np.sum(w * s.beta_cant)
            rows.append(rec)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "summary_round2_emlr_pooled.csv", index=False)
    print(out.round(3).to_string())

    e08 = pd.read_csv(RES / "e08_cstar_controls.csv")
    e08["cstar_ci_gt0"] = e08.delta_cstar_lo > 0
    e08["aou_ci_gt0"] = e08.delta_aou_lo > 0
    e08["talk_trend_ci_excl0"] = (e08.beta_talk_lo > 0) | (e08.beta_talk_hi < 0)
    e08.to_csv(RES / "summary_round2_cstar.csv", index=False)
    up = e08[e08.layer.map(upper) <= 27.0]
    print("\nC* upper layers: CI>0", int(up.cstar_ci_gt0.sum()), "/", len(up), "; AOU-version CI>0", int(up.aou_ci_gt0.sum()))
    print("TAlk trend CI excludes 0 in", int(e08.talk_trend_ci_excl0.sum()), "/", len(e08), "groups; sign:",
          np.sign(e08[e08.talk_trend_ci_excl0].beta_talk).value_counts().to_dict())


if __name__ == "__main__":
    main()
