"""Summaries for H1, H2 (robustness) and H4 from E02 outputs.

Writes research/results/summary_h1.csv, summary_robustness.csv, summary_h4.csv and prints them.
Classification rules (fixed here, applied mechanically):
  H1 verdict (layers whose upper gamma edge <= shallow_gamma_max):
    supported    : 95% CI of R overlaps [0.85, 1.15]
    contradicted : CI entirely outside the band
  Robust deviation: primary delta CI excludes 0 AND the sign of delta is identical in every
    sensitivity variant except B1 (baseline) and S1 (no TRACE).
  Above detection: |delta| >= detection limit (4 umol/kg / median within-cell span).
  H4: |beta_DIC(relaxed) - beta_DIC(strict)| > half-width of the primary CI, deep layers (gamma >= 27.8).
"""

from __future__ import annotations

import glob

import numpy as np
import pandas as pd

from research.data_io import ROOT, load_config


def upper_edge(layer: str) -> float:
    return float(layer.strip("(]").split(",")[1])


def main() -> None:
    cfg = load_config()
    res = ROOT / "research" / "results"
    prim = pd.read_csv(res / "e02_layer_trends_primary.csv")
    lo_b, hi_b = cfg["h1_band"]

    h1 = prim[prim.layer.map(upper_edge) <= cfg["shallow_gamma_max"]].copy()
    h1["verdict"] = np.where((h1.R_ci_hi < lo_b) | (h1.R_ci_lo > hi_b), "contradicted", "supported")
    h1 = h1[["basin", "layer", "beta_tco2", "beta_canth", "R", "R_ci_lo", "R_ci_hi", "verdict"]]
    h1.to_csv(res / "summary_h1.csv", index=False)

    var = pd.concat([pd.read_csv(f) for f in glob.glob(str(res / "e02_layer_trends_S*.csv"))])
    var = var[var.variant != "S1_relaxed_qc"]
    sign = var.groupby(["basin", "layer"]).delta.agg(lambda s: int(np.sign(s).nunique() == 1)).rename("sign_stable")
    smin = var.groupby(["basin", "layer"]).delta.min().rename("delta_min_variants")
    smax = var.groupby(["basin", "layer"]).delta.max().rename("delta_max_variants")
    rob = prim.set_index(["basin", "layer"])[["delta", "delta_ci_lo", "delta_ci_hi", "detection_limit", "R"]] \
        .join([sign, smin, smax]).reset_index()
    rob["ci_excludes_0"] = (rob.delta_ci_lo > 0) | (rob.delta_ci_hi < 0)
    rob["robust_deviation"] = rob.ci_excludes_0 & (rob.sign_stable == 1) & \
        (np.sign(rob.delta_min_variants) == np.sign(rob.delta))
    rob["above_detection_limit"] = rob.delta.abs() >= rob.detection_limit
    rob.to_csv(res / "summary_robustness.csv", index=False)

    rel = pd.read_csv(res / "e02_layer_trends_S1_relaxed_qc.csv").set_index(["basin", "layer"])
    h4 = prim.set_index(["basin", "layer"])[["beta_tco2", "ci_lo_tco2", "ci_hi_tco2"]] \
        .join(rel[["beta_tco2"]].rename(columns={"beta_tco2": "beta_relaxed"})).reset_index()
    h4 = h4[h4.layer.map(lambda s: float(s.strip("(]").split(",")[0])) >= 27.8]
    h4["abs_change"] = (h4.beta_relaxed - h4.beta_tco2).abs()
    h4["ci_half_width"] = (h4.ci_hi_tco2 - h4.ci_lo_tco2) / 2
    h4["qc_sensitive"] = h4.abs_change > h4.ci_half_width
    h4.to_csv(res / "summary_h4.csv", index=False)

    pd.set_option("display.width", 250)
    print("H1:\n", h1.round(3).to_string(index=False), "\n", h1.verdict.value_counts().to_dict())
    print("\nRobustness:\n", rob.round(3).to_string(index=False))
    print("\nH4:\n", h4.round(3).to_string(index=False), "\n qc_sensitive:", int(h4.qc_sensitive.sum()), "/", len(h4))


if __name__ == "__main__":
    main()
