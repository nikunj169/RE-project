"""E02: within-cell trends of observed DIC (and TRACE Cant if available) per basin x layer.

Primary analysis for H1/H2. If research/data/trace_canth.csv.gz exists, Cant is fitted with
the identical estimator on the same samples and joint bootstrap CIs are reported for
R = beta_DIC/beta_Cant and delta = beta_DIC - beta_Cant.

Usage:
  python -m research.experiments.e02_layer_trends [--variant primary] [--boot 500]
Variants implement the pre-registered sensitivity analyses S1-S10 (see methodology).
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data
from research.models.within_cell_trend import add_cells, fit_trend

VARIANTS = {
    "primary": {},
    "S1_relaxed_qc": {"strict": False},
    "S2_cell2": {"cell_deg": 2}, "S2_cell6": {"cell_deg": 6}, "S2_cell8": {"cell_deg": 8},
    "S3_span5": {"min_span": 5}, "S3_span15": {"min_span": 15},
    "S4_linear": {"form": "linear"},
    "S5_no_aou": {"variables": ("theta", "salinity")},
    "S6_no_theta_s": {"variables": ("aou",)},
    "S7_post1990": {"min_year": 1990},
    "S9_seasonal": {"seasonal": True},
    "S10_phosphate": {"variables": ("theta", "salinity", "phosphate")},
    "B1_no_conditioning": {"variables": ()},
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="primary", choices=list(VARIANTS))
    ap.add_argument("--boot", type=int, default=None)
    args = ap.parse_args()
    v = VARIANTS[args.variant]
    cfg = load_config()
    n_boot = args.boot if args.boot is not None else cfg["bootstrap"]["n"]

    df = load_data(cfg, strict=v.get("strict", True))
    if "cell_deg" in v:
        df["cell"] = add_cells(df, v["cell_deg"])
    if "min_year" in v:
        df = df[df.year >= v["min_year"]]
    variables = tuple(v.get("variables", cfg["conditioning"]["variables"]))
    df = df.dropna(subset=list(variables)) if variables else df
    extra = None
    if v.get("seasonal"):
        ang = 2 * np.pi * (df["month"].fillna(7) - 0.5) / 12
        df = df.assign(sin_m=np.sin(ang), cos_m=np.cos(ang))
        extra = ["sin_m", "cos_m"]

    trace_path = ROOT / "research" / "data" / "trace_canth.csv.gz"
    responses, ratio = ["tco2"], None
    if trace_path.exists() and v.get("strict", True):
        tr = pd.read_csv(trace_path).set_index("row_id")
        df = df.join(tr[["canth", "canth_plus", "canth_minus"]], how="inner")
        df = df[np.isfinite(df.canth)]
        responses, ratio = ["tco2", "canth"], ("tco2", "canth")

    rows = []
    t0 = time.time()
    for basin, layer, g in groups(df, cfg, v.get("min_span")):
        fit = fit_trend(g, responses, form=v.get("form", cfg["conditioning"]["form"]),
                        variables=variables, extra=extra, n_boot=n_boot,
                        seed=cfg["bootstrap"]["seed"], ratio=ratio)
        rec = dict(variant=args.variant, basin=basin, layer=layer, n=fit.n, n_cells=fit.n_cells,
                   n_cruises=fit.n_cruises, median_cell_span_yr=fit.median_span,
                   detection_limit=cfg["dic_offset_umol_kg"] / fit.median_span)
        for r in responses:
            rec[f"beta_{r}"] = fit.beta[r]
            rec[f"se_cr1_{r}"] = fit.se_cr1[r]
            if n_boot:
                rec[f"ci_lo_{r}"], rec[f"ci_hi_{r}"] = fit.ci(r)
        if ratio and n_boot:
            rec["R"] = fit.beta["tco2"] / fit.beta["canth"]
            rec["R_ci_lo"], rec["R_ci_hi"] = fit.ci("tco2/canth")
            rec["delta"] = fit.beta["tco2"] - fit.beta["canth"]
            rec["delta_ci_lo"], rec["delta_ci_hi"] = fit.ci("tco2-canth")
            rec["direct_dCant_dt"] = float((g.canth_plus - g.canth_minus).mean())
        rows.append(rec)
        print(f"{basin:15s} {layer:14s} n={fit.n:6d} beta_DIC={fit.beta['tco2']:+.3f} "
              f"({time.time() - t0:.0f}s)", flush=True)

    out = pd.DataFrame(rows)
    out_dir = ROOT / "research" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_dir / f"e02_layer_trends_{args.variant}.csv", index=False)
    meta = {"variant": args.variant, "settings": {k: str(x) for k, x in v.items()},
            "n_boot": n_boot, "seed": cfg["bootstrap"]["seed"], "responses": responses,
            "raw_sha256": cfg["data"]["raw_sha256"]}
    (out_dir / f"e02_layer_trends_{args.variant}.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
