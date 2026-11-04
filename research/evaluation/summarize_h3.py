"""H3 summary: fraction of temporal-holdout bias removed by a time or Cant covariate, and
model-ranking stability across validation designs.

Criterion (pre-registered): >= 70 % of mean holdout bias removed, i.e.
    1 - |bias_with| / |bias_static| >= 0.70, per basin x model x split (mean model excluded).

Usage: python -m research.evaluation.summarize_h3 --tag strict_none-year-canth
"""

from __future__ import annotations

import argparse

import pandas as pd

from research.data_io import ROOT

RES = ROOT / "research" / "results"
MODELS = ["linear", "rank1", "rank2", "full_quadratic", "cubic"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="strict_none-year-canth")
    args = ap.parse_args()
    t = pd.read_csv(RES / f"e03_temporal_{args.tag}.csv")
    t = t[t.model.isin(MODELS) & t.split.isin(["validation", "external"])]
    w = t.pivot_table(index=["basin", "model", "split"], columns="covariate", values=["Bias", "RMSE"])
    rows = []
    for cov in [c for c in ("year", "canth") if c in t.covariate.unique()]:
        frac = 1 - w[("Bias", cov)].abs() / w[("Bias", "none")].abs()
        d_rmse = w[("RMSE", cov)] - w[("RMSE", "none")]
        part = pd.DataFrame({"covariate": cov, "bias_static": w[("Bias", "none")], "bias_with": w[("Bias", cov)],
                             "frac_bias_removed": frac, "rmse_static": w[("RMSE", "none")],
                             "rmse_with": w[("RMSE", cov)], "rmse_change": d_rmse}).reset_index()
        rows.append(part)
    out = pd.concat(rows)
    out["meets_70pct"] = out.frac_bias_removed >= 0.70
    out.to_csv(RES / f"summary_h3_{args.tag}.csv", index=False)

    # ranking per validation design (static models)
    b = pd.read_csv(RES / f"e03_blocked_summary_{args.tag}.csv")
    rk = []
    for cov in out.covariate.unique().tolist() + ["none"]:
        tt = t[(t.covariate == cov)]
        for (basin, split), g in tt.groupby(["basin", "split"]):
            rk.append(dict(covariate=cov, basin=basin, design=split,
                           ranking=" < ".join(g.sort_values("RMSE").model)))
        for (basin, blk), g in b[(b.covariate == cov) & b.model.isin(MODELS)].groupby(["basin", "block"]):
            rk.append(dict(covariate=cov, basin=basin, design=f"{blk}-block",
                           ranking=" < ".join(g.sort_values("RMSE_mean").model)))
    rk = pd.DataFrame(rk).drop_duplicates()
    rk.to_csv(RES / f"summary_h3_rankings_{args.tag}.csv", index=False)

    pd.set_option("display.width", 250)
    print(out.round(2).to_string(index=False))
    print(out.groupby(["covariate", "split"]).agg(n=("meets_70pct", "size"), n_meet=("meets_70pct", "sum"),
                                                   median_frac=("frac_bias_removed", "median"),
                                                   median_rmse_change=("rmse_change", "median")).round(3))
    print(rk[rk.covariate == "none"].to_string(index=False))


if __name__ == "__main__":
    main()
