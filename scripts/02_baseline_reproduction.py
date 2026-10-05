"""
Step 02: Reproduce the old baseline results.

Fits the rank-1 quadratic and full quadratic models per basin,
computes all metrics, and saves a baseline reproduction report.

This is the critical verification step — if results don't match the old
paper, we log discrepancies rather than silently adjusting.
"""

import sys
import os
import logging
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.config_loader import load_config, get_path
from src.basins import assign_basin
from src.splits import split_temporal
from src.models.rank1_quadratic import Rank1Quadratic
from src.models.full_quadratic import FullQuadratic
from src.models.linear import LinearModel
from src.analysis.metrics import compute_metrics
from src.analysis.depth_analysis import depth_stratified_metrics

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

# Old paper reference values (for discrepancy logging)
OLD_PAPER_VALUES = {
    "Atlantic": {"R2_test": 0.907, "RMSE_test": 20.17, "RMSE_poly": 23.41,
                  "alpha": 1.4859, "beta": -0.2377, "gamma": 2.2196,
                  "delta": -36.62, "epsilon": 1923.5},
    "Indian": {"R2_test": 0.985, "RMSE_test": 14.49, "RMSE_poly": 13.73,
                "alpha": 0.4992, "beta": -0.1693, "gamma": 1.2860,
                "delta": 10.07, "epsilon": 1433.1},
    "Pacific": {"R2_test": 0.985, "RMSE_test": 16.00, "RMSE_poly": 15.50,
                 "alpha": 1.1594, "beta": -0.2696, "gamma": 1.7753,
                 "delta": -20.15, "epsilon": 1790.3},
    "Southern Ocean": {"R2_test": 0.968, "RMSE_test": 10.58, "RMSE_poly": 10.46,
                        "alpha": 0.7289, "beta": -0.2307, "gamma": 1.6378,
                        "delta": -6.20, "epsilon": 1811.2},
}


def main():
    logger.info("=" * 70)
    logger.info("STEP 02: Baseline Reproduction")
    logger.info("=" * 70)

    config = load_config()
    output_dir = get_path("results_processed")
    os.makedirs(output_dir, exist_ok=True)

    # Load processed data
    data_dir = get_path("processed_data")
    df = pd.read_csv(data_dir / "glodap_qc_filtered.csv")
    logger.info(f"Loaded {len(df):,} observations")

    # Assign basins and splits
    df["basin"] = assign_basin(df, config.get("basins"))
    splits = split_temporal(df,
                            train_end=config["temporal_split"]["train_end_year"],
                            val_end=config["temporal_split"]["validation_end_year"])

    all_results = []
    depth_results = []
    discrepancies = []
    coefficients = []

    for basin_name in ["Atlantic", "Indian", "Pacific", "Southern Ocean"]:
        logger.info(f"\n{'─' * 60}")
        logger.info(f"BASIN: {basin_name}")

        mask = (df["basin"] == basin_name).values
        sub = df[mask]

        tr = splits["train"][mask]
        te = splits["validation"][mask]
        ex = splits["external"][mask]

        S = sub["salinity"].values
        T = sub["temperature"].values
        A = sub["aou"].values
        Y = sub["tco2"].values
        depth = sub["depth"].values

        n_train = tr.sum()
        n_val = te.sum()
        n_ext = ex.sum()
        logger.info(f"  n_train={n_train:,}  n_val={n_val:,}  n_ext={n_ext:,}")

        if n_val < 100:
            logger.warning(f"  Skipping {basin_name}: too few validation samples")
            continue

        # ── Rank-1 Quadratic ────────────────────────────────────────────
        logger.info("  Fitting rank-1 quadratic...")
        r1 = Rank1Quadratic()
        r1.fit(S[tr], T[tr], A[tr], Y[tr])

        Yp_val_r1 = r1.predict(S[te], T[te], A[te])
        m_val_r1 = compute_metrics(Y[te], Yp_val_r1)

        Yp_ext_r1 = r1.predict(S[ex], T[ex], A[ex])
        m_ext_r1 = compute_metrics(Y[ex], Yp_ext_r1)

        logger.info(f"  R1 Val: R²={m_val_r1['R2']:.4f}  RMSE={m_val_r1['RMSE']:.2f}")
        logger.info(f"  R1 Ext: R²={m_ext_r1['R2']:.4f}  RMSE={m_ext_r1['RMSE']:.2f}")

        # ── Full Quadratic Benchmark ────────────────────────────────────
        logger.info("  Fitting full quadratic benchmark...")
        fq = FullQuadratic()
        fq.fit(S[tr], T[tr], A[tr], Y[tr])

        Yp_val_fq = fq.predict(S[te], T[te], A[te])
        m_val_fq = compute_metrics(Y[te], Yp_val_fq)

        improvement = (m_val_fq["RMSE"] - m_val_r1["RMSE"]) / m_val_fq["RMSE"] * 100
        logger.info(f"  FQ Val: R²={m_val_fq['R2']:.4f}  RMSE={m_val_fq['RMSE']:.2f}")
        logger.info(f"  Improvement: {improvement:+.1f}%")

        # ── Coefficient comparison ──────────────────────────────────────
        p = r1.params_
        coefficients.append({
            "Basin": basin_name,
            "alpha": p["alpha"], "beta": p["beta"], "gamma": p["gamma"],
            "delta": p["delta"], "epsilon": p["epsilon"],
        })

        # ── Check against old paper ────────────────────────────────────
        if basin_name in OLD_PAPER_VALUES:
            old = OLD_PAPER_VALUES[basin_name]
            rmse_diff = abs(m_val_r1["RMSE"] - old["RMSE_test"])
            if rmse_diff > 1.0:  # tolerance of 1 μmol/kg
                discrepancy = {
                    "Basin": basin_name,
                    "Metric": "RMSE_test",
                    "Old_value": old["RMSE_test"],
                    "Reproduced": round(m_val_r1["RMSE_test"], 2) if "RMSE_test" in m_val_r1 else round(m_val_r1["RMSE"], 2),
                    "Difference": rmse_diff,
                    "Likely_reason": "Different QC, split implementation, or fitting precision",
                    "Resolution": "Under investigation",
                }
                discrepancies.append(discrepancy)
                logger.warning(f"  DISCREPANCY: RMSE old={old['RMSE_test']:.2f} vs new={m_val_r1['RMSE']:.2f}")

        # ── Store results ──────────────────────────────────────────────
        all_results.append({
            "Basin": basin_name,
            "N_train": int(n_train), "N_val": int(n_val), "N_ext": int(n_ext),
            "R2_val_r1": round(m_val_r1["R2"], 4),
            "RMSE_val_r1": round(m_val_r1["RMSE"], 3),
            "MAE_val_r1": round(m_val_r1["MAE"], 3),
            "Bias_val_r1": round(m_val_r1["Bias"], 3),
            "R2_ext_r1": round(m_ext_r1["R2"], 4),
            "RMSE_ext_r1": round(m_ext_r1["RMSE"], 3),
            "R2_val_fq": round(m_val_fq["R2"], 4),
            "RMSE_val_fq": round(m_val_fq["RMSE"], 3),
            "MAE_val_fq": round(m_val_fq["MAE"], 3),
            "Improvement_pct": round(improvement, 2),
        })

        # ── Depth-stratified ───────────────────────────────────────────
        logger.info("  Computing depth-stratified metrics...")
        dm = depth_stratified_metrics(Y[te], Yp_val_r1, depth[te])
        dm["Basin"] = basin_name

        # Also compute for full quadratic
        dm_fq = depth_stratified_metrics(Y[te], Yp_val_fq, depth[te])
        dm_fq = dm_fq.rename(columns={"RMSE": "RMSE_fq", "R2": "R2_fq", "MAE": "MAE_fq"})

        dm_merged = dm.merge(dm_fq[["depth_bin", "RMSE_fq", "R2_fq"]],
                              on="depth_bin", how="left")
        depth_results.append(dm_merged)

    # ── Save results ────────────────────────────────────────────────────
    logger.info("\n" + "=" * 70)
    logger.info("Saving baseline reproduction results...")

    df_results = pd.DataFrame(all_results)
    df_results.to_csv(output_dir / "baseline_reproduction.csv", index=False)
    logger.info("  Saved baseline_reproduction.csv")

    if depth_results:
        df_depth = pd.concat(depth_results, ignore_index=True)
        df_depth.to_csv(output_dir / "baseline_depth_stratified.csv", index=False)
        logger.info("  Saved baseline_depth_stratified.csv")

    df_coeffs = pd.DataFrame(coefficients)
    df_coeffs.to_csv(output_dir / "baseline_coefficients.csv", index=False)
    logger.info("  Saved baseline_coefficients.csv")

    # ── Discrepancy log ────────────────────────────────────────────────
    if discrepancies:
        df_disc = pd.DataFrame(discrepancies)
        df_disc.to_csv(output_dir / "discrepancy_log.csv", index=False)
        logger.warning(f"  Saved {len(discrepancies)} discrepancies to discrepancy_log.csv")
    else:
        logger.info("  No significant discrepancies found ✓")

    # ── Summary report ─────────────────────────────────────────────────
    report_path = output_dir / "baseline_reproduction_report.md"
    with open(report_path, "w") as f:
        f.write("# Baseline Reproduction Report\n\n")
        f.write(f"Generated from {len(df):,} QC-passing GLODAP observations.\n\n")
        f.write("## Per-Basin Validation Metrics\n\n")
        f.write(df_results.to_markdown(index=False))
        f.write("\n\n## Coefficients\n\n")
        f.write(df_coeffs.to_markdown(index=False))
        if discrepancies:
            f.write("\n\n## Discrepancies\n\n")
            f.write(pd.DataFrame(discrepancies).to_markdown(index=False))
        else:
            f.write("\n\n## Discrepancies\n\nNo significant discrepancies found.\n")

    logger.info(f"  Saved baseline_reproduction_report.md")
    logger.info("\nStep 02 complete ✓")


if __name__ == "__main__":
    main()
