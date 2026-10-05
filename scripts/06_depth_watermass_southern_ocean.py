"""
Step 06: Depth, water-mass, and Southern Ocean failure analyses.

Runs all the detailed analyses that characterize where the models
succeed and fail.
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
from src.analysis.depth_analysis import depth_stratified_comparison
from src.analysis.watermass import classify_water_mass, watermass_performance
from src.analysis.derivatives import derivative_summary
from src.analysis.southern_ocean import analyze_abyssal_failure

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def main():
    logger.info("=" * 70)
    logger.info("STEP 06: Depth, Water-Mass, and Southern Ocean Analysis")
    logger.info("=" * 70)

    config = load_config()
    output_dir = get_path("results_processed")
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    data_dir = get_path("processed_data")
    df = pd.read_csv(data_dir / "glodap_qc_filtered.csv")
    df["basin"] = assign_basin(df, config.get("basins"))
    splits = split_temporal(df,
                            train_end=config["temporal_split"]["train_end_year"],
                            val_end=config["temporal_split"]["validation_end_year"])

    # ── Depth-stratified comparison ─────────────────────────────────────
    logger.info("\nDepth-stratified analysis...")
    depth_results = []

    for basin_name in ["Atlantic", "Indian", "Pacific", "Southern Ocean"]:
        mask = (df["basin"] == basin_name).values
        sub = df[mask]
        tr = splits["train"][mask]
        te = splits["validation"][mask]

        S, T, A, Y = sub["salinity"].values, sub["temperature"].values, \
                      sub["aou"].values, sub["tco2"].values
        depth = sub["depth"].values

        if te.sum() < 100:
            continue

        # Fit rank-1 and full quadratic
        r1 = Rank1Quadratic()
        r1.fit(S[tr], T[tr], A[tr], Y[tr])

        fq = FullQuadratic()
        fq.fit(S[tr], T[tr], A[tr], Y[tr])

        Yp_r1 = r1.predict(S[te], T[te], A[te])
        Yp_fq = fq.predict(S[te], T[te], A[te])

        dc = depth_stratified_comparison(
            Y[te], Yp_r1, Yp_fq, depth[te],
            label_a="Rank1", label_b="FullQuad")
        dc["Basin"] = basin_name
        depth_results.append(dc)

    if depth_results:
        df_depth = pd.concat(depth_results, ignore_index=True)
        df_depth.to_csv(output_dir / "depth_stratified_comparison.csv", index=False)
        logger.info("  Saved depth_stratified_comparison.csv")

    # ── Water-mass analysis (Atlantic) ──────────────────────────────────
    logger.info("\nWater-mass analysis (Atlantic)...")
    atl_mask = (df["basin"] == "Atlantic").values
    atl = df[atl_mask]
    tr_atl = splits["train"][atl_mask]
    te_atl = splits["validation"][atl_mask]

    if te_atl.sum() > 100:
        S, T, A, Y = atl["salinity"].values, atl["temperature"].values, \
                      atl["aou"].values, atl["tco2"].values

        r1 = Rank1Quadratic()
        r1.fit(S[tr_atl], T[tr_atl], A[tr_atl], Y[tr_atl])
        Yp_r1 = r1.predict(S[te_atl], T[te_atl], A[te_atl])

        fq = FullQuadratic()
        fq.fit(S[tr_atl], T[tr_atl], A[tr_atl], Y[tr_atl])
        Yp_fq = fq.predict(S[te_atl], T[te_atl], A[te_atl])

        wm_r1 = watermass_performance(Y[te_atl], Yp_r1, atl[te_atl])
        wm_r1["Model"] = "Rank1"
        wm_r1["Basin"] = "Atlantic"

        wm_fq = watermass_performance(Y[te_atl], Yp_fq, atl[te_atl])
        wm_fq["Model"] = "FullQuad"
        wm_fq["Basin"] = "Atlantic"

        wm = pd.concat([wm_r1, wm_fq], ignore_index=True)
        wm.to_csv(output_dir / "watermass_atlantic.csv", index=False)
        logger.info("  Saved watermass_atlantic.csv")

    # ── Derivative analysis ─────────────────────────────────────────────
    logger.info("\nDerivative analysis...")
    deriv_results = []

    for basin_name in ["Atlantic", "Indian", "Pacific", "Southern Ocean"]:
        mask = (df["basin"] == basin_name).values
        sub = df[mask]
        tr = splits["train"][mask]

        S, T, A, Y = sub["salinity"].values, sub["temperature"].values, \
                      sub["aou"].values, sub["tco2"].values

        r1 = Rank1Quadratic()
        r1.fit(S[tr], T[tr], A[tr], Y[tr])
        p = r1.params_

        ds = derivative_summary(
            S, T, A, p["alpha"], p["beta"], p["gamma"], p["delta"])
        ds["basin"] = basin_name
        deriv_results.append(ds)

    df_deriv = pd.DataFrame(deriv_results)
    df_deriv.to_csv(output_dir / "derivative_analysis.csv", index=False)
    logger.info("  Saved derivative_analysis.csv")

    # ── Southern Ocean abyssal failure ──────────────────────────────────
    logger.info("\nSouthern Ocean abyssal failure analysis...")
    so_mask = (df["basin"] == "Southern Ocean").values
    so = df[so_mask]

    if len(so) > 100:
        result = analyze_abyssal_failure(
            so["salinity"].values, so["temperature"].values,
            so["aou"].values, so["tco2"].values,
            so["depth"].values, so["year"].values,
            train_end=config["temporal_split"]["train_end_year"],
            val_end=config["temporal_split"]["validation_end_year"])

        logger.info(f"  Global model: R²={result['global_r2']:.4f}, RMSE={result['global_rmse']:.2f}")
        logger.info(f"  Local model:  R²={result['local_r2']:.4f}, RMSE={result['local_rmse']:.2f}")
        logger.info(f"  Failure is coefficient-transfer: {result['failure_is_coefficient_transfer']}")

        # Save
        pd.DataFrame([result]).to_csv(output_dir / "southern_ocean_abyssal.csv", index=False)
        logger.info("  Saved southern_ocean_abyssal.csv")

    logger.info("\nStep 06 complete ✓")


if __name__ == "__main__":
    main()
