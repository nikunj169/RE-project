"""
Step 05: Equivalence testing with dependence-aware methods.

For each basin:
  1. Run TOST at multiple margins (1, 2, 5, 10 μmol/kg)
  2. Run cluster bootstrap for paired AE differences
  3. Report both naive and dependence-aware results
"""

import sys
import os
import logging
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.config_loader import load_config, get_path
from src.basins import assign_basin
from src.splits import split_temporal, assign_spatial_blocks
from src.models.rank1_quadratic import Rank1Quadratic
from src.models.full_quadratic import FullQuadratic
from src.validation.equivalence import tost_multiple_margins, bootstrap_paired_ae_diff

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def main():
    logger.info("=" * 70)
    logger.info("STEP 05: Equivalence Analysis")
    logger.info("=" * 70)

    config = load_config()
    output_dir = get_path("results_processed")
    os.makedirs(output_dir, exist_ok=True)

    margins = config["equivalence"]["margins_umol_kg"]
    alpha = config["equivalence"]["alpha"]
    n_boot = config["bootstrap"]["n_replicates"]

    # Load data
    data_dir = get_path("processed_data")
    df = pd.read_csv(data_dir / "glodap_qc_filtered.csv")
    df["basin"] = assign_basin(df, config.get("basins"))
    splits = split_temporal(df,
                            train_end=config["temporal_split"]["train_end_year"],
                            val_end=config["temporal_split"]["validation_end_year"])

    all_results = []

    for basin_name in ["Atlantic", "Indian", "Pacific", "Southern Ocean"]:
        logger.info(f"\n{'─' * 60}")
        logger.info(f"BASIN: {basin_name}")

        mask = (df["basin"] == basin_name).values
        sub = df[mask]
        tr = splits["train"][mask]
        te = splits["validation"][mask]

        S = sub["salinity"].values
        T = sub["temperature"].values
        A = sub["aou"].values
        Y = sub["tco2"].values

        if te.sum() < 100:
            continue

        # Fit both models
        r1 = Rank1Quadratic()
        r1.fit(S[tr], T[tr], A[tr], Y[tr])

        fq = FullQuadratic()
        fq.fit(S[tr], T[tr], A[tr], Y[tr])

        Yp_r1 = r1.predict(S[te], T[te], A[te])
        Yp_fq = fq.predict(S[te], T[te], A[te])

        # ── Naive TOST ─────────────────────────────────────────────────
        logger.info("  Running naive TOST...")
        tost_results = tost_multiple_margins(Y[te], Yp_r1, Yp_fq, margins, alpha)

        for tr_ in tost_results:
            tr_["basin"] = basin_name
            tr_["method"] = "naive_tost"
            all_results.append(tr_)

        # ── Cluster bootstrap ──────────────────────────────────────────
        logger.info("  Running cluster bootstrap...")
        spatial_groups = assign_spatial_blocks(sub[te]).values

        boot_result = bootstrap_paired_ae_diff(
            Y[te], Yp_r1, Yp_fq,
            n_bootstrap=n_boot, confidence=1 - alpha,
            groups=spatial_groups, random_state=config["primary_seed"])

        if boot_result:
            boot_result["basin"] = basin_name
            boot_result["method"] = "cluster_bootstrap"
            all_results.append(boot_result)
            logger.info(f"    Bootstrap CI: [{boot_result['ci_lo']:.4f}, {boot_result['ci_hi']:.4f}]")

    # ── Save ────────────────────────────────────────────────────────────
    df_results = pd.DataFrame(all_results)
    df_results.to_csv(output_dir / "equivalence_analysis.csv", index=False)
    logger.info(f"\nSaved equivalence_analysis.csv")

    logger.info("\nStep 05 complete ✓")


if __name__ == "__main__":
    main()
