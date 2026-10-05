"""
Step 01: Load GLODAP data, apply QC, assign basins, create temporal splits.
Produces the processed dataset used by all downstream analyses.
"""

import sys
import os
import logging

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.config_loader import load_config, get_path
from src.data_loader import load_and_qc, save_processed
from src.basins import assign_basin, verify_mutual_exclusivity, get_basin_counts
from src.splits import split_temporal, verify_no_overlap

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def main():
    logger.info("=" * 70)
    logger.info("STEP 01: Data Preprocessing")
    logger.info("=" * 70)

    config = load_config()

    # ── Load and QC ────────────────────────────────────────────────────────
    logger.info("Loading and filtering GLODAP data...")
    df, qc_stats = load_and_qc()

    logger.info(f"QC statistics: {qc_stats}")

    # ── Assign basins ──────────────────────────────────────────────────────
    logger.info("Assigning basins...")
    df["basin"] = assign_basin(df, config.get("basins"))

    # Verify mutual exclusivity
    exclusivity = verify_mutual_exclusivity(df, config.get("basins"))
    logger.info(f"Mutual exclusivity: {exclusivity['mutually_exclusive']}")
    if not exclusivity["mutually_exclusive"]:
        logger.error("Basin definitions are NOT mutually exclusive!")
        raise ValueError("Basin overlap detected")

    # ── Create temporal splits ─────────────────────────────────────────────
    logger.info("Creating temporal splits...")
    splits = split_temporal(df,
                            train_end=config["temporal_split"]["train_end_year"],
                            val_end=config["temporal_split"]["validation_end_year"])
    df["split"] = "train"
    df.loc[splits["validation"], "split"] = "validation"
    df.loc[splits["external"], "split"] = "external"

    # Verify no overlap
    overlap = verify_no_overlap(splits)
    logger.info(f"Split overlap check: {overlap}")

    # ── Basin counts ───────────────────────────────────────────────────────
    basin_counts = get_basin_counts(df)
    logger.info("Basin counts:")
    logger.info(f"\n{basin_counts.to_string(index=False)}")

    # ── Save processed data ────────────────────────────────────────────────
    output_dir = get_path("processed_data")
    os.makedirs(output_dir, exist_ok=True)

    output_path = output_dir / "glodap_qc_filtered.csv"
    df.to_csv(output_path, index=False)
    logger.info(f"Saved processed data: {output_path}")
    logger.info(f"  Total rows: {len(df):,}")

    # Save basin counts
    basin_counts.to_csv(output_dir / "basin_counts.csv", index=False)
    logger.info("Saved basin counts to basin_counts.csv")

    # Save QC stats
    import json
    with open(output_dir / "qc_stats.json", "w") as f:
        json.dump(qc_stats, f, indent=2)
    logger.info("Saved QC stats to qc_stats.json")

    logger.info("Step 01 complete ✓")


if __name__ == "__main__":
    main()
