"""
Step 04: Run PySR Pareto-optimal search across basins and operator sets.

For each basin and each search space (A-D):
  1. Run PySR
  2. Extract the complete Pareto-optimal candidate library
  3. Evaluate each candidate on the validation set
  4. Save candidates to CSV

Multiple seeds per configuration for stability analysis.
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
from src.models.symbolic import SymbolicRegressionSearch, PYSR_AVAILABLE
from src.analysis.metrics import compute_metrics
from src.analysis.pareto import classify_equation_family

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def main():
    if not PYSR_AVAILABLE:
        logger.error("PySR not installed. Install with: pip install pysr && python -m pysr install")
        logger.error("Julia is required: https://julialang.org/downloads/")
        return

    logger.info("=" * 70)
    logger.info("STEP 04: PySR Pareto-Optimal Search")
    logger.info("=" * 70)

    config = load_config()
    pysr_config = load_config("pysr_config.yaml")
    output_dir = get_path("results_raw") / "symbolic_regression"
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    data_dir = get_path("processed_data")
    df = pd.read_csv(data_dir / "glodap_qc_filtered.csv")
    df["basin"] = assign_basin(df, config.get("basins"))
    splits = split_temporal(df,
                            train_end=config["temporal_split"]["train_end_year"],
                            val_end=config["temporal_split"]["validation_end_year"])

    seeds = config["random_seeds"]
    search_spaces = pysr_config["search_spaces"]

    all_candidates = []

    for basin_name in ["Atlantic", "Indian", "Pacific", "Southern Ocean"]:
        mask = (df["basin"] == basin_name).values
        sub = df[mask]
        tr = splits["train"][mask]
        te = splits["validation"][mask]

        S_train = sub["salinity"].values[tr]
        T_train = sub["temperature"].values[tr]
        A_train = sub["aou"].values[tr]
        Y_train = sub["tco2"].values[tr]

        S_val = sub["salinity"].values[te]
        T_val = sub["temperature"].values[te]
        A_val = sub["aou"].values[te]
        Y_val = sub["tco2"].values[te]

        subsample_max = pysr_config.get("common", {}).get("subsample_max", 50000)

        for search_name, search_cfg in search_spaces.items():
            logger.info(f"\n{'─' * 60}")
            logger.info(f"Basin: {basin_name}  Search: {search_name}")

            # Use primary seed for main search
            primary_seed = config["primary_seed"]
            logger.info(f"  Running with seed={primary_seed}...")

            try:
                sr = SymbolicRegressionSearch(search_cfg, random_state=primary_seed)
                sr.fit(S_train, T_train, A_train, Y_train,
                       subsample_max=subsample_max)

                candidates = sr.get_candidates()
                if candidates is not None and len(candidates) > 0:
                    # Evaluate each candidate on validation set
                    for idx, row in candidates.iterrows():
                        try:
                            y_pred = sr.predict_from_equation(idx, S_val, T_val, A_val)
                            m = compute_metrics(Y_val, y_pred)
                            rmse_val = m["RMSE"]
                            r2_val = m["R2"]
                            mae_val = m["MAE"]
                        except Exception:
                            rmse_val = np.nan
                            r2_val = np.nan
                            mae_val = np.nan

                        all_candidates.append({
                            "basin": basin_name,
                            "search_space": search_name,
                            "seed": primary_seed,
                            "candidate_idx": int(idx),
                            "complexity": int(row["complexity"]),
                            "loss": float(row["loss"]),
                            "equation": str(row["equation"]),
                            "family": classify_equation_family(str(row["equation"])),
                            "RMSE_val": rmse_val,
                            "R2_val": r2_val,
                            "MAE_val": mae_val,
                        })

                    logger.info(f"  Extracted {len(candidates)} candidates")
            except Exception as e:
                logger.error(f"  Search failed: {e}")

    # Save all candidates
    if all_candidates:
        df_all = pd.DataFrame(all_candidates)
        output_path = output_dir / "symbolic_candidates_all.csv"
        df_all.to_csv(output_path, index=False)
        logger.info(f"\nSaved {len(df_all)} candidates to {output_path}")

        # Per-basin files
        for basin_name in df_all["basin"].unique():
            basin_df = df_all[df_all["basin"] == basin_name]
            basin_path = output_dir / f"symbolic_candidates_{basin_name.lower().replace(' ', '_')}.csv"
            basin_df.to_csv(basin_path, index=False)

    logger.info("\nStep 04 complete ✓")


if __name__ == "__main__":
    main()
