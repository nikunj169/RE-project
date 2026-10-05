"""
Step 03: Fit all model families and compare.

Models:
  0 - Baseline mean
  1 - Linear regression
  2 - Rank-1 quadratic
  3 - Rank-2 quadratic
  4 - Full second-order polynomial
  5 - Cubic polynomial

Produces the model-family comparison table used throughout the paper.
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
from src.models.mean_baseline import MeanBaseline
from src.models.linear import LinearModel
from src.models.rank1_quadratic import Rank1Quadratic
from src.models.rank2_quadratic import Rank2Quadratic
from src.models.full_quadratic import FullQuadratic
from src.models.cubic import CubicPolynomial
from src.analysis.metrics import compute_metrics

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

MODEL_CLASSES = [
    ("Model0_mean", MeanBaseline),
    ("Model1_linear", LinearModel),
    ("Model2_rank1_quad", Rank1Quadratic),
    ("Model3_rank2_quad", Rank2Quadratic),
    ("Model4_full_quad", FullQuadratic),
    ("Model5_cubic", CubicPolynomial),
]


def main():
    logger.info("=" * 70)
    logger.info("STEP 03: Model Family Comparison")
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

    all_results = []
    all_coefficients = []

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

        for model_name, model_class in MODEL_CLASSES:
            logger.info(f"  Fitting {model_name}...")
            try:
                model = model_class()
                model.fit(S[tr], T[tr], A[tr], Y[tr])

                for parameter, value in model.params_.items():
                    all_coefficients.append({
                        "Basin": basin_name,
                        "Model": model_name,
                        "Model_name": model.name,
                        "Parameter": parameter,
                        "Coefficient": float(value),
                    })

                # Validation metrics
                Yp_val = model.predict(S[te], T[te], A[te])
                m_val = compute_metrics(Y[te], Yp_val)

                # External metrics
                Yp_ext = model.predict(S[ex], T[ex], A[ex])
                m_ext = compute_metrics(Y[ex], Yp_ext)

                desc = model.get_description()

                all_results.append({
                    "Basin": basin_name,
                    "Model": model_name,
                    "Model_name": desc["name"],
                    "Category": desc["category"],
                    "Parameters": desc["parameter_count"],
                    "Effective_rank": desc["effective_rank"],
                    "Complexity": desc["complexity"],
                    "N_val": int(te.sum()),
                    "N_ext": int(ex.sum()),
                    "RMSE_val": round(m_val["RMSE"], 3),
                    "MAE_val": round(m_val["MAE"], 3),
                    "R2_val": round(m_val["R2"], 4),
                    "Bias_val": round(m_val["Bias"], 3),
                    "RMSE_ext": round(m_ext["RMSE"], 3),
                    "MAE_ext": round(m_ext["MAE"], 3),
                    "R2_ext": round(m_ext["R2"], 4),
                    "Bias_ext": round(m_ext["Bias"], 3),
                })

                logger.info(f"    Val: RMSE={m_val['RMSE']:.2f}  R²={m_val['R2']:.4f}")
            except Exception as e:
                logger.error(f"    FAILED: {e}")
                all_results.append({
                    "Basin": basin_name,
                    "Model": model_name,
                    "Model_name": model_class.name,
                    "Error": str(e),
                })

    # ── Save ────────────────────────────────────────────────────────────
    df_results = pd.DataFrame(all_results)
    df_results.to_csv(output_dir / "model_family_comparison.csv", index=False)
    logger.info(f"\nSaved model_family_comparison.csv")

    df_coefficients = pd.DataFrame(all_coefficients)
    df_coefficients.to_csv(output_dir / "model_family_coefficients.csv", index=False)
    logger.info("Saved model_family_coefficients.csv")

    # ── Accuracy gain per parameter ─────────────────────────────────────
    logger.info("\nAccuracy gain per additional parameter:")
    for basin_name in ["Atlantic", "Indian", "Pacific", "Southern Ocean"]:
        basin_df = df_results[df_results["Basin"] == basin_name].copy()
        if len(basin_df) < 2:
            continue
        basin_df = basin_df.dropna(subset=["RMSE_val", "Parameters"])
        basin_df = basin_df.sort_values("Parameters")

        ref_rmse = basin_df.iloc[0]["RMSE_val"]
        ref_params = basin_df.iloc[0]["Parameters"]

        for _, row in basin_df.iloc[1:].iterrows():
            delta_rmse = ref_rmse - row["RMSE_val"]
            delta_params = row["Parameters"] - ref_params
            gain = delta_rmse / delta_params if delta_params > 0 else 0
            logger.info(f"  {basin_name}: {row['Model']}: "
                        f"+{delta_params:.0f} params → {-delta_rmse:+.2f} RMSE "
                        f"({gain:+.3f} per param)")

    logger.info("\nStep 03 complete ✓")


if __name__ == "__main__":
    main()
