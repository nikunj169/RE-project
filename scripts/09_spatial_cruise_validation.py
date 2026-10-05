"""Run leakage-safe spatial and cruise validation for the final model hierarchy."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.analysis.metrics import compute_metrics
from src.basins import assign_basin
from src.config_loader import get_path, load_config
from src.models.frozen_symbolic import FrozenSymbolicModel
from src.models.full_quadratic import FullQuadratic
from src.models.linear import LinearModel
from src.models.mean_baseline import MeanBaseline
from src.models.rank1_quadratic import Rank1Quadratic
from src.models.rank2_quadratic import Rank2Quadratic
from src.splits import (
    assign_cruise_blocks,
    assign_spatial_blocks,
    non_external_mask,
    validate_group_disjointness,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

BASINS = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]
MODEL_FACTORIES = {
    "linear": LinearModel,
    "rank1": Rank1Quadratic,
    "rank2": Rank2Quadratic,
    "full_quadratic": FullQuadratic,
}


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_frozen_candidates(df, candidates):
    """Select one symbolic candidate per basin using temporal validation only."""
    selected = []
    for basin in BASINS:
        subset = candidates[candidates["basin"] == basin].copy()
        if subset.empty:
            raise ValueError(f"No symbolic candidates available for {basin}")
        subset = subset.sort_values(
            ["RMSE_val", "complexity", "search_space", "seed", "candidate_idx"]
        )
        row = subset.iloc[0]
        selected.append(
            {
                "basin": basin,
                "candidate_id": (
                    f"{row['basin']}|{row['search_space']}|{int(row['seed'])}|"
                    f"{int(row['candidate_idx'])}"
                ),
                "search_space": row["search_space"],
                "seed": int(row["seed"]),
                "candidate_idx": int(row["candidate_idx"]),
                "equation": row["equation"],
                "sympy_format": row.get("sympy_format", row["equation"]),
                "complexity": int(row["complexity"]),
                "selection_rule": "minimum temporal validation RMSE; ties by complexity then identifiers",
                "selection_data": "year 2015-2017 validation only",
                "external_holdout_locked": True,
            }
        )
    return selected


def build_model(name, selected_row=None):
    if name == "frozen_symbolic":
        return FrozenSymbolicModel(
            equation=selected_row["equation"],
            sympy_format=selected_row.get("sympy_format"),
            complexity=selected_row.get("complexity"),
        )
    return MODEL_FACTORIES[name]()


def run_blocked_validation(df, basin, block_type, model_names, selected_row, n_splits=5):
    subset = df[df["basin"] == basin].copy()
    subset = subset[subset["year"] < 2018].copy()
    if subset.empty:
        raise ValueError(f"No eligible pre-2018 rows for {basin}")
    if subset["year"].max() >= 2018:
        raise AssertionError("External holdout leakage in blocked-validation frame")

    if block_type == "spatial":
        groups = assign_spatial_blocks(subset).to_numpy()
    elif block_type == "cruise":
        groups = assign_cruise_blocks(subset).to_numpy()
    else:
        raise ValueError(f"Unknown block type: {block_type}")

    n_groups = len(np.unique(groups))
    actual_splits = min(n_splits, n_groups)
    if actual_splits < 2:
        raise ValueError(f"Not enough {block_type} groups for {basin}: {n_groups}")

    X_dummy = np.zeros((len(subset), 1))
    y = subset["tco2"].to_numpy(float)
    gkf = GroupKFold(n_splits=actual_splits)
    fold_rows = []
    prediction_rows = []

    for fold, (train_idx, test_idx) in enumerate(gkf.split(X_dummy, y, groups)):
        train_groups = groups[train_idx]
        test_groups = groups[test_idx]
        group_diag = validate_group_disjointness(train_groups, test_groups)
        train = subset.iloc[train_idx]
        test = subset.iloc[test_idx]
        if test["year"].max() >= 2018 or train["year"].max() >= 2018:
            raise AssertionError("External holdout appeared in blocked-validation fold")

        for model_name in model_names:
            model = build_model(model_name, selected_row)
            model.fit(
                train["salinity"].to_numpy(float),
                train["temperature"].to_numpy(float),
                train["aou"].to_numpy(float),
                train["tco2"].to_numpy(float),
            )
            prediction = model.predict(
                test["salinity"].to_numpy(float),
                test["temperature"].to_numpy(float),
                test["aou"].to_numpy(float),
            )
            metrics = compute_metrics(test["tco2"].to_numpy(float), prediction)
            fold_rows.append(
                {
                    "basin": basin,
                    "block_type": block_type,
                    "model": model_name,
                    "fold": fold,
                    "n_test": len(test),
                    "n_train": len(train),
                    "n_groups_total": n_groups,
                    "n_train_groups": group_diag["n_train_groups"],
                    "n_test_groups": group_diag["n_test_groups"],
                    "min_test_year": int(test["year"].min()),
                    "max_test_year": int(test["year"].max()),
                    **metrics,
                }
            )
            prediction_rows.append(
                pd.DataFrame(
                    {
                        "original_index": test.index.to_numpy(),
                        "basin": basin,
                        "block_type": block_type,
                        "fold": fold,
                        "group_id": test_groups,
                        "model": model_name,
                        "year": test["year"].to_numpy(),
                        "cruise": test["cruise"].to_numpy(),
                        "latitude": test["latitude"].to_numpy(),
                        "longitude": test["longitude"].to_numpy(),
                        "tco2": test["tco2"].to_numpy(),
                        "prediction": prediction,
                        "residual": test["tco2"].to_numpy() - prediction,
                    }
                )
            )

    fold_df = pd.DataFrame(fold_rows)
    prediction_df = pd.concat(prediction_rows, ignore_index=True)
    return fold_df, prediction_df


def summarize(fold_df):
    group_cols = ["basin", "block_type", "model"]
    summary = (
        fold_df.groupby(group_cols, as_index=False)
        .agg(
            n_folds=("fold", "nunique"),
            n_test=("n_test", "sum"),
            n_groups_total=("n_groups_total", "max"),
            RMSE_mean=("RMSE", "mean"),
            RMSE_sd=("RMSE", "std"),
            MAE_mean=("MAE", "mean"),
            MAE_sd=("MAE", "std"),
            R2_mean=("R2", "mean"),
            R2_sd=("R2", "std"),
            Bias_mean=("Bias", "mean"),
            Bias_sd=("Bias", "std"),
        )
    )
    return summary


def main():
    root = Path(__file__).resolve().parents[1]
    config = load_config()
    processed_path = get_path("processed_data") / "glodap_qc_filtered.csv"
    sr_path = get_path("results_raw") / "symbolic_regression" / "symbolic_candidates_all.csv"
    output_dir = get_path("results_processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(processed_path)
    df["basin"] = assign_basin(df, config["basins"])
    external_mask = non_external_mask(df, 2018)
    if external_mask.all():
        raise AssertionError("Expected post-2018 external rows were not found")
    candidates = pd.read_csv(sr_path)
    selected = select_frozen_candidates(df, candidates)
    selected_path = root / "config" / "selected_symbolic_candidate.yaml"
    selected_path.write_text(
        yaml.safe_dump(
            {
                "selection_data_external_excluded": True,
                "external_start_year": 2018,
                "candidates": selected,
            },
            sort_keys=False,
        )
    )

    all_folds = []
    all_predictions = []
    model_names = ["linear", "rank1", "rank2", "full_quadratic", "frozen_symbolic"]
    cruise_available = "cruise" in df.columns and df["cruise"].notna().all()
    validation_status = {
        "spatial": "performed",
        "cruise": "performed" if cruise_available else "not_performed_missing_reliable_cruise_identifier",
    }

    for basin in BASINS:
        selected_row = next(row for row in selected if row["basin"] == basin)
        for block_type in ["spatial", "cruise"] if cruise_available else ["spatial"]:
            fold_df, prediction_df = run_blocked_validation(
                df, basin, block_type, model_names, selected_row, n_splits=5
            )
            all_folds.append(fold_df)
            all_predictions.append(prediction_df)

    fold_df = pd.concat(all_folds, ignore_index=True)
    prediction_df = pd.concat(all_predictions, ignore_index=True)
    summary_df = summarize(fold_df)
    fold_df.to_csv(output_dir / "spatial_cruise_validation_per_fold.csv", index=False)
    prediction_df.to_csv(output_dir / "spatial_cruise_validation_predictions.csv", index=False)
    summary_df.to_csv(output_dir / "spatial_cruise_validation_summary.csv", index=False)

    manifest = {
        "external_start_year": 2018,
        "external_rows_excluded": int((~external_mask).sum()),
        "spatial_block_definition": "5 degree latitude x 10 degree longitude",
        "cruise_identifier": "cruise" if cruise_available else None,
        "validation_status": validation_status,
        "models": model_names,
        "candidate_selection": "temporal validation RMSE only; external holdout locked",
        "processed_data_sha256": sha256_file(processed_path),
        "selected_candidate_config": str(selected_path.relative_to(root)),
    }
    (output_dir / "spatial_cruise_validation_manifest.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False)
    )

    # Generalization summary combines existing temporal/external results with blocked summaries.
    temporal = pd.read_csv(output_dir / "model_family_comparison.csv")
    temporal = temporal.rename(
        columns={"Model": "model", "Basin": "basin", "RMSE_val": "temporal_RMSE", "MAE_val": "temporal_MAE", "R2_val": "temporal_R2", "Bias_val": "temporal_Bias", "RMSE_ext": "external_RMSE", "MAE_ext": "external_MAE", "R2_ext": "external_R2", "Bias_ext": "external_Bias"}
    )
    temporal = temporal[["basin", "model", "temporal_RMSE", "temporal_MAE", "temporal_R2", "temporal_Bias", "external_RMSE", "external_MAE", "external_R2", "external_Bias"]]
    general = temporal.merge(
        summary_df.pivot(index=["basin", "model"], columns="block_type", values="RMSE_mean").reset_index().rename(columns={"spatial": "spatial_RMSE", "cruise": "cruise_RMSE"}),
        on=["basin", "model"],
        how="left",
    )
    general.to_csv(output_dir / "VALIDATION_GENERALIZATION_SUMMARY.csv", index=False)
    general.to_markdown(output_dir / "VALIDATION_GENERALIZATION_SUMMARY.md", index=False)
    logger.info("Blocked validation and generalization summary complete")


if __name__ == "__main__":
    main()
