"""Generate final model hierarchy, recurrence, Pareto, and provenance outputs."""

from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.config_loader import get_path, load_config
from src.analysis.pareto import compute_pareto_frontier

ROOT = Path(__file__).resolve().parents[1]
BASINS = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]
SEARCHES = ["SearchA_basic", "SearchB_explog", "SearchC_sqrt_square", "SearchD_broad"]
SEEDS = [42, 123, 456]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(command):
    try:
        return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT).strip()
    except Exception as exc:
        return f"unavailable: {exc}"


def aic_bic(y, prediction, parameter_count):
    residual = np.asarray(y, dtype=float) - np.asarray(prediction, dtype=float)
    n = len(residual)
    rss = float(np.sum(residual**2))
    sigma2 = max(rss / n, np.finfo(float).tiny)
    loglik = -0.5 * n * (np.log(2 * np.pi * sigma2) + 1)
    return {
        "AIC": -2 * loglik + 2 * parameter_count,
        "BIC": -2 * loglik + parameter_count * np.log(n),
    }


def recurrence_table(sr, pareto):
    rows = []
    for basin in BASINS:
        for search in SEARCHES:
            for seed in SEEDS:
                cell = sr[(sr.basin == basin) & (sr.search_space == search) & (sr.seed == seed)]
                r1 = cell[cell.family == "rank-1 quadratic"]
                best = cell.loc[cell.RMSE_val.idxmin()] if not cell.empty else None
                pareto_cell = pareto[(pareto.pareto_basin == basin) & (pareto.pareto_search == search)]
                pareto_r1 = pareto_cell[pareto_cell.family == "rank-1 quadratic"]
                rows.append(
                    {
                        "basin": basin,
                        "search_space": search,
                        "seed": seed,
                        "rank1_present": bool(len(r1)),
                        "rank1_n_candidates": int(len(r1)),
                        "rank1_best_RMSE_val": float(r1.RMSE_val.min()) if len(r1) else np.nan,
                        "rank1_best_complexity": int(r1.loc[r1.RMSE_val.idxmin(), "complexity"]) if len(r1) else np.nan,
                        "rank1_pareto_present": bool(len(pareto_r1)),
                        "best_unrestricted_RMSE_val": float(best.RMSE_val) if best is not None else np.nan,
                        "best_unrestricted_family": best.family if best is not None else None,
                        "rank1_minus_best_RMSE": float(r1.RMSE_val.min() - best.RMSE_val) if len(r1) and best is not None else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def pareto_summary(sr):
    rows = []
    rank1_rmse = {"Atlantic": 20.172, "Indian": 14.483, "Pacific": 16.002, "Southern Ocean": 10.575}
    for basin in BASINS:
        cell = sr[sr.basin == basin].dropna(subset=["RMSE_val"])
        pareto = compute_pareto_frontier(cell, "complexity", "RMSE_val")
        best = pareto.loc[pareto.RMSE_val.idxmin()]
        for tolerance in [0.01, 0.02, 0.05]:
            eligible = pareto[pareto.RMSE_val <= best.RMSE_val * (1 + tolerance)]
            simplest = eligible.loc[eligible.complexity.idxmin()]
            rows.append(
                {
                    "basin": basin,
                    "best_RMSE": best.RMSE_val,
                    "best_complexity": best.complexity,
                    "rank1_RMSE": rank1_rmse[basin],
                    "rank1_complexity": 8,
                    "tolerance": tolerance,
                    "simplest_within_tolerance_RMSE": simplest.RMSE_val,
                    "simplest_within_tolerance_complexity": simplest.complexity,
                    "rank1_within_tolerance": rank1_rmse[basin] <= best.RMSE_val * (1 + tolerance),
                }
            )
    return pd.DataFrame(rows)


def main():
    config = load_config()
    results_raw = get_path("results_raw") / "symbolic_regression"
    results_processed = get_path("results_processed")
    results_tables = get_path("results_tables")
    results_tables.mkdir(parents=True, exist_ok=True)

    sr = pd.read_csv(results_raw / "symbolic_candidates_all.csv")
    pareto = pd.read_csv(results_raw / "pareto_frontiers.csv")
    recurrence = recurrence_table(sr, pareto)
    recurrence.to_csv(results_tables / "equation_family_recurrence.csv", index=False)
    pareto_summary(sr).to_csv(results_tables / "pareto_summary.csv", index=False)

    mfc = pd.read_csv(results_processed / "model_family_comparison.csv")
    mfc["model"] = mfc["Model"].map({
        "Model0_mean": "mean",
        "Model1_linear": "linear",
        "Model2_rank1_quad": "rank1",
        "Model3_rank2_quad": "rank2",
        "Model4_full_quad": "full_quadratic",
        "Model5_cubic": "cubic",
    })
    hierarchy = mfc.copy()
    hierarchy["RMSE_delta_vs_rank1"] = np.nan
    hierarchy["RMSE_pct_delta_vs_rank1"] = np.nan
    for basin in hierarchy.Basin.unique():
        reference = hierarchy[(hierarchy.Basin == basin) & (hierarchy.model == "rank1")]["RMSE_val"].iloc[0]
        mask = hierarchy.Basin == basin
        hierarchy.loc[mask, "RMSE_delta_vs_rank1"] = hierarchy.loc[mask, "RMSE_val"] - reference
        hierarchy.loc[mask, "RMSE_pct_delta_vs_rank1"] = (hierarchy.loc[mask, "RMSE_val"] - reference) / reference * 100
    hierarchy.to_csv(results_processed / "parametric_model_hierarchy.csv", index=False)

    rich_hierarchy = results_processed / "final_model_hierarchy.csv"
    if rich_hierarchy.exists():
        rich = pd.read_csv(rich_hierarchy)
        if "AIC_train" in rich.columns:
            rich[["basin", "model", "parameter_count", "expression_complexity", "AIC_train", "BIC_train", "aic_bic_comparable"]].to_csv(
                results_processed / "model_information_criteria.csv", index=False
            )
    else:
        hierarchy.to_csv(rich_hierarchy, index=False)

    generalization_path = results_processed / "VALIDATION_GENERALIZATION_SUMMARY.csv"
    blocked_path = results_processed / "spatial_cruise_validation_summary.csv"
    if generalization_path.exists() and blocked_path.exists():
        model_map = {
            "Model0_mean": "mean",
            "Model1_linear": "linear",
            "Model2_rank1_quad": "rank1",
            "Model3_rank2_quad": "rank2",
            "Model4_full_quad": "full_quadratic",
            "Model5_cubic": "cubic",
        }
        temporal = pd.read_csv(results_processed / "model_family_comparison.csv")
        temporal["model"] = temporal["Model"].map(model_map)
        temporal = temporal.rename(columns={
            "Basin": "basin",
            "RMSE_val": "temporal_RMSE",
            "MAE_val": "temporal_MAE",
            "R2_val": "temporal_R2",
            "Bias_val": "temporal_Bias",
            "RMSE_ext": "external_RMSE",
            "MAE_ext": "external_MAE",
            "R2_ext": "external_R2",
            "Bias_ext": "external_Bias",
        })
        temporal = temporal[["basin", "model", "temporal_RMSE", "temporal_MAE", "temporal_R2", "temporal_Bias", "external_RMSE", "external_MAE", "external_R2", "external_Bias"]]
        rich = pd.read_csv(results_processed / "final_model_hierarchy.csv")
        frozen = rich[rich["model"] == "frozen_symbolic"].rename(columns={
            "validation_RMSE": "temporal_RMSE",
            "validation_MAE": "temporal_MAE",
            "validation_R2": "temporal_R2",
            "validation_Bias": "temporal_Bias",
            "external_RMSE": "external_RMSE",
            "external_MAE": "external_MAE",
            "external_R2": "external_R2",
            "external_Bias": "external_Bias",
        })[["basin", "model", "temporal_RMSE", "temporal_MAE", "temporal_R2", "temporal_Bias", "external_RMSE", "external_MAE", "external_R2", "external_Bias"]]
        temporal = pd.concat([temporal, frozen], ignore_index=True)
        blocked = pd.read_csv(blocked_path)
        blocked["metric_regime"] = blocked["block_type"] + "_RMSE"
        blocked_wide = blocked.pivot_table(index=["basin", "model"], columns="block_type", values="RMSE_mean").reset_index()
        blocked_wide = blocked_wide.rename(columns={"spatial": "spatial_RMSE", "cruise": "cruise_RMSE"})
        general = temporal.merge(blocked_wide, on=["basin", "model"], how="left")
        general.to_csv(generalization_path, index=False)
        general.to_markdown(results_processed / "VALIDATION_GENERALIZATION_SUMMARY.md", index=False)

    (results_processed / "final_model_hierarchy.md").write_text(
        pd.read_csv(rich_hierarchy).to_markdown(index=False) if rich_hierarchy.exists() else hierarchy.to_markdown(index=False)
    )

    raw_path = get_path("raw_data")
    provenance = {
        "project_root": str(ROOT),
        "raw_data": str(raw_path),
        "raw_data_sha256": sha256(raw_path) if raw_path.exists() else None,
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "julia": command_output(["julia", "--version"]),
        "pysr": command_output([sys.executable, "-c", "import pysr; print(pysr.__version__)"]),
        "external_holdout_start_year": config["temporal_split"]["validation_end_year"],
        "selection_data": "pre-2018 only; validation/model selection is 2015-2017",
        "external_holdout_locked": True,
        "symbolic_manifest_sha256": sha256(results_raw / "experiment_manifest.yaml"),
    }
    (ROOT / "results" / "REPRODUCIBILITY_FINAL.yaml").write_text(yaml.safe_dump(provenance, sort_keys=False))
    (ROOT / "results" / "REPRODUCIBILITY_FINAL.md").write_text(
        "# Reproducibility Final Record\n\n"
        "This record freezes the software, data, split, and holdout provenance for the final analysis.\n\n"
        + yaml.safe_dump(provenance, sort_keys=False)
        + "\nExact final-output command sequence:\n\n"
        "```bash\n"
        "python scripts/01_preprocess.py\n"
        "python scripts/02_baseline_reproduction.py\n"
        "python scripts/03_model_families.py\n"
        "python experiments/run_symbolic_search.py\n"
        "python scripts/09_spatial_cruise_validation.py\n"
        "python scripts/10_finalize_analysis_outputs.py\n"
        "python scripts/07_generate_figures.py\n"
        "python scripts/08_generate_tables.py\n"
        "```\n"
    )


if __name__ == "__main__":
    main()
