"""
Complete symbolic-regression experiment.

Runs PySR across all basins × search spaces × seeds.
Saves complete candidate libraries, not just best equations.
Generates Pareto frontier data and equation family classifications.
"""

import sys
import os
import logging
import time
import json
import yaml
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.config_loader import load_config, get_path
from src.basins import assign_basin
from src.splits import split_temporal
from src.models.symbolic import SymbolicRegressionSearch, PYSR_AVAILABLE
from src.analysis.metrics import compute_metrics
from src.analysis.pareto import compute_pareto_frontier
from src.analysis.family_classification import classify_structure

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("logs/symbolic_search.log"),
    ])
logger = logging.getLogger(__name__)


def run_single_search(basin_name, search_name, search_cfg, seed,
                       S_train, T_train, A_train, Y_train,
                       S_val, T_val, A_val, Y_val,
                       subsample_max=30000):
    """Run a single PySR search and evaluate all candidates."""
    start_time = time.time()
    result = {
        "basin": basin_name,
        "search_space": search_name,
        "seed": seed,
        "status": "started",
    }

    try:
        sr = SymbolicRegressionSearch(search_cfg, random_state=seed)
        sr.fit(S_train, T_train, A_train, Y_train, subsample_max=subsample_max)

        runtime = time.time() - start_time
        candidates = sr.get_candidates()

        if candidates is None or len(candidates) == 0:
            result["status"] = "no_candidates"
            return result, None

        # Evaluate each candidate on validation set
        eval_results = []
        for idx, row in candidates.iterrows():
            try:
                y_pred = sr.predict_from_equation(idx, S_val, T_val, A_val)
                m = compute_metrics(Y_val, y_pred)
                eval_entry = {
                    "basin": basin_name,
                    "search_space": search_name,
                    "seed": seed,
                    "candidate_idx": int(idx),
                    "complexity": int(row["complexity"]),
                    "loss": float(row["loss"]),
                    "score": float(row.get("score", np.nan)),
                    "equation": str(row["equation"]),
                    "family": classify_structure(str(row["equation"])),
                    "RMSE_val": m["RMSE"],
                    "MAE_val": m["MAE"],
                    "R2_val": m["R2"],
                    "Bias_val": m["Bias"],
                    "runtime_sec": runtime,
                }
                if "sympy_format" in row.index:
                    eval_entry["sympy_format"] = str(row["sympy_format"])
                eval_results.append(eval_entry)
            except Exception as e:
                logger.debug(f"  Candidate {idx} evaluation failed: {e}")

        result["status"] = "success"
        result["n_candidates"] = len(eval_results)
        result["runtime_sec"] = runtime

        logger.info(f"  {basin_name}/{search_name}/seed={seed}: "
                     f"{len(eval_results)} candidates in {runtime:.0f}s")

        return result, pd.DataFrame(eval_results) if eval_results else None

    except Exception as e:
        result["status"] = "failed"
        result["error"] = str(e)
        result["runtime_sec"] = time.time() - start_time
        logger.error(f"  FAILED: {basin_name}/{search_name}/seed={seed}: {e}")
        return result, None


def main():
    if not PYSR_AVAILABLE:
        logger.error("PySR not installed!")
        return

    logger.info("=" * 70)
    logger.info("COMPLETE SYMBOLIC REGRESSION EXPERIMENT")
    logger.info("=" * 70)

    config = load_config()
    pysr_config = load_config("pysr_config.yaml")

    output_dir = get_path("results_raw") / "symbolic_regression"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    # Load data
    data_dir = get_path("processed_data")
    df = pd.read_csv(data_dir / "glodap_qc_filtered.csv")
    df["basin"] = assign_basin(df, config.get("basins"))
    splits = split_temporal(df,
                            train_end=config["temporal_split"]["train_end_year"],
                            val_end=config["temporal_split"]["validation_end_year"])

    # Use reduced search settings for feasibility
    # (Full settings would be: populations=30, population_size=50, niterations=100)
    # These settings still produce meaningful Pareto fronts
    search_spaces = pysr_config["search_spaces"]
    seeds = [42, 123, 456]  # Start with 3 seeds, expand later if needed
    subsample_max = 30000  # Subsample large training sets

    all_candidates = []
    all_run_logs = []

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

        for search_name, search_cfg_raw in search_spaces.items():
            # Use moderate settings for computational feasibility
            search_cfg = dict(search_cfg_raw)
            search_cfg["populations"] = min(search_cfg.get("populations", 30), 20)
            search_cfg["population_size"] = min(search_cfg.get("population_size", 50), 40)
            search_cfg["niterations"] = min(search_cfg.get("niterations", 100), 50)
            search_cfg["maxsize"] = min(search_cfg.get("maxsize", 30), 25)

            for seed in seeds:
                logger.info(f"\n{'─' * 60}")
                logger.info(f"Starting: {basin_name} | {search_name} | seed={seed}")

                run_result, candidates_df = run_single_search(
                    basin_name, search_name, search_cfg, seed,
                    S_train, T_train, A_train, Y_train,
                    S_val, T_val, A_val, Y_val,
                    subsample_max=subsample_max)

                all_run_logs.append(run_result)

                if candidates_df is not None:
                    all_candidates.append(candidates_df)

    # ── Save all results ─────────────────────────────────────────────────
    logger.info(f"\n{'=' * 70}")
    logger.info("Saving results...")

    if all_candidates:
        df_all = pd.concat(all_candidates, ignore_index=True)
        output_path = output_dir / "symbolic_candidates_all.csv"
        df_all.to_csv(output_path, index=False)
        logger.info(f"  Saved {len(df_all)} total candidates to {output_path}")

        # Per-basin files
        for basin_name in df_all["basin"].unique():
            basin_df = df_all[df_all["basin"] == basin_name]
            basin_path = output_dir / f"symbolic_candidates_{basin_name.lower().replace(' ', '_')}.csv"
            basin_df.to_csv(basin_path, index=False)

        # Per-search-space files
        for search_name in df_all["search_space"].unique():
            search_df = df_all[df_all["search_space"] == search_name]
            search_path = output_dir / f"symbolic_candidates_{search_name}.csv"
            search_df.to_csv(search_path, index=False)

        # ── Compute Pareto frontiers ─────────────────────────────────────
        logger.info("Computing Pareto frontiers...")
        pareto_results = []
        for basin_name in df_all["basin"].unique():
            for search_name in df_all["search_space"].unique():
                subset = df_all[(df_all["basin"] == basin_name) &
                                (df_all["search_space"] == search_name)]
                if len(subset) == 0:
                    continue
                pareto = compute_pareto_frontier(subset, "complexity", "RMSE_val")
                pareto = pareto.copy()
                pareto["pareto_basin"] = basin_name
                pareto["pareto_search"] = search_name
                pareto_results.append(pareto)

        if pareto_results:
            df_pareto = pd.concat(pareto_results, ignore_index=True)
            df_pareto.to_csv(output_dir / "pareto_frontiers.csv", index=False)
            logger.info(f"  Saved Pareto frontiers: {len(df_pareto)} entries")

        # ── Family frequency analysis ────────────────────────────────────
        logger.info("Computing equation family frequencies...")
        freq_results = []
        for basin_name in df_all["basin"].unique():
            for search_name in df_all["search_space"].unique():
                for seed in seeds:
                    subset = df_all[(df_all["basin"] == basin_name) &
                                    (df_all["search_space"] == search_name) &
                                    (df_all["seed"] == seed)]
                    if len(subset) == 0:
                        continue
                    family_counts = subset["family"].value_counts()
                    for family, count in family_counts.items():
                        freq_results.append({
                            "basin": basin_name,
                            "search_space": search_name,
                            "seed": seed,
                            "family": family,
                            "count": int(count),
                            "fraction": count / len(subset),
                        })

        if freq_results:
            df_freq = pd.DataFrame(freq_results)
            df_freq.to_csv(output_dir / "equation_family_frequency.csv", index=False)
            logger.info(f"  Saved equation family frequencies")

    # ── Save run log ─────────────────────────────────────────────────────
    df_log = pd.DataFrame(all_run_logs)
    df_log.to_csv(output_dir / "experiment_run_log.csv", index=False)

    # Save manifest
    manifest = {
        "experiment": "symbolic_regression_search",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "python_version": sys.version,
        "pysr_version": __import__("pysr").__version__,
        "n_basins": 4,
        "n_search_spaces": len(search_spaces),
        "n_seeds": len(seeds),
        "subsample_max": subsample_max,
        "seeds": seeds,
        "search_spaces": list(search_spaces.keys()),
        "total_candidates": len(df_all) if all_candidates else 0,
        "n_successful_runs": sum(1 for r in all_run_logs if r.get("status") == "success"),
        "n_failed_runs": sum(1 for r in all_run_logs if r.get("status") == "failed"),
    }
    with open(output_dir / "experiment_manifest.yaml", "w") as f:
        yaml.dump(manifest, f, default_flow_style=False)

    logger.info(f"\n{'=' * 70}")
    logger.info("SYMBOLIC REGRESSION EXPERIMENT COMPLETE")
    logger.info(f"  Total candidates: {manifest['total_candidates']}")
    logger.info(f"  Successful runs: {manifest['n_successful_runs']}")
    logger.info(f"  Failed runs: {manifest['n_failed_runs']}")
    logger.info(f"{'=' * 70}")


if __name__ == "__main__":
    main()
