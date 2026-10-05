"""
Analyze symbolic regression results: Pareto frontiers, family frequencies,
operator-space dependence, and rank-1 competitiveness.

Run AFTER run_symbolic_search.py completes.
"""

import sys
import os
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.config_loader import load_config, get_path
from src.analysis.pareto import compute_pareto_frontier
from src.models.rank1_quadratic import Rank1Quadratic
from src.models.full_quadratic import FullQuadratic
from src.models.linear import LinearModel
from src.models.cubic import CubicPolynomial
from src.models.mean_baseline import MeanBaseline
from src.analysis.metrics import compute_metrics

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def load_symbolic_results():
    """Load the symbolic regression candidate library."""
    results_dir = get_path("results_raw") / "symbolic_regression"
    path = results_dir / "symbolic_candidates_all.csv"
    if not path.exists():
        logger.error(f"Results not found: {path}")
        return None
    df = pd.read_csv(path)
    logger.info(f"Loaded {len(df)} candidates from {path}")
    return df


def load_reference_models(basin_name, df_data, splits):
    """Fit reference models and return their validation metrics."""
    from src.basins import assign_basin

    if "basin" not in df_data.columns or df_data["basin"].isna().all():
        df_data["basin"] = assign_basin(df_data)

    mask = (df_data["basin"] == basin_name).values
    if mask.sum() == 0:
        return {}

    sub = df_data[mask]
    tr = splits["train"][mask]
    te = splits["validation"][mask]

    if tr.sum() == 0 or te.sum() == 0:
        return {}

    S, T, A, Y = sub["salinity"].values, sub["temperature"].values, \
                  sub["aou"].values, sub["tco2"].values

    models = {
        "Mean": MeanBaseline(),
        "Linear": LinearModel(),
        "Rank-1 Quad": Rank1Quadratic(),
        "Full Quad": FullQuadratic(),
        "Cubic": CubicPolynomial(),
    }

    results = {}
    for name, model in models.items():
        try:
            model.fit(S[tr], T[tr], A[tr], Y[tr])
            y_pred = model.predict(S[te], T[te], A[te])
            m = compute_metrics(Y[te], y_pred)
            results[name] = {
                "RMSE": m["RMSE"],
                "complexity": model.complexity,
                "parameters": model.parameter_count,
                "category": model.category,
            }
        except Exception as e:
            logger.warning(f"  {name} failed for {basin_name}: {e}")

    return results


def plot_pareto_with_reference(df_sr, ref_models, basin_name, output_dir):
    """Plot Pareto frontier with reference model overlay."""
    fig, ax = plt.subplots(figsize=(10, 7))

    # All symbolic candidates
    ax.scatter(df_sr["complexity"], df_sr["RMSE_val"],
               alpha=0.3, s=25, color="gray", label="SR candidates", zorder=2)

    # Pareto frontier
    pareto = compute_pareto_frontier(df_sr, "complexity", "RMSE_val")
    pareto = pareto.sort_values("complexity")
    ax.plot(pareto["complexity"], pareto["RMSE_val"],
            "o-", color="red", markersize=8, linewidth=2.5,
            label="Pareto frontier", zorder=3)

    # Reference models
    colors = {"Mean": "black", "Linear": "blue", "Rank-1 Quad": "green",
              "Full Quad": "orange", "Cubic": "purple"}
    markers = {"Mean": "s", "Linear": "^", "Rank-1 Quad": "D",
               "Full Quad": "v", "Cubic": "p"}

    for name, m in ref_models.items():
        if "RMSE" in m:
            ax.scatter(m["complexity"], m["RMSE"],
                       s=150, marker=markers.get(name, "o"),
                       color=colors.get(name, "gray"),
                       edgecolors="black", linewidth=1.5,
                       label=f'{name} (p={m["parameters"]})', zorder=4)

    ax.set_xlabel("Expression Complexity", fontsize=13)
    ax.set_ylabel("Validation RMSE (μmol/kg)", fontsize=13)
    ax.set_title(f"Pareto Frontier: {basin_name}", fontsize=15)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_dir / f"pareto_{basin_name.lower().replace(' ', '_')}.png", dpi=300)
    plt.close(fig)
    logger.info(f"  Saved Pareto plot for {basin_name}")


def analyze_family_frequency(df):
    """Analyze equation family recurrence across seeds and search spaces."""
    results = []
    for basin in df["basin"].unique():
        for search in df["search_space"].unique():
            subset = df[(df["basin"] == basin) & (df["search_space"] == search)]
            if len(subset) == 0:
                continue
            # Count families per seed
            seed_families = {}
            for seed in subset["seed"].unique():
                seed_sub = subset[subset["seed"] == seed]
                for family, count in seed_sub["family"].value_counts().items():
                    if family not in seed_families:
                        seed_families[family] = []
                    seed_families[family].append(count)

            for family, counts in seed_families.items():
                results.append({
                    "basin": basin,
                    "search_space": search,
                    "family": family,
                    "n_seeds_present": len(counts),
                    "mean_count": np.mean(counts),
                    "total_count": sum(counts),
                })

    return pd.DataFrame(results)


def main():
    logger.info("=" * 70)
    logger.info("PARETO FRONTIER AND FAMILY ANALYSIS")
    logger.info("=" * 70)

    config = load_config()
    output_dir = get_path("results_processed")
    fig_dir = get_path("results_figures")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)

    # Load SR results
    df_sr = load_symbolic_results()
    if df_sr is None or len(df_sr) == 0:
        logger.error("No symbolic regression results found. Run run_symbolic_search.py first.")
        return

    # Load data for reference models
    data_dir = get_path("processed_data")
    df_data = pd.read_csv(data_dir / "glodap_qc_filtered.csv")
    df_data["basin"] = df_data.get("basin", None)
    if df_data["basin"].isna().all():
        from src.basins import assign_basin
        df_data["basin"] = assign_basin(df_data, config.get("basins"))

    from src.splits import split_temporal
    splits = split_temporal(df_data,
                            train_end=config["temporal_split"]["train_end_year"],
                            val_end=config["temporal_split"]["validation_end_year"])

    # ── Per-basin Pareto analysis ────────────────────────────────────────
    for basin_name in df_sr["basin"].unique():
        logger.info(f"\n{'─' * 60}")
        logger.info(f"BASIN: {basin_name}")

        basin_df = df_sr[df_sr["basin"] == basin_name]

        # Reference models
        ref_models = load_reference_models(basin_name, df_data, splits)
        logger.info("  Reference models:")
        for name, m in ref_models.items():
            if "RMSE" in m:
                logger.info(f"    {name}: RMSE={m['RMSE']:.2f}, complexity={m['complexity']}")

        # Best SR candidate
        best_idx = basin_df["RMSE_val"].idxmin()
        best = basin_df.loc[best_idx]
        logger.info(f"  Best SR: RMSE={best['RMSE_val']:.2f}, complexity={best['complexity']}")
        logger.info(f"    Equation: {best['equation']}")

        # Plot Pareto frontier with reference models
        plot_pareto_with_reference(basin_df, ref_models, basin_name, fig_dir)

    # ── Family frequency analysis ───────────────────────────────────────
    logger.info(f"\n{'─' * 60}")
    logger.info("Family frequency analysis...")
    df_freq = analyze_family_frequency(df_sr)
    df_freq.to_csv(output_dir / "equation_family_frequency.csv", index=False)
    logger.info(f"  Saved equation_family_frequency.csv ({len(df_freq)} entries)")

    # ── Combined Pareto plot ─────────────────────────────────────────────
    logger.info("Generating combined Pareto plot...")
    basins = df_sr["basin"].unique()
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for i, basin_name in enumerate(basins[:4]):
        ax = axes[i]
        basin_df = df_sr[df_sr["basin"] == basin_name]

        ax.scatter(basin_df["complexity"], basin_df["RMSE_val"],
                   alpha=0.3, s=20, color="gray")

        pareto = compute_pareto_frontier(basin_df, "complexity", "RMSE_val")
        pareto = pareto.sort_values("complexity")
        ax.plot(pareto["complexity"], pareto["RMSE_val"],
                "o-", color="red", markersize=6, linewidth=2)

        # Rank-1 reference line
        ref = load_reference_models(basin_name, df_data, splits)
        if "Rank-1 Quad" in ref and "RMSE" in ref["Rank-1 Quad"]:
            ax.axhline(y=ref["Rank-1 Quad"]["RMSE"], color="green",
                       linestyle="--", alpha=0.7, label="Rank-1 quad")
        if "Full Quad" in ref and "RMSE" in ref["Full Quad"]:
            ax.axhline(y=ref["Full Quad"]["RMSE"], color="orange",
                       linestyle="--", alpha=0.7, label="Full quad")

        ax.set_xlabel("Complexity")
        ax.set_ylabel("Validation RMSE")
        ax.set_title(basin_name)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle("Symbolic Regression Pareto Frontiers", fontsize=15)
    plt.tight_layout()
    fig.savefig(fig_dir / "fig1_pareto_frontier_combined.png", dpi=300)
    plt.close(fig)
    logger.info("  Saved fig1_pareto_frontier_combined.png")

    logger.info("\nAnalysis complete ✓")


if __name__ == "__main__":
    main()
