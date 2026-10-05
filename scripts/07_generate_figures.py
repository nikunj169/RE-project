"""
Step 07: Generate publication-quality figures.

All figures are generated from saved analysis outputs — no hard-coded values.
"""

import sys
import os
import logging
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.config_loader import load_config, get_path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def figure_data_overview(output_dir, fig_dir):
    """Figure 1: QC data and basin distributions."""
    data_path = get_path("processed_data") / "glodap_qc_filtered.csv"
    if not data_path.exists():
        return
    df = pd.read_csv(data_path)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    basin_order = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]
    for basin, ax in zip(basin_order, axes.flat):
        sub = df[df["basin"] == basin]
        ax.hexbin(sub["temperature"], sub["salinity"], C=sub["tco2"], gridsize=40, mincnt=1, cmap="viridis")
        ax.set_title(f"{basin} (n={len(sub):,})")
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Salinity")
        ax.grid(alpha=0.2)
    fig.suptitle("QC-filtered GLODAP observations by basin", fontsize=14)
    plt.tight_layout()
    fig.savefig(fig_dir / "fig1_data_overview.png")
    plt.close(fig)


def figure_model_comparison(output_dir, fig_dir):
    """Figure 2: Model-family comparison across basins."""
    path = output_dir / "model_family_comparison.csv"
    if not path.exists():
        logger.warning(f"  Skipping Figure 2: {path} not found")
        return

    df = pd.read_csv(path)
    df = df.dropna(subset=["RMSE_val"])

    basins = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]
    models = df["Model"].unique()

    fig, axes = plt.subplots(1, 4, figsize=(16, 5), sharey=True)

    for i, basin in enumerate(basins):
        ax = axes[i]
        bdf = df[df["Basin"] == basin].sort_values("Parameters")
        if len(bdf) == 0:
            continue

        x = range(len(bdf))
        bars = ax.bar(x, bdf["RMSE_val"], color="steelblue", alpha=0.8, edgecolor="black", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(bdf["Model"].str.replace("Model\\d+_", "", regex=True),
                           rotation=45, ha="right", fontsize=8)
        ax.set_title(basin)
        ax.set_ylabel("Validation RMSE (μmol/kg)" if i == 0 else "")
        ax.grid(axis="y", alpha=0.3)

        # Annotate parameter counts
        for j, (_, row) in enumerate(bdf.iterrows()):
            ax.text(j, row["RMSE_val"] + 0.5, f'p={row["Parameters"]}',
                    ha="center", fontsize=7, color="gray")

    fig.suptitle("Model Family Comparison: Validation RMSE by Basin", fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(fig_dir / "fig2_model_comparison.png")
    plt.close(fig)
    logger.info("  Saved fig2_model_comparison.png")


def figure_pareto_frontier(output_dir, fig_dir):
    """Figure 1: Pareto frontier for symbolic regression candidates."""
    path = output_dir.parent / "raw" / "symbolic_regression" / "symbolic_candidates_all.csv"
    if not path.exists():
        logger.warning(f"  Skipping Figure 1: {path} not found")
        return

    df = pd.read_csv(path)
    basins = df["basin"].unique()

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, basin in enumerate(basins[:4]):
        ax = axes[i]
        bdf = df[df["basin"] == basin].dropna(subset=["RMSE_val"])

        if len(bdf) == 0:
            continue

        # Plot all candidates
        ax.scatter(bdf["complexity"], bdf["RMSE_val"],
                   alpha=0.3, s=20, color="gray", label="All candidates")

        # Highlight Pareto frontier
        from src.analysis.pareto import compute_pareto_frontier
        pareto = compute_pareto_frontier(bdf, "complexity", "RMSE_val")
        pareto = pareto.sort_values("complexity")
        ax.plot(pareto["complexity"], pareto["RMSE_val"],
                "o-", color="red", markersize=6, linewidth=2, label="Pareto frontier")

        # Mark rank-1 quadratic complexity
        ax.axvline(x=8, color="blue", linestyle="--", alpha=0.5, label="Rank-1 quadratic")

        ax.set_xlabel("Expression Complexity")
        ax.set_ylabel("Validation RMSE (μmol/kg)")
        ax.set_title(basin)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle("Symbolic Regression Pareto Frontier", fontsize=14)
    plt.tight_layout()
    fig.savefig(fig_dir / "fig1_pareto_frontier.png")
    plt.close(fig)
    logger.info("  Saved fig1_pareto_frontier.png")


def figure_depth_stratified(output_dir, fig_dir):
    """Figure 5: Depth-stratified RMSE and model improvement."""
    path = output_dir / "depth_stratified_comparison.csv"
    if not path.exists():
        logger.warning(f"  Skipping Figure 5: {path} not found")
        return

    df = pd.read_csv(path)
    basins = df["Basin"].unique()

    fig, axes = plt.subplots(1, len(basins), figsize=(5 * len(basins), 5), sharey=True)
    if len(basins) == 1:
        axes = [axes]

    for i, basin in enumerate(basins):
        ax = axes[i]
        bdf = df[df["Basin"] == basin]

        x = range(len(bdf))
        w = 0.35
        ax.bar([xi - w/2 for xi in x], bdf["RMSE_Rank1"], w,
               label="Rank-1 Quadratic", color="steelblue")
        ax.bar([xi + w/2 for xi in x], bdf["RMSE_FullQuad"], w,
               label="Full Quadratic", color="lightcoral")

        ax.set_xticks(x)
        ax.set_xticklabels(bdf["depth_bin"], rotation=45, ha="right", fontsize=8)
        ax.set_title(basin)
        ax.set_ylabel("RMSE (μmol/kg)" if i == 0 else "")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Depth-Stratified Performance", fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(fig_dir / "fig5_depth_stratified.png")
    plt.close(fig)
    logger.info("  Saved fig5_depth_stratified.png")


def figure_recurrence(output_dir, fig_dir):
    """Figure 4: rank-1 recurrence by basin, operator space, and seed."""
    path = get_path("results_tables") / "equation_family_recurrence.csv"
    if not path.exists():
        logger.warning(f"  Skipping recurrence figure: {path} not found")
        return
    df = pd.read_csv(path)
    pivot = df.pivot_table(index=["basin", "search_space"], columns="seed", values="rank1_present", aggfunc="max").astype(float)
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(pivot.to_numpy(), cmap="Greens", vmin=0, vmax=1, aspect="auto")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([f"{b} | {s}" for b, s in pivot.index], fontsize=8)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([str(s) for s in pivot.columns])
    ax.set_xlabel("Random seed")
    ax.set_ylabel("Basin | operator space")
    ax.set_title("Rank-1 recurrence across basin × operator space × seed")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            ax.text(j, i, "yes" if pivot.iloc[i, j] else "no", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label="Rank-1 present")
    plt.tight_layout()
    fig.savefig(fig_dir / "fig4_rank1_recurrence.png")
    plt.close(fig)


def figure_generalization(output_dir, fig_dir):
    """Figure 5: temporal, spatial, cruise, and external generalization."""
    path = output_dir / "VALIDATION_GENERALIZATION_SUMMARY.csv"
    if not path.exists():
        logger.warning(f"  Skipping generalization figure: {path} not found")
        return
    df = pd.read_csv(path)
    df = df[df["model"].isin(["linear", "rank1", "rank2", "full_quadratic", "frozen_symbolic"])]
    basins = df["basin"].unique()
    fig, axes = plt.subplots(1, len(basins), figsize=(5 * len(basins), 5), sharey=True)
    if len(basins) == 1:
        axes = [axes]
    metrics = [("temporal_RMSE", "Temporal"), ("spatial_RMSE", "Spatial"), ("cruise_RMSE", "Cruise"), ("external_RMSE", "External")]
    for i, basin in enumerate(basins):
        ax = axes[i]
        bdf = df[df.basin == basin].set_index("model")
        models = [m for m in ["linear", "rank1", "rank2", "full_quadratic", "frozen_symbolic"] if m in bdf.index]
        x = np.arange(len(models))
        width = 0.18
        for k, (column, label) in enumerate(metrics):
            if column in bdf.columns:
                values = bdf.loc[models, column]
                ax.bar(x + (k - 1.5) * width, values, width, label=label)
        ax.set_xticks(x)
        ax.set_xticklabels(models, rotation=45, ha="right", fontsize=8)
        ax.set_title(basin)
        ax.set_ylabel("RMSE (μmol/kg)" if i == 0 else "")
        ax.grid(axis="y", alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.suptitle("Generalization across temporal, spatial, cruise, and external regimes")
    plt.tight_layout()
    fig.savefig(fig_dir / "fig5_generalization.png")
    plt.close(fig)


def figure_southern_abyss(output_dir, fig_dir):
    """Figure 8: Southern Ocean abyssal global versus local fit."""
    path = output_dir / "southern_ocean_abyssal.csv"
    if not path.exists():
        return
    row = pd.read_csv(path).iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(["Global", "Local"], [row.global_r2, row.local_r2], color=["#9467bd", "#2ca02c"])
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_ylabel("R²")
    axes[0].set_title("Abyssal Southern Ocean R²")
    axes[1].bar(["Global", "Local"], [row.global_rmse, row.local_rmse], color=["#9467bd", "#2ca02c"])
    axes[1].set_ylabel("RMSE (μmol/kg)")
    axes[1].set_title("Abyssal Southern Ocean RMSE")
    plt.tight_layout()
    fig.savefig(fig_dir / "fig8_southern_abyss.png")
    plt.close(fig)


def figure_watermass(output_dir, fig_dir):
    """Figure 7: Atlantic water-mass RMSE and bias."""
    path = output_dir / "watermass_atlantic.csv"
    if not path.exists():
        return
    df = pd.read_csv(path)
    masses = df.water_mass.unique()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for model, sub in df.groupby("Model"):
        sub = sub.set_index("water_mass").reindex(masses)
        axes[0].plot(masses, sub.RMSE, marker="o", label=model)
        axes[1].plot(masses, sub.Bias, marker="o", label=model)
    axes[0].set_ylabel("RMSE (μmol/kg)")
    axes[1].set_ylabel("Bias (μmol/kg)")
    for ax in axes:
        ax.tick_params(axis="x", rotation=45)
        ax.grid(alpha=0.3)
        ax.legend()
    fig.suptitle("Atlantic performance by threshold-based water mass")
    plt.tight_layout()
    fig.savefig(fig_dir / "fig7_watermass_performance.png")
    plt.close(fig)


def figure_derivatives(output_dir, fig_dir):
    """Derivative analysis figure."""
    path = output_dir / "derivative_analysis.csv"
    if not path.exists():
        logger.warning(f"  Skipping derivative figure: {path} not found")
        return

    df = pd.read_csv(path)

    fig, ax = plt.subplots(figsize=(10, 6))

    basins = df["basin"].tolist()
    x = np.arange(len(basins))
    w = 0.25

    ax.bar(x - w, df["frac_dS_positive_pct"], w, label="dTCO2/dS > 0", color="blue", alpha=0.7)
    ax.bar(x, df["frac_dT_negative_pct"], w, label="dTCO2/dT < 0", color="red", alpha=0.7)
    ax.bar(x + w, df["frac_dA_positive_pct"], w, label="dTCO2/dA > 0", color="green", alpha=0.7)

    ax.set_xticks(x)
    ax.set_xticklabels(basins)
    ax.set_ylabel("Fraction with expected sign (%)")
    ax.set_title("Derivative Sign Analysis Across Observed Domain")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(95, 101)

    plt.tight_layout()
    fig.savefig(fig_dir / "fig_derivative_analysis.png")
    plt.close(fig)
    logger.info("  Saved fig_derivative_analysis.png")


def main():
    logger.info("=" * 70)
    logger.info("STEP 07: Generate Figures")
    logger.info("=" * 70)

    config = load_config()
    output_dir = get_path("results_processed")
    fig_dir = get_path("results_figures")
    os.makedirs(fig_dir, exist_ok=True)

    figure_data_overview(output_dir, fig_dir)
    figure_model_comparison(output_dir, fig_dir)
    figure_pareto_frontier(output_dir, fig_dir)
    figure_depth_stratified(output_dir, fig_dir)
    figure_recurrence(output_dir, fig_dir)
    figure_generalization(output_dir, fig_dir)
    figure_watermass(output_dir, fig_dir)
    figure_southern_abyss(output_dir, fig_dir)
    figure_derivatives(output_dir, fig_dir)

    logger.info("\nStep 07 complete ✓")


if __name__ == "__main__":
    main()
