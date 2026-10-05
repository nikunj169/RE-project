"""
Step 08: Generate publication-quality tables.

All tables are generated from saved result CSVs — no hard-coded values.
"""

import sys
import os
import logging
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.config_loader import load_config, get_path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def generate_latex_table(df, caption, label, columns=None, float_format=".3f"):
    """Convert a DataFrame to a LaTeX table."""
    if columns:
        df = df[columns]

    fmt = lambda x: f"{x:{float_format}}" if isinstance(x, float) else str(x)
    latex = df.to_latex(
        index=False,
        caption=caption,
        label=label,
        float_format=fmt,
        escape=True,
        column_format="l" + "r" * (len(df.columns) - 1),
    )
    latex = latex.replace(
        "\\begin{tabular}",
        "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}",
        1,
    )
    latex = latex.replace("\n\\end{tabular}", "\n\\end{tabular}%\n}", 1)
    return latex


def generate_longtable(df, caption, label, columns=None, float_format=".3f"):
    """Convert a DataFrame to a multipage longtable."""
    if columns:
        df = df[columns]
    fmt = lambda x: f"{x:{float_format}}" if isinstance(x, float) else str(x)
    return df.to_latex(
        index=False,
        caption=caption,
        label=label,
        float_format=fmt,
        escape=True,
        longtable=True,
        column_format="l" + "r" * (len(df.columns) - 1),
    )


def main():
    logger.info("=" * 70)
    logger.info("STEP 08: Generate Tables")
    logger.info("=" * 70)

    config = load_config()
    output_dir = get_path("results_processed")
    table_dir = get_path("results_tables")
    os.makedirs(table_dir, exist_ok=True)

    def write_table(filename, df, caption, label, columns=None, float_format=".3f"):
        latex = generate_latex_table(df, caption, label, columns, float_format)
        with open(table_dir / filename, "w") as f:
            f.write(latex)
        logger.info("  Wrote %s", filename)

    def write_longtable(filename, df, caption, label, columns=None, float_format=".3f"):
        latex = generate_longtable(df, caption, label, columns, float_format)
        with open(table_dir / filename, "w") as f:
            f.write(latex)
        logger.info("  Wrote %s", filename)

    # ── Table 1: Sample sizes ───────────────────────────────────────────
    path = output_dir / "basin_counts.csv"
    if not path.exists():
        path = get_path("processed_data") / "basin_counts.csv"
    if path.exists():
        write_table(
            "table1_sample_sizes.tex",
            pd.read_csv(path),
            "Sample sizes by basin and temporal split.",
            "tab:sample_sizes",
        )

    # ── Table 2: Predictor distributions ────────────────────────────────
    data_path = get_path("processed_data") / "glodap_qc_filtered.csv"
    if data_path.exists():
        data = pd.read_csv(data_path)
        rows = []
        for basin, sub in data.groupby("basin"):
            for variable in ["salinity", "temperature", "aou", "tco2", "depth"]:
                values = sub[variable].dropna()
                rows.append({
                    "Basin": basin,
                    "Variable": variable,
                    "N": len(values),
                    "Mean": values.mean(),
                    "SD": values.std(),
                    "Q25": values.quantile(0.25),
                    "Median": values.median(),
                    "Q75": values.quantile(0.75),
                })
        write_table(
            "table2_predictor_distributions.tex",
            pd.DataFrame(rows),
            "Predictor and target distributions by basin.",
            "tab:predictor_distributions",
        )

    # ── Table 3: PySR configuration ─────────────────────────────────────
    pysr_config = load_config("pysr_config.yaml")
    config_rows = []
    for name, values in pysr_config["search_spaces"].items():
        config_rows.append({
            "Search": name,
            "Binary operators": ", ".join(values.get("binary_operators", [])),
            "Unary operators": ", ".join(values.get("unary_operators", [])),
            "Max size": values.get("maxsize"),
            "Populations": values.get("populations"),
            "Population size": values.get("population_size"),
            "Iterations": values.get("niterations"),
        })
    write_table(
        "table3_pysr_configuration.tex",
        pd.DataFrame(config_rows),
        "Symbolic-regression search configurations.",
        "tab:pysr_configuration",
    )

    # ── Table 4: Pareto candidates ──────────────────────────────────────
    path = get_path("results_raw") / "symbolic_regression" / "pareto_frontiers.csv"
    if path.exists():
        df = pd.read_csv(path)
        cols = ["pareto_basin", "pareto_search", "seed", "complexity", "RMSE_val", "MAE_val", "R2_val", "family", "equation"]
        cols = [c for c in cols if c in df.columns]
        write_table(
            "table4_pareto_equations.tex",
            df[cols].sort_values(["pareto_basin", "complexity"]),
            "Pareto-optimal symbolic-regression candidates.",
            "tab:pareto_equations",
        )

    # ── Table 5: Model family comparison ────────────────────────────────
    path = output_dir / "model_family_comparison.csv"
    if path.exists():
        df = pd.read_csv(path)
        cols = ["Basin", "Model_name", "Parameters", "RMSE_val", "MAE_val",
                "R2_val", "RMSE_ext", "R2_ext"]
        cols = [c for c in cols if c in df.columns]
        latex = generate_latex_table(
            df[cols],
            caption="Model family comparison. Validation and external holdout metrics.",
            label="tab:model_comparison")
        with open(table_dir / "table5_model_comparison.tex", "w") as f:
            f.write(latex)
        logger.info("  Table 5: model comparison")

    # ── Table 9: Depth-stratified results ───────────────────────────────
    path = output_dir / "depth_stratified_comparison.csv"
    if path.exists():
        df = pd.read_csv(path)
        latex = generate_latex_table(
            df,
            caption="Depth-stratified model comparison by basin.",
            label="tab:depth_stratified")
        with open(table_dir / "table9_depth_stratified.tex", "w") as f:
            f.write(latex)
        logger.info("  Table 9: depth-stratified")

    # ── Tables 11 and 12: Equivalence analysis ─────────────────────────
    path = output_dir / "equivalence_analysis.csv"
    if path.exists():
        df = pd.read_csv(path)
        write_table(
            "table11_equivalence.tex",
            df,
            "Dependence-aware model-difference and equivalence results.",
            "tab:equivalence",
        )
        margins = df[df["method"] == "naive_tost"].copy()
        if not margins.empty:
            write_table(
                "table12_equivalence_margins.tex",
                margins,
                "Equivalence-margin sensitivity across basins.",
                "tab:equivalence_margins",
            )

    # ── Table 13: Derivative analysis ───────────────────────────────────
    path = output_dir / "derivative_analysis.csv"
    if path.exists():
        df = pd.read_csv(path)
        write_table(
            "table13_derivatives.tex",
            df,
            "Derivative sign analysis across the observed predictor domain.",
            "tab:derivatives",
            float_format=".2f",
        )

    # ── Final hierarchy and information criteria ────────────────────────
    path = output_dir / "final_model_hierarchy.csv"
    if path.exists():
        df = pd.read_csv(path)
        cols = ["basin", "model", "parameter_count", "expression_complexity", "validation_RMSE", "validation_MAE", "validation_R2", "external_RMSE", "external_MAE", "external_R2", "validation_delta_vs_rank1"]
        cols = [c for c in cols if c in df.columns]
        write_table("table6_model_hierarchy.tex", df[cols], "Final model hierarchy across temporal validation and external holdout.", "tab:model_hierarchy")

    path = output_dir / "spatial_cruise_validation_summary.csv"
    if path.exists():
        df = pd.read_csv(path)
        cols = ["basin", "block_type", "model", "n_groups_total", "RMSE_mean", "RMSE_sd", "R2_mean"]
        display = df[cols].rename(columns={
            "basin": "Basin",
            "block_type": "Block",
            "model": "Model",
            "n_groups_total": "Groups",
            "RMSE_mean": "RMSE",
            "RMSE_sd": "RMSE SD",
            "R2_mean": "R2",
        })
        write_longtable("table8_spatial_validation.tex", display, "Spatial- and cruise-blocked validation summary. Full fold diagnostics are in the CSV output.", "tab:spatial_validation")

    path = output_dir / "watermass_atlantic.csv"
    if path.exists():
        df = pd.read_csv(path)
        write_table("table10_watermass_results.tex", df, "Atlantic performance by threshold-based water mass.", "tab:watermass")

    path = table_dir / "equation_family_recurrence.csv"
    if path.exists():
        recurrence = pd.read_csv(path)
        cols = ["basin", "search_space", "seed", "rank1_present", "rank1_pareto_present", "rank1_n_candidates", "rank1_best_RMSE_val", "rank1_best_complexity"]
        display = recurrence[cols].rename(columns={
            "basin": "Basin",
            "search_space": "Search",
            "seed": "Seed",
            "rank1_present": "R1",
            "rank1_pareto_present": "Pareto",
            "rank1_n_candidates": "N",
            "rank1_best_RMSE_val": "R1 RMSE",
            "rank1_best_complexity": "R1 C",
        })
        write_longtable("table7_recurrence.tex", display, "Rank-1 recurrence by basin, operator space, and seed. R1 C denotes PySR expression complexity.", "tab:recurrence")

    path = table_dir / "pareto_summary.csv"
    if path.exists():
        write_table("table4_pareto_summary.tex", pd.read_csv(path), "Pareto-frontier summaries and tolerance thresholds.", "tab:pareto_summary")

    logger.info("\nStep 08 complete ✓")


if __name__ == "__main__":
    main()
