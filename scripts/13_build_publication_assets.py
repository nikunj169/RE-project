#!/usr/bin/env python3
"""Build publication figures and tables from saved analysis outputs."""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from src.models.cubic import CubicPolynomial
from src.models.full_quadratic import FullQuadratic
from src.models.linear import LinearModel
from src.models.mean_baseline import MeanBaseline
from src.models.rank1_quadratic import Rank1Quadratic
from src.models.rank2_quadratic import Rank2Quadratic

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "results" / "processed"
RAW_SR = ROOT / "results" / "raw" / "symbolic_regression"
TABLES_SRC = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures"
LATEX_FIG = ROOT / "final_latex" / "figures"
LATEX_TAB = ROOT / "final_latex" / "tables"
SUPP_FIG = ROOT / "final_latex" / "supplementary" / "supplementary_figures"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

BASINS = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]
MODEL_ORDER = ["mean", "linear", "rank1", "rank2", "full_quadratic", "cubic", "frozen_symbolic"]
MODEL_LABELS = {
    "mean": "Mean",
    "linear": "Linear",
    "rank1": "Rank-1",
    "rank2": "Rank-2",
    "full_quadratic": "Full quad.",
    "cubic": "Cubic",
    "frozen_symbolic": "Best SR",
}


def ensure_dirs():
    for d in (FIG_DIR, LATEX_FIG, LATEX_TAB, SUPP_FIG):
        d.mkdir(parents=True, exist_ok=True)


def savefig(fig, name):
    for dest in (FIG_DIR, LATEX_FIG):
        fig.savefig(dest / name)
    plt.close(fig)


def write_tex(path: Path, content: str):
    path.write_text(content)
    (LATEX_TAB / path.name).write_text(content)


def fig1_data_overview():
    df = pd.read_csv(ROOT / "data" / "processed" / "glodap_qc_filtered.csv")
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), sharex=False, sharey=False)
    vmin, vmax = 1900, 2400
    for ax, basin in zip(axes.flat, BASINS):
        sub = df[df["basin"] == basin]
        hb = ax.hexbin(
            sub["temperature"],
            sub["salinity"],
            C=sub["tco2"],
            reduce_C_function=np.mean,
            gridsize=45,
            mincnt=1,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(f"{basin} ($n={len(sub):,}$)")
        ax.set_xlabel("Temperature (°C)")
        ax.set_ylabel("Salinity")
    fig.colorbar(hb, ax=axes.ravel().tolist(), shrink=0.72, label=r"TCO$_2$ ($\mu$mol kg$^{-1}$)")
    fig.suptitle("QC-filtered GLODAP v2.2023 observations by basin", y=1.01)
    savefig(fig, "figure1_data_overview.png")


def fig2_hierarchy():
    df = pd.read_csv(PROC / "final_model_hierarchy.csv")
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 4.4), sharey=True)
    models = ["mean", "linear", "rank1", "rank2", "full_quadratic", "cubic", "frozen_symbolic"]
    colors = ["#9e9e9e", "#4c78a8", "#f58518", "#54a24b", "#e45756", "#b279a2", "#72b7b2"]
    for ax, basin in zip(axes, BASINS):
        b = df[df.basin == basin].set_index("model")
        vals = [b.loc[m, "validation_RMSE"] for m in models]
        ks = [int(b.loc[m, "parameter_count"]) if pd.notna(b.loc[m, "parameter_count"]) else None for m in models]
        bars = ax.bar(range(len(models)), vals, color=colors, edgecolor="black", linewidth=0.4)
        ax.set_xticks(range(len(models)))
        ax.set_xticklabels([MODEL_LABELS[m] for m in models], rotation=55, ha="right")
        ax.set_title(basin)
        if ax is axes[0]:
            ax.set_ylabel(r"Temporal-validation RMSE ($\mu$mol kg$^{-1}$)")
        for j, (v, k) in enumerate(zip(vals, ks)):
            lab = f"K={k}" if k is not None else "SR"
            ax.text(j, v + 1.2, lab, ha="center", fontsize=6.5, color="0.3")
    fig.tight_layout()
    savefig(fig, "figure2_model_hierarchy.png")


def fig3_pareto():
    df = pd.read_csv(RAW_SR / "symbolic_candidates_all.csv")
    hier = pd.read_csv(PROC / "final_model_hierarchy.csv")
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.6))
    from src.analysis.pareto import compute_pareto_frontier

    for ax, basin in zip(axes.flat, BASINS):
        bdf = df[df.basin == basin].dropna(subset=["RMSE_val", "complexity"])
        ax.scatter(bdf["complexity"], bdf["RMSE_val"], s=16, c="0.65", alpha=0.45, label="SR candidates")
        pareto = compute_pareto_frontier(bdf, "complexity", "RMSE_val").sort_values("complexity")
        ax.plot(pareto["complexity"], pareto["RMSE_val"], "o-", color="#d62728", ms=5, lw=1.5, label="Pareto front")
        h = hier[hier.basin == basin]
        r1 = h[h.model == "rank1"].iloc[0]
        ax.scatter([r1.expression_complexity], [r1.validation_RMSE], s=90, c="#f58518", zorder=5, label="Rank-1 (K=5, C=8)")
        ax.axvline(8, color="#f58518", ls="--", lw=0.8, alpha=0.6)
        ax.set_title(basin)
        ax.set_xlabel("PySR expression complexity")
        ax.set_ylabel(r"Validation RMSE ($\mu$mol kg$^{-1}$)")
        ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    savefig(fig, "figure3_pareto_frontier.png")


def fig4_recurrence():
    rec = pd.read_csv(TABLES_SRC / "equation_family_recurrence.csv")
    searches = ["SearchA_basic", "SearchB_explog", "SearchC_sqrt_square", "SearchD_broad"]
    seeds = [42, 123, 456]
    short = {"SearchA_basic": "A", "SearchB_explog": "B", "SearchC_sqrt_square": "C", "SearchD_broad": "D"}
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.6))
    for ax, basin in zip(axes, BASINS):
        mat = np.zeros((len(searches), len(seeds)))
        for i, s in enumerate(searches):
            for j, seed in enumerate(seeds):
                row = rec[(rec.basin == basin) & (rec.search_space == s) & (rec.seed == seed)]
                mat[i, j] = float(row.rank1_present.iloc[0])
        ax.imshow(mat, cmap="Greens", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(3))
        ax.set_xticklabels(seeds)
        ax.set_yticks(range(4))
        ax.set_yticklabels([short[s] for s in searches])
        ax.set_xlabel("Seed")
        ax.set_title(f"{basin}\n({int(mat.sum())}/12 cells)")
        for i in range(4):
            for j in range(3):
                ax.text(j, i, "yes" if mat[i, j] else "no", ha="center", va="center", fontsize=8,
                        color="white" if mat[i, j] else "0.2")
        if ax is axes[0]:
            ax.set_ylabel("Operator space")
    fig.suptitle("Rank-1 recurrence across 48 basin × operator-space × seed cells", y=1.05)
    fig.tight_layout()
    savefig(fig, "figure4_rank1_recurrence.png")


def fig5_generalization():
    g = pd.read_csv(PROC / "VALIDATION_GENERALIZATION_SUMMARY.csv")
    models = ["linear", "rank1", "rank2", "full_quadratic", "frozen_symbolic"]
    metrics = [
        ("temporal_RMSE", "Temporal val."),
        ("spatial_RMSE", "Spatial block"),
        ("cruise_RMSE", "Cruise block"),
        ("external_RMSE", "External ≥2018"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(12.8, 4.4), sharey=True)
    width = 0.18
    for ax, basin in zip(axes, BASINS):
        b = g[g.basin == basin].set_index("model")
        x = np.arange(len(models))
        for k, (col, lab) in enumerate(metrics):
            ax.bar(x + (k - 1.5) * width, b.loc[models, col], width, label=lab)
        ax.set_xticks(x)
        ax.set_xticklabels([MODEL_LABELS[m] for m in models], rotation=40, ha="right")
        ax.set_title(basin)
        if ax is axes[0]:
            ax.set_ylabel(r"RMSE ($\mu$mol kg$^{-1}$)")
            ax.legend(fontsize=6.5, loc="upper right")
    fig.tight_layout()
    savefig(fig, "figure5_generalization.png")


def fig6_depth():
    d = pd.read_csv(PROC / "depth_stratified_comparison.csv")
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 4.3), sharey=True)
    for ax, basin in zip(axes, BASINS):
        b = d[d.Basin == basin]
        x = np.arange(len(b))
        w = 0.38
        ax.bar(x - w / 2, b.RMSE_Rank1, w, label="Rank-1", color="#f58518")
        ax.bar(x + w / 2, b.RMSE_FullQuad, w, label="Full quadratic", color="#e45756")
        ax.set_xticks(x)
        ax.set_xticklabels(b.depth_bin, rotation=40, ha="right", fontsize=7)
        ax.set_title(basin)
        if ax is axes[0]:
            ax.set_ylabel(r"RMSE ($\mu$mol kg$^{-1}$)")
            ax.legend(fontsize=7)
    fig.tight_layout()
    savefig(fig, "figure6_depth_stratified.png")


def fig7_watermass():
    w = pd.read_csv(PROC / "watermass_atlantic.csv")
    masses = ["Surface", "Thermocline", "AAIW", "NADW", "AABW", "Unclassified"]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    for model, color in [("Rank1", "#f58518"), ("FullQuad", "#e45756")]:
        sub = w[w.Model == model].set_index("water_mass").reindex(masses)
        axes[0].plot(masses, sub.RMSE, "o-", color=color, label=model.replace("FullQuad", "Full quadratic").replace("Rank1", "Rank-1"))
        axes[1].plot(masses, sub.Bias, "o-", color=color, label=model.replace("FullQuad", "Full quadratic").replace("Rank1", "Rank-1"))
    axes[0].set_ylabel(r"RMSE ($\mu$mol kg$^{-1}$)")
    axes[1].set_ylabel(r"Bias ($\mu$mol kg$^{-1}$)")
    axes[1].axhline(0, color="k", lw=0.7)
    for ax in axes:
        ax.tick_params(axis="x", rotation=30)
        ax.legend(fontsize=8)
    fig.suptitle("Atlantic water-mass performance (temporal-validation subset)")
    fig.tight_layout()
    savefig(fig, "figure7_watermass.png")


def fig8_abyss():
    row = pd.read_csv(PROC / "southern_ocean_abyssal.csv").iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8))
    axes[0].bar(["Basin-scale\nrank-1", "Local abyssal\nrank-1"], [row.global_r2, row.local_r2], color=["#9467bd", "#2ca02c"])
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].set_ylabel(r"$R^2$")
    axes[0].set_title("Southern Ocean abyss $R^2$")
    axes[1].bar(["Basin-scale\nrank-1", "Local abyssal\nrank-1"], [row.global_rmse, row.local_rmse], color=["#9467bd", "#2ca02c"])
    axes[1].set_ylabel(r"RMSE ($\mu$mol kg$^{-1}$)")
    axes[1].set_title("Southern Ocean abyss RMSE")
    fig.tight_layout()
    savefig(fig, "figure8_southern_abyss.png")


def fig9_derivatives():
    d = pd.read_csv(PROC / "derivative_analysis.csv")
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    x = np.arange(len(d))
    w = 0.25
    ax.bar(x - w, d.frac_dS_positive_pct, w, label=r"$\partial\mathrm{TCO}_2/\partial S>0$", color="#4c78a8")
    ax.bar(x, d.frac_dT_negative_pct, w, label=r"$\partial\mathrm{TCO}_2/\partial T<0$", color="#e45756")
    ax.bar(x + w, d.frac_dA_positive_pct, w, label=r"$\partial\mathrm{TCO}_2/\partial\mathrm{AOU}>0$", color="#54a24b")
    ax.set_xticks(x)
    ax.set_xticklabels(d.basin)
    ax.set_ylabel("Fraction of observations (%)")
    ax.set_ylim(99.9, 100.05)
    ax.legend(fontsize=8)
    ax.set_title("Empirical derivative-sign fractions of the fitted rank-1 surface")
    fig.tight_layout()
    savefig(fig, "figure9_derivatives.png")


def fig10_obs_pred():
    df = pd.read_csv(ROOT / "data" / "processed" / "glodap_qc_filtered.csv")
    models = [
        ("mean", "Mean", MeanBaseline()),
        ("linear", "Linear", LinearModel()),
        ("rank1", "Rank-1", Rank1Quadratic()),
        ("rank2", "Rank-2", Rank2Quadratic()),
        ("full_quadratic", "Full quadratic", FullQuadratic()),
        ("cubic", "Cubic", CubicPolynomial()),
    ]
    validation = df[df.year >= 2015]
    latest_year = int(validation.year.max())
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.4), sharex=True, sharey=True)
    limits = [validation.tco2.min(), validation.tco2.max()]

    for ax, (_, label, model) in zip(axes.flat, models):
        observed = []
        predicted = []
        for basin in BASINS:
            train = df[(df.basin == basin) & (df.year < 2015)]
            val = validation[validation.basin == basin]
            model.fit(train.salinity, train.temperature, train.aou, train.tco2)
            observed.append(val.tco2.to_numpy())
            predicted.append(model.predict(val.salinity, val.temperature, val.aou))
        observed = np.concatenate(observed)
        predicted = np.concatenate(predicted)
        ax.hexbin(observed, predicted, gridsize=55, mincnt=1, cmap="magma", bins="log")
        ax.plot(limits, limits, "w--", lw=1.2)
        ax.set_title(label)
        ax.set_xlabel(r"Observed TCO$_2$ ($\mu$mol kg$^{-1}$)")
        ax.set_ylabel(r"Predicted TCO$_2$ ($\mu$mol kg$^{-1}$)")
        ax.set_xlim(limits)
        ax.set_ylim(limits)

    fig.suptitle(f"Equation-family predictions (2015–{latest_year})", y=1.01)
    fig.tight_layout()
    savefig(fig, "figure10_obs_vs_pred.png")


def fig11_yearwise_obs_pred():
    df = pd.read_csv(ROOT / "data" / "processed" / "glodap_qc_filtered.csv")
    models = [
        ("mean", "Mean", MeanBaseline()),
        ("linear", "Linear", LinearModel()),
        ("rank1", "Rank-1", Rank1Quadratic()),
        ("rank2", "Rank-2", Rank2Quadratic()),
        ("full_quadratic", "Full quadratic", FullQuadratic()),
        ("cubic", "Cubic", CubicPolynomial()),
    ]
    evaluation = df[df.year >= 2015]
    years = sorted(evaluation.year.astype(int).unique())
    limits = [evaluation.tco2.min(), evaluation.tco2.max()]

    for model_name, label, model in models:
        predictions = {}
        for basin in BASINS:
            train = df[(df.basin == basin) & (df.year < 2015)]
            model.fit(train.salinity, train.temperature, train.aou, train.tco2)
            for year in years:
                subset = evaluation[(evaluation.basin == basin) & (evaluation.year == year)]
                if subset.empty:
                    continue
                predictions.setdefault(year, []).append(
                    (subset.tco2.to_numpy(), model.predict(subset.salinity, subset.temperature, subset.aou))
                )

        ncols = 3
        nrows = int(np.ceil(len(years) / ncols))
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(13.2, 4.0 * nrows), sharex=True, sharey=True
        )
        axes = np.atleast_1d(axes).ravel()
        for ax, year in zip(axes, years):
            observed = np.concatenate([pair[0] for pair in predictions[year]])
            predicted = np.concatenate([pair[1] for pair in predictions[year]])
            ax.hexbin(observed, predicted, gridsize=55, mincnt=1, cmap="magma", bins="log")
            ax.plot(limits, limits, "w--", lw=1.2)
            ax.set_title(str(year))
            ax.set_xlabel(r"Observed TCO$_2$ ($\mu$mol kg$^{-1}$)")
            ax.set_ylabel(r"Predicted TCO$_2$ ($\mu$mol kg$^{-1}$)")
            ax.set_xlim(limits)
            ax.set_ylim(limits)
        for ax in axes[len(years):]:
            ax.set_visible(False)
        fig.suptitle(f"{label} equation: year-wise predictions (2015–{years[-1]})", y=1.01)
        fig.tight_layout()
        savefig(fig, f"figure11_yearwise_{model_name}_obs_vs_pred.png")


def compute_prediction_intervals():
    df = pd.read_csv(ROOT / "data" / "processed" / "glodap_qc_filtered.csv")
    coef = pd.read_csv(PROC / "baseline_coefficients.csv")
    rows = []
    for basin in BASINS:
        c = coef[coef.Basin == basin].iloc[0]
        b = df[df.basin == basin].copy()

        def pred(frame):
            z = c.alpha * frame.salinity + c.beta * frame.temperature + c.gamma * (frame.aou / 100.0) + c.delta
            return z ** 2 + c.epsilon

        train = b[b.year < 2015]
        val = b[(b.year >= 2015) & (b.year < 2018)]
        ext = b[b.year >= 2018]
        resid = train.tco2 - pred(train)
        lo, hi = np.quantile(resid, [0.025, 0.975])
        for split, frame in [("validation", val), ("external", ext)]:
            r = frame.tco2 - pred(frame)
            cover = ((r >= lo) & (r <= hi)).mean() * 100
            rows.append(
                {
                    "basin": basin,
                    "split": split,
                    "n": len(frame),
                    "q025_umol": lo,
                    "q975_umol": hi,
                    "interval_width": hi - lo,
                    "coverage_pct": cover,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(PROC / "prediction_interval_coverage.csv", index=False)
    return out


def tables(pi):
    sizes = pd.read_csv(ROOT / "data" / "processed" / "basin_counts.csv")
    write_tex(
        TABLES_SRC / "table1_sample_sizes.tex",
        r"""\begin{table}[htbp]
\centering
\caption{Quality-controlled GLODAP v2.2023 sample sizes after basin assignment. Training years are $<2015$; temporal validation years are $2015$--$2017$; the locked external holdout is $\mathrm{year}\ge 2018$.}
\label{tab:sample_sizes}
\begin{tabular}{lrrrr}
\toprule
Basin & $N$ & Training & Temporal validation & External holdout \\
\midrule
Atlantic & 147{,}448 & 125{,}762 & 11{,}436 & 10{,}250 \\
Indian & 40{,}391 & 32{,}498 & 3{,}764 & 4{,}129 \\
Pacific & 187{,}443 & 138{,}202 & 31{,}791 & 17{,}450 \\
Southern Ocean & 97{,}248 & 82{,}045 & 6{,}483 & 8{,}720 \\
Total & 472{,}530 & 378{,}507 & 53{,}474 & 40{,}549 \\
\bottomrule
\end{tabular}
\end{table}
""",
    )

    write_tex(
        TABLES_SRC / "table2_model_families.tex",
        r"""\begin{table}[htbp]
\centering
\caption{Model families compared in this study. Fitted parameter count $K$ is distinct from PySR expression complexity $C$.}
\label{tab:model_families}
\small
\begin{tabular}{llp{6.6cm}rr}
\toprule
Family & Form & Notes & $K$ & $C$ \\
\midrule
Mean & basin training mean & Null baseline & 1 & 1 \\
Linear & $\alpha S+\beta T+\gamma A_{100}+\delta$ & Additive surface & 4 & 4 \\
Rank-1 quadratic & $(\alpha S+\beta T+\gamma A_{100}+\delta)^2+\epsilon$ & Hessian $H=2\mathbf{v}\mathbf{v}^\top$ & 5 & 8 \\
Rank-2 quadratic & $z_1^2+z_2^2+\epsilon$ & Second quadratic direction & 9 & 16 \\
Full quadratic & all 10 monomials of degree $\le 2$ & Unconstrained second-order surface & 10 & 20 \\
Cubic & all monomials of degree $\le 3$ & Flexible polynomial benchmark & 20 & 30 \\
Best SR & basin-specific PySR expression & Frozen on temporal validation RMSE & --- & 15--25 \\
\bottomrule
\end{tabular}
\end{table}
""",
    )

    hier = pd.read_csv(PROC / "final_model_hierarchy.csv")
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Temporal-validation performance (years 2015--2017). Units of RMSE, MAE and bias are $\mu\mathrm{mol}\,\mathrm{kg}^{-1}$.}",
        r"\label{tab:temporal}",
        r"\small",
        r"\begin{tabular}{llrrrrr}",
        r"\toprule",
        r"Basin & Model & $N$ & RMSE & MAE & $R^2$ & Bias \\",
        r"\midrule",
    ]
    for basin in BASINS:
        b = hier[hier.basin == basin]
        first = True
        for _, r in b.iterrows():
            name = MODEL_LABELS.get(r.model, r.model)
            basin_cell = basin if first else ""
            first = False
            lines.append(
                f"{basin_cell} & {name} & {int(r.n_validation):,} & {r.validation_RMSE:.2f} & {r.validation_MAE:.2f} & {r.validation_R2:.3f} & {r.validation_Bias:+.2f} \\\\"
            )
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table3_temporal_validation.tex", "\n".join(lines))

    spat = pd.read_csv(PROC / "spatial_cruise_validation_summary.csv")
    keep = ["linear", "rank1", "rank2", "full_quadratic", "frozen_symbolic"]
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Blocked generalization on pre-2018 observations. Spatial blocks are $5^\circ$ latitude $\times$ $10^\circ$ longitude; cruise blocks use the GLODAP cruise identifier. RMSE is the mean across GroupKFold folds.}",
        r"\label{tab:spatial}",
        r"\small",
        r"\begin{tabular}{llllrrrr}",
        r"\toprule",
        r"Basin & Block & Model & Groups & $N$ & RMSE & $R^2$ & Bias \\",
        r"\midrule",
    ]
    for basin in BASINS:
        for btype in ["spatial", "cruise"]:
            sub = spat[(spat.basin == basin) & (spat.block_type == btype) & (spat.model.isin(keep))]
            first = True
            for _, r in sub.iterrows():
                lines.append(
                    f"{basin if first else ''} & {btype} & {MODEL_LABELS.get(r.model, r.model)} & {int(r.n_groups_total)} & {int(r.n_test):,} & {r.RMSE_mean:.2f} & {r.R2_mean:.3f} & {r.Bias_mean:+.2f} \\\\"
                )
                first = False
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table4_spatial_cruise.tex", "\n".join(lines))

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Locked external holdout (year $\ge 2018$). This split was not used for candidate search, operator choice, or model selection.}",
        r"\label{tab:external}",
        r"\small",
        r"\begin{tabular}{llrrrrr}",
        r"\toprule",
        r"Basin & Model & $N$ & RMSE & MAE & $R^2$ & Bias \\",
        r"\midrule",
    ]
    for basin in BASINS:
        b = hier[hier.basin == basin]
        first = True
        for _, r in b.iterrows():
            name = MODEL_LABELS.get(r.model, r.model)
            lines.append(
                f"{basin if first else ''} & {name} & {int(r.n_external):,} & {r.external_RMSE:.2f} & {r.external_MAE:.2f} & {r.external_R2:.3f} & {r.external_Bias:+.2f} \\\\"
            )
            first = False
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table5_external_holdout.tex", "\n".join(lines))

    write_tex(
        TABLES_SRC / "table6_rank_comparison.tex",
        r"""\begin{table}[htbp]
\centering
\caption{Rank-1 versus rank-2 and full quadratic models on temporal validation. $\Delta$RMSE is relative to rank-1 (negative values indicate improvement over rank-1). Units are $\mu\mathrm{mol}\,\mathrm{kg}^{-1}$.}
\label{tab:rank_compare}
\begin{tabular}{lrrrrrr}
\toprule
Basin & Rank-1 RMSE & Rank-2 RMSE & $\Delta$ rank-2 & Full-quad RMSE & $\Delta$ full quad & Rank-1 vs full quad (\%) \\
\midrule
Atlantic & 20.17 & 22.16 & $+1.98$ & 23.41 & $+3.23$ & $-13.8$ \\
Indian & 14.48 & 14.74 & $+0.25$ & 13.73 & $-0.75$ & $+5.5$ \\
Pacific & 16.00 & 15.71 & $-0.29$ & 15.50 & $-0.50$ & $+3.2$ \\
Southern Ocean & 10.58 & 10.54 & $-0.04$ & 10.46 & $-0.11$ & $+1.1$ \\
\bottomrule
\end{tabular}
\end{table}
""",
    )

    rec = pd.read_csv(TABLES_SRC / "equation_family_recurrence.csv")
    n_cells = int(rec.rank1_present.sum())
    write_tex(
        TABLES_SRC / "table7_recurrence.tex",
        r"""\begin{table}[htbp]
\centering
\caption{Symbolic-regression library and rank-1 recurrence. A cell is a basin $\times$ operator-space $\times$ seed combination. Rank-1 counts use the heuristic family classifier and are not algebraically deduplicated.}
\label{tab:recurrence}
\begin{tabular}{lr}
\toprule
Quantity & Value \\
\midrule
Completed searches & 48 \\
Extracted/evaluated candidates & 599 \\
Pareto-optimal entries & 252 \\
Rank-1 classified candidates & 86 \\
Cells with $\ge 1$ rank-1 candidate & 24/48 \\
Basins with rank-1 recovery & 4/4 \\
Atlantic cells with rank-1 & 3/12 (Search C only) \\
Indian cells with rank-1 & 7/12 \\
Pacific cells with rank-1 & 5/12 \\
Southern Ocean cells with rank-1 & 9/12 \\
\bottomrule
\end{tabular}
\end{table}
""".replace("24/48", f"{n_cells}/48"),
    )

    coef = pd.read_csv(PROC / "baseline_coefficients.csv")
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Basin-specific rank-1 coefficients fitted on training data (year $<2015$). The model is $\mathrm{TCO}_2=(\alpha S+\beta T+\gamma\mathrm{AOU}/100+\delta)^2+\epsilon$. These five numbers are fitted parameters, not PySR complexity.}",
        r"\label{tab:coefficients}",
        r"\begin{tabular}{lrrrrr}",
        r"\toprule",
        r"Basin & $\alpha$ & $\beta$ & $\gamma$ & $\delta$ & $\epsilon$ \\",
        r"\midrule",
    ]
    for _, r in coef.iterrows():
        lines.append(f"{r.Basin} & {r.alpha:.4f} & {r.beta:.4f} & {r.gamma:.4f} & {r.delta:.2f} & {r.epsilon:.2f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table8_rank1_coefficients.tex", "\n".join(lines))

    d = pd.read_csv(PROC / "depth_stratified_comparison.csv")
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Depth-stratified rank-1 and full-quadratic RMSE on the temporal-validation subset. Units are $\mu\mathrm{mol}\,\mathrm{kg}^{-1}$.}",
        r"\label{tab:depth}",
        r"\small",
        r"\begin{tabular}{llrrrrr}",
        r"\toprule",
        r"Basin & Depth & $N$ & Rank-1 RMSE & Rank-1 $R^2$ & Full-quad RMSE & Rank-1 bias \\",
        r"\midrule",
    ]
    for basin in BASINS:
        b = d[d.Basin == basin]
        first = True
        for _, r in b.iterrows():
            lines.append(
                f"{basin if first else ''} & {r.depth_bin} & {int(r.n):,} & {r.RMSE_Rank1:.2f} & {r.R2_Rank1:.3f} & {r.RMSE_FullQuad:.2f} & {r.Bias_Rank1:+.2f} \\\\"
            )
            first = False
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table9_depth_stratified.tex", "\n".join(lines))

    w = pd.read_csv(PROC / "watermass_atlantic.csv")
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Atlantic threshold-based water-mass performance on the temporal-validation subset.}",
        r"\label{tab:watermass}",
        r"\small",
        r"\begin{tabular}{llrrrrr}",
        r"\toprule",
        r"Water mass & Model & $N$ & RMSE & MAE & $R^2$ & Bias \\",
        r"\midrule",
    ]
    for mass in ["Surface", "Thermocline", "AAIW", "NADW", "AABW", "Unclassified"]:
        sub = w[w.water_mass == mass]
        first = True
        for _, r in sub.iterrows():
            lab = "Rank-1" if r.Model == "Rank1" else "Full quadratic"
            lines.append(
                f"{mass if first else ''} & {lab} & {int(r.n):,} & {r.RMSE:.2f} & {r.MAE:.2f} & {r.R2:.3f} & {r.Bias:+.2f} \\\\"
            )
            first = False
        lines.append(r"\addlinespace")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table10_watermass.tex", "\n".join(lines))

    eq = pd.read_csv(PROC / "equivalence_analysis.csv")
    boot = eq[eq.method == "cluster_bootstrap"]
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Dependence-aware paired absolute-error comparison of rank-1 versus full quadratic (cluster bootstrap, 1000 resamples of $5^\circ\times 10^\circ$ blocks). Positive mean differences indicate larger rank-1 absolute error. Intervals are 95\%. Units are $\mu\mathrm{mol}\,\mathrm{kg}^{-1}$.}",
        r"\label{tab:equivalence}",
        r"\begin{tabular}{lrrrr}",
        r"\toprule",
        r"Basin & Mean $\Delta|\mathrm{error}|$ & Bootstrap SD & 95\% CI low & 95\% CI high \\",
        r"\midrule",
    ]
    for _, r in boot.iterrows():
        lines.append(f"{r.basin} & {r.mean_diff:.3f} & {r.boot_std:.3f} & {r.ci_lo:.3f} & {r.ci_hi:.3f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table11_equivalence.tex", "\n".join(lines))

    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Rank-1 residual quantile prediction intervals (2.5th--97.5th percentiles of training residuals) evaluated on temporal validation and the locked external holdout. This is a nonparametric residual interval, not a model of observation dependence.}",
        r"\label{tab:pi}",
        r"\small",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Basin & Split & $N$ & Interval width & Empirical coverage (\%) \\",
        r"\midrule",
    ]
    for _, r in pi.iterrows():
        lines.append(
            f"{r.basin} & {r.split} & {int(r.n):,} & {r.interval_width:.2f} & {r.coverage_pct:.1f} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    write_tex(TABLES_SRC / "table12_prediction_intervals.tex", "\n".join(lines))

    write_tex(
        TABLES_SRC / "table13_pysr_config.tex",
        r"""\begin{table}[htbp]
\centering
\caption{PySR operator spaces actually used. Common settings: MSE loss, deterministic serial execution, 64-bit precision, training subsample of 30{,}000 observations per search, PySR 2.3.0, Julia 1.12.7.}
\label{tab:pysr}
\small
\begin{tabular}{llp{4.8cm}rrr}
\toprule
Space & Binary & Unary & Max size & Populations $\times$ size & Iterations \\
\midrule
A & $+,-,\times,\div$ & none & 30 & $30\times 50$ & 100 \\
B & $+,-,\times,\div$ & $\exp,\log$ & 30 & $30\times 50$ & 100 \\
C & $+,-,\times,\div$ & $\mathrm{square},\sqrt{\cdot},\exp,\log$ & 30 & $30\times 50$ & 100 \\
D & $+,-,\times,\div$ & $\mathrm{square},\sqrt{\cdot},\exp,\log,|\cdot|,\mathrm{cube},\sqrt[3]{\cdot}$ & 35 & $40\times 60$ & 150 \\
\bottomrule
\end{tabular}
\end{table}
""",
    )


def main():
    os.chdir(ROOT)
    import sys

    sys.path.insert(0, str(ROOT))
    ensure_dirs()
    fig1_data_overview()
    fig2_hierarchy()
    fig3_pareto()
    fig4_recurrence()
    fig5_generalization()
    fig6_depth()
    fig7_watermass()
    fig8_abyss()
    fig9_derivatives()
    fig10_obs_pred()
    fig11_yearwise_obs_pred()
    pi = compute_prediction_intervals()
    tables(pi)
    print("Publication assets written.")
    print(pi.to_string(index=False))


if __name__ == "__main__":
    main()
