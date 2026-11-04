"""Publication figures for the research program (all from saved results; no refitting).

Usage: python -m research.visualization.figures
Writes PNG (300 dpi) and PDF to research/figures/.
"""

from __future__ import annotations

import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from research.data_io import ROOT, groups, load_config, load_data

OUT = ROOT / "research" / "figures"
RES = ROOT / "research" / "results"
BASINS = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]
COL = {"Atlantic": "#1f77b4", "Indian": "#2ca02c", "Pacific": "#d95f02", "Southern Ocean": "#7570b3"}
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "savefig.dpi": 300, "savefig.bbox": "tight"})


def mid(layer: str) -> float:
    a, b = (float(x) for x in layer.strip("(]").split(","))
    a = 25.6 if a == 0 else a
    b = 28.35 if b > 29 else b
    return (a + b) / 2


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)


def fig1_map(cfg):
    df = load_data(cfg)
    used = pd.concat([g for _, _, g in groups(df, cfg)])
    prof = used.drop_duplicates("profile_id")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [2.2, 1]})
    lon = np.where(prof.longitude > 180, prof.longitude - 360, prof.longitude)
    for b in BASINS:
        m = (prof.basin == b).to_numpy()
        ax[0].scatter(lon[m], prof.latitude[m], s=1.5, color=COL[b], label=b, rasterized=True)
    ax[0].set(xlabel="Longitude (°E)", ylabel="Latitude (°N)", xlim=(-180, 180), ylim=(-80, 70),
              title="(a) Profiles in qualifying 4°×4° cells (strict QC)")
    ax[0].legend(markerscale=6, fontsize=7, loc="lower left")
    yr = used.groupby(["year", "basin"]).size().unstack(fill_value=0)[BASINS]
    ax[1].stackplot(yr.index, yr.T.values / 1000, colors=[COL[b] for b in BASINS])
    ax[1].set(xlabel="Year", ylabel="Samples (thousands)", title="(b) Temporal coverage")
    save(fig, "fig1_study_area")


def fig2_profiles():
    d = pd.read_csv(RES / "e02_layer_trends_primary.csv")
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, b in zip(axs, BASINS):
        s = d[d.basin == b].copy()
        s["y"] = s.layer.map(mid)
        s = s.sort_values("y")
        ax.errorbar(s.beta_tco2, s.y - 0.012, xerr=[s.beta_tco2 - s.ci_lo_tco2, s.ci_hi_tco2 - s.beta_tco2],
                    fmt="o", color="k", ms=4, capsize=2, label="Observed DIC")
        ax.errorbar(s.beta_canth, s.y + 0.012, xerr=[s.beta_canth - s.ci_lo_canth, s.ci_hi_canth - s.beta_canth],
                    fmt="s", color=COL[b], ms=4, capsize=2, label="TRACE C$_{ant}$")
        ax.axvline(0, color="0.7", lw=0.8)
        ax.set(title=b, xlabel="Trend at fixed θ, S, AOU\n(µmol kg$^{-1}$ yr$^{-1}$)", xlim=(-0.15, 1.3))
    axs[0].invert_yaxis()
    axs[0].set_ylabel("Neutral density γ$^n$ (kg m$^{-3}$; layer midpoint)")
    axs[0].legend(fontsize=7, loc="lower right")
    save(fig, "fig2_trend_profiles")


def fig3_ratio(cfg):
    h = pd.read_csv(RES / "summary_h1.csv")
    fig, ax = plt.subplots(figsize=(6, 3.8))
    lo, hi = cfg["h1_band"]
    ax.axvspan(lo, hi, color="0.88", label="TRACE ±15 % band")
    ax.axvline(1, color="0.5", lw=0.8)
    for i, b in enumerate(BASINS):
        s = h[h.basin == b].copy()
        s["y"] = s.layer.map(mid) + (i - 1.5) * 0.02
        for _, r in s.iterrows():
            mk = "o" if r.verdict == "supported" else "X"
            ax.errorbar(r.R, r.y, xerr=[[r.R - r.R_ci_lo], [r.R_ci_hi - r.R]], fmt=mk, color=COL[b], ms=5, capsize=2)
        ax.plot([], [], "o", color=COL[b], label=b)
    ax.plot([], [], "X", color="k", label="contradicted (H1)")
    ax.invert_yaxis()
    ax.set(xscale="log", xlim=(0.5, 15), xlabel="R = β$_{DIC}$ / β$_{Cant}$ (log scale)",
           ylabel="γ$^n$ layer midpoint (kg m$^{-3}$)", title="H1: observed vs TRACE accumulation (γ$^n$ < 27.8)")
    ax.legend(fontsize=7, loc="lower right")
    save(fig, "fig3_ratio_h1")


def fig4_delta():
    rob = pd.read_csv(RES / "summary_robustness.csv")
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, b in zip(axs, BASINS):
        s = rob[rob.basin == b].copy()
        s["y"] = s.layer.map(mid)
        ax.axvspan(-s.detection_limit.mean(), s.detection_limit.mean(), color="0.9",
                   label="offset detection limit")
        ax.hlines(s.y, s.delta_min_variants, s.delta_max_variants, color=COL[b], lw=5, alpha=0.35,
                  label="range over variants")
        ax.errorbar(s.delta, s.y, xerr=[s.delta - s.delta_ci_lo, s.delta_ci_hi - s.delta], fmt="o",
                    color="k", ms=4, capsize=2, label="primary (95 % CI)")
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set(title=b, xlabel="δ = β$_{DIC}$ − β$_{Cant}$\n(µmol kg$^{-1}$ yr$^{-1}$)", xlim=(-0.4, 0.5))
    axs[0].invert_yaxis()
    axs[0].set_ylabel("γ$^n$ layer midpoint (kg m$^{-3}$)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.12), frameon=False)
    save(fig, "fig4_delta_robustness")


def fig5_h3():
    files = sorted(glob.glob(str(RES / "e03_temporal_strict_*.csv")))
    if not files:
        print("fig5 skipped: no strict E03 results")
        return
    t = pd.read_csv(files[0])
    t = t[(t.split == "external") & (t.model != "mean")]
    covs = [c for c in ["none", "year", "canth"] if c in t.covariate.unique()]
    lab = {"none": "static", "year": "+ linear time", "canth": "+ TRACE C$_{ant}$"}
    models = ["linear", "rank1", "rank2", "full_quadratic", "cubic"]
    fig, axs = plt.subplots(2, 4, figsize=(11, 5), sharex=True)
    w = 0.8 / len(covs)
    for j, b in enumerate(BASINS):
        s = t[t.basin == b]
        for k, c in enumerate(covs):
            v = s[s.covariate == c].set_index("model").loc[models]
            x = np.arange(len(models)) + (k - (len(covs) - 1) / 2) * w
            axs[0, j].bar(x, v.Bias, w, label=lab[c])
            axs[1, j].bar(x, v.RMSE, w, label=lab[c])
        axs[0, j].axhline(0, color="k", lw=0.6)
        axs[0, j].set_title(b)
        axs[1, j].set_xticks(range(len(models)), ["lin", "r1", "r2", "quad", "cub"])
    axs[0, 0].set_ylabel("External (≥2018) bias\n(µmol kg$^{-1}$)")
    axs[1, 0].set_ylabel("External RMSE\n(µmol kg$^{-1}$)")
    axs[0, 0].legend(fontsize=7)
    save(fig, "fig5_h3_holdout_bias")


def fig6_emlr():
    p = RES / "e05_emlr_comparison.csv"
    if not p.exists():
        return
    d = pd.read_csv(p)
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, b in zip(axs, BASINS):
        s = d[d.basin == b].copy()
        s["y"] = s.layer.map(mid)
        s = s.sort_values("y")
        ax.errorbar(s.beta_dic, s.y, xerr=[s.beta_dic - s.beta_dic_lo, s.beta_dic_hi - s.beta_dic],
                    fmt="o", color="k", ms=4, capsize=2, label="Observed DIC (1989–2012)")
        ax.plot(s.trace_rate_94_07, s.y, "s-", color=COL[b], ms=4, label="TRACE C$_{ant}$ 1994–2007")
        ax.errorbar(s.emlr_rate, s.y, xerr=[s.emlr_rate - s.emlr_ens_min, s.emlr_ens_max - s.emlr_rate],
                    fmt="D--", color="0.45", ms=4, capsize=2, label="eMLR(C*) 1994–2007 (ens. range)")
        ax.axvline(0, color="0.8", lw=0.8)
        ax.set(title=b, xlabel="Rate (µmol kg$^{-1}$ yr$^{-1}$)", xlim=(-0.1, 1.45))
    axs[0].invert_yaxis()
    axs[0].set_ylabel("γ$^n$ layer midpoint (kg m$^{-3}$)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.08), frameon=False)
    save(fig, "fig6_emlr_comparison")


def figS1_synthetic():
    d = pd.read_csv(RES / "e01_synthetic_recovery.csv")
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.2))
    for s, g in d.groupby("sigma_cruise"):
        ax[0].scatter(g.beta_true + 0.01 * s, g.bias, label=f"σ$_{{cruise}}$={s:g}", s=14)
        ax[1].scatter(g.beta_true + 0.01 * s, g.coverage_pct, s=14)
    ax[0].axhspan(-0.02, 0.02, color="0.9")
    ax[0].set(xlabel="True trend (µmol kg$^{-1}$ yr$^{-1}$)", ylabel="Mean bias", title="(a) Estimator bias")
    ax[1].axhspan(90, 99, color="0.9")
    ax[1].axhline(95, color="0.5", lw=0.8)
    ax[1].set(xlabel="True trend (µmol kg$^{-1}$ yr$^{-1}$)", ylabel="95 % CI coverage (%)",
              title="(b) Bootstrap coverage (100 reps)")
    ax[0].legend(fontsize=7)
    save(fig, "figS1_synthetic_recovery")


def main():
    cfg = load_config()
    fig1_map(cfg)
    fig2_profiles()
    fig3_ratio(cfg)
    fig4_delta()
    fig5_h3()
    fig6_emlr()
    figS1_synthetic()
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
