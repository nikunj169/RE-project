"""Round-2 publication figures (from saved results only; no refitting).

Usage: python -m research.visualization.figures_round2
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from research.data_io import groups, load_config, load_data
from research.visualization.figures import BASINS, COL, OUT, RES, mid, save


def fig1_map_coast(cfg):
    """Study area with coastlines (cartopy) + temporal coverage; replaces fig1 if cartopy is present."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    df = load_data(cfg)
    used = pd.concat([g for _, _, g in groups(df, cfg)])
    prof = used.drop_duplicates("profile_id")
    lon = np.where(prof.longitude > 180, prof.longitude - 360, prof.longitude)
    fig = plt.figure(figsize=(10.5, 3.8))
    ax = fig.add_axes([0.0, 0.08, 0.64, 0.84], projection=ccrs.Robinson(central_longitude=-150))
    ax.add_feature(cfeature.LAND, facecolor="0.85", zorder=0)
    ax.coastlines(lw=0.4)
    ax.set_global()
    for b in BASINS:
        m = (prof.basin == b).to_numpy()
        ax.scatter(lon[m], prof.latitude[m], s=1.2, color=COL[b], label=b, transform=ccrs.PlateCarree(),
                   rasterized=True)
    ax.set_title("(a) Profiles in qualifying 4°×4° cells (strict QC)")
    ax.legend(markerscale=6, fontsize=7, loc="lower left", bbox_to_anchor=(0.0, -0.05))
    ax2 = fig.add_axes([0.71, 0.18, 0.28, 0.68])
    yr = used.groupby(["year", "basin"]).size().unstack(fill_value=0)[BASINS]
    ax2.stackplot(yr.index, yr.T.values / 1000, colors=[COL[b] for b in BASINS])
    ax2.set(xlabel="Year", ylabel="Samples (thousands)", title="(b) Temporal coverage")
    save(fig, "fig1_study_area")


def fig7_differential():
    """E02 absolute excess vs E06 drift-immune excess (relative to deep reference)."""
    e02 = pd.read_csv(RES / "e02_layer_trends_primary.csv")
    e06 = pd.read_csv(RES / "e06_differential_trends.csv")
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, b in zip(axs, BASINS):
        a = e02[e02.basin == b].assign(y=lambda s: s.layer.map(mid))
        d = e06[e06.basin == b].assign(y=lambda s: s.layer.map(mid))
        ax.errorbar(a.delta, a.y - 0.015, xerr=[a.delta - a.delta_ci_lo, a.delta_ci_hi - a.delta], fmt="o",
                    color="0.55", ms=4, capsize=2, label="δ, single layer (E02; drift-sensitive)")
        ax.errorbar(d.ddelta, d.y + 0.015, xerr=[d.ddelta - d.ddelta_lo, d.ddelta_hi - d.ddelta], fmt="s",
                    color=COL[b], ms=4, capsize=2, label="Δδ vs deep reference, cruise FE (E06)")
        ax.axhspan(27.8, 28.6, color="0.93", zorder=0)
        ax.text(-0.33, 28.3, "reference\nlayers", fontsize=7, color="0.4", va="center")
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set(title=b, xlim=(-0.35, 0.4), xlabel="Excess of DIC over TRACE\n(µmol kg$^{-1}$ yr$^{-1}$)")
    axs[0].invert_yaxis()
    axs[0].set_ylabel("γ$^n$ layer midpoint (kg m$^{-3}$)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.1), frameon=False)
    save(fig, "fig7_drift_immune_excess")


def fig8_scale():
    e06 = pd.read_csv(RES / "e06_differential_trends.csv")
    sc = pd.read_csv(RES / "e11_scale_factor.csv").set_index("set")
    fig, ax = plt.subplots(figsize=(5.2, 4))
    x = np.linspace(0, 0.95, 10)
    m, lo, hi = sc.loc["all", ["m", "m_lo", "m_hi"]]
    ax.fill_between(x, lo * x, hi * x, color="0.88", label=f"all basins: m = {m:.2f} [{lo:.2f}, {hi:.2f}]")
    ax.plot(x, m * x, color="0.3", lw=1)
    for b in BASINS:
        s = e06[e06.basin == b]
        ax.errorbar(s.dbeta_cant, s.ddelta, yerr=[s.ddelta - s.ddelta_lo, s.ddelta_hi - s.ddelta], fmt="o",
                    color=COL[b], ms=4, capsize=2, label=f"{b}: m = {sc.loc[b, 'm']:.2f}")
    ax.axhline(0, color="0.6", lw=0.7)
    ax.set(xlabel="TRACE C$_{ant}$ rate relative to deep reference\n(µmol kg$^{-1}$ yr$^{-1}$)",
           ylabel="Drift-immune excess Δδ (µmol kg$^{-1}$ yr$^{-1}$)",
           title="Excess scales with the TRACE rate")
    ax.legend(fontsize=7, loc="upper left")
    save(fig, "fig8_scale_factor")


def fig9_trace_sensitivity():
    d = pd.read_csv(RES / "e07_trace_sensitivity.csv")
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    mk = {"default": ("o", "k"), "dg1.0": ("<", "C0"), "dg1.6": (">", "C0"), "diseq0.9": ("v", "C3"),
          "diseq1.1": ("^", "C3"), "sf0.8": ("d", "C2"), "sf1.2": ("D", "C2"), "theta": ("x", "C4")}
    lab = {"default": "default", "dg1.0": "Δ/Γ = 1.0", "dg1.6": "Δ/Γ = 1.6", "diseq0.9": "diseq × 0.9",
           "diseq1.1": "diseq × 1.1", "sf0.8": "mean age × 0.8", "sf1.2": "mean age × 1.2", "theta": "θ input"}
    for ax, b in zip(axs, BASINS):
        s = d[d.basin == b].assign(y=lambda q: q.layer.map(mid))
        for i, (v, (m, c)) in enumerate(mk.items()):
            q = s[s.variant == v]
            ax.plot(q.delta, q.y + (i - 3.5) * 0.008, m, color=c, ms=4, label=lab[v], mfc="none" if v != "default" else c)
        e12 = RES / "e12_trace_unadjusted.csv"
        if e12.exists():
            u = pd.read_csv(e12)
            u = u[u.basin == b].assign(y=lambda q: q.layer.map(mid))
            ax.plot(u.delta_atm, u.y + 0.04, "*", color="C1", ms=6, label="unadjusted atm. CO$_2$ (E12)")
        dd = s[s.variant == "default"]
        ax.hlines(dd.y, dd.delta_lo, dd.delta_hi, color="k", lw=0.8)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set(title=b, xlim=(-0.2, 0.32), xlabel="δ = β$_{DIC}$ − β$_{Cant}$\n(µmol kg$^{-1}$ yr$^{-1}$)")
    axs[0].invert_yaxis()
    axs[0].set_ylabel("γ$^n$ layer midpoint (kg m$^{-3}$)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.14), frameon=False)
    save(fig, "fig9_trace_sensitivity")


def fig10_emlr_products():
    p = RES / "e10_emlr_products.csv"
    if not p.exists():
        return
    d = pd.read_csv(p)
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, b in zip(axs, BASINS):
        s = d[d.basin == b].assign(y=lambda q: q.layer.map(mid)).sort_values("y")
        s = s[s.trace_m23_94_04 > 0.05]
        for key, m, lab in (("g19_94_07", "o", "Gruber 2019, 1994–2007"), ("m23_94_04", "s", "Müller 2023, 1994–2004"),
                            ("m23_04_14", "D", "Müller 2023, 2004–2014")):
            ax.plot(s[f"ratio_emlr_trace_{key}"], s.y, m + "-", ms=4, label=lab)
        ax.plot(s.R_obs, s.y, "k*", ms=7, label="observed DIC / TRACE (E02)")
        ax.axvspan(0.85, 1.15, color="0.9", zorder=0)
        ax.axvline(1, color="0.5", lw=0.8)
        ax.set(title=b, xlim=(0.4, 2.0), xlabel="Rate ratio to TRACE (same samples)")
    axs[0].invert_yaxis()
    axs[0].set_ylabel("γ$^n$ layer midpoint (kg m$^{-3}$)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.08), frameon=False)
    save(fig, "fig10_emlr_products")


def fig11_cstar():
    p = RES / "e08_cstar_controls.csv"
    if not p.exists():
        return
    d = pd.read_csv(p)
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, b in zip(axs, BASINS):
        s = d[d.basin == b].assign(y=lambda q: q.layer.map(mid))
        ax.errorbar(s.delta_aou, s.y - 0.012, xerr=[s.delta_aou - s.delta_aou_lo, s.delta_aou_hi - s.delta_aou],
                    fmt="o", color="k", ms=4, capsize=2, label="DIC | θ, S, AOU")
        ax.errorbar(s.delta_cstar, s.y + 0.012, xerr=[s.delta_cstar - s.delta_cstar_lo, s.delta_cstar_hi - s.delta_cstar],
                    fmt="s", color=COL[b], ms=4, capsize=2, label="C* | θ, S")
        ax.plot(-0.5 * s.beta_talk, s.y, "x", color="0.4", label="−½ × TAlk trend | θ, S")
        ax.axvline(0, color="0.5", lw=0.8)
        ax.set(title=b, xlim=(-0.3, 0.45), xlabel="Excess over TRACE / control\n(µmol kg$^{-1}$ yr$^{-1}$)")
    axs[0].invert_yaxis()
    axs[0].set_ylabel("γ$^n$ layer midpoint (kg m$^{-3}$)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.08), frameon=False)
    save(fig, "fig11_cstar_controls")


def fig12_perfect_model():
    """E13: rates in MRI-ESM2-0 sampled like GLODAP, beside the observed rates."""
    p = RES / "e13_model_layer_trends.csv"
    if not p.exists():
        return
    m = pd.read_csv(p)
    o = pd.read_csv(RES / "e02_layer_trends_primary.csv")
    fig, axs = plt.subplots(1, 4, figsize=(11, 3.8), sharey=True)
    for ax, b in zip(axs, BASINS):
        s = m[m.basin == b].assign(y=lambda q: q.layer.map(mid)).sort_values("y")
        t = o[o.basin == b].assign(y=lambda q: q.layer.map(mid)).sort_values("y")
        ax.plot(s.beta_canth_true, s.y, "-", color="k", lw=2, label="model true C$_{ant}$")
        ax.plot(s.beta_tco2, s.y, "o", color=COL[b], ms=4, label="model DIC | θ, S, AOU")
        ax.plot(s.beta_canth_trace, s.y, "s--", color="0.5", ms=3, mfc="none", label="TRACE on model T, S")
        pa = RES / "e13_model_layer_trends_atm.csv"
        if pa.exists():
            a = pd.read_csv(pa)
            a = a[a.basin == b].assign(y=lambda q: q.layer.map(mid)).sort_values("y")
            ax.plot(a.beta_canth_trace, a.y, "^:", color="C1", ms=3, label="TRACE on model T, S, unadjusted atm. CO$_2$ (E13d)")
        ax.plot(t.beta_tco2, t.y, "*", color="C3", ms=6, label="observed DIC | θ, S, AOU (E02)")
        ax.plot(t.beta_canth, t.y, "x", color="C3", ms=4, label="TRACE on observed T, S (E02)")
        ax.axvline(0, color="0.6", lw=0.7)
        ax.set(title=b, xlabel="Trend (µmol kg$^{-1}$ yr$^{-1}$)")
    axs[0].invert_yaxis()
    axs[0].set_ylabel("γ$^n$ layer midpoint (kg m$^{-3}$)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.16), frameon=False)
    save(fig, "fig12_perfect_model")


def fig13_models():
    """Upper-ocean (gamma <= 27.0) ratios by basin: observed DIC/TRACE vs model truth/TRACE (E13, E15)."""
    obs = pd.read_csv(RES / "e09_pooled.csv")
    series = [("observed DIC / TRACE (E09)", obs, None, "k", "*")]
    for f, lab, c, mk in (("e13_model_pooled_to2014.csv", "MRI-ESM2-0, true / TRACE", "C0", "o"),
                          ("e13_model_pooled_to2014_atm.csv", "MRI-ESM2-0, true / TRACE unadjusted", "C0", "^"),
                          ("e13_model_pooled_noresm2_lm.csv", "NorESM2-LM, true / TRACE", "C2", "s"),
                          ("e13_model_pooled_noresm2_lm_atm.csv", "NorESM2-LM, true / TRACE unadjusted", "C2", "v")):
        if (RES / f).exists():
            series.append((lab, pd.read_csv(RES / f), "b_true_vs_trace", c, mk))
    if len(series) < 5:
        return
    cats = ["all"] + BASINS
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    for k, (lab, d, test, c, mk) in enumerate(series):
        d = d[d.layer_set.str.startswith("upper")]
        if test:
            d = d[d.test == test]
        d = d.set_index("basin").reindex(cats)
        x = np.arange(len(cats)) + (k - 2) * 0.13
        ax.errorbar(x, d.R, yerr=[d.R - d.R_lo, d.R_hi - d.R], fmt=mk, color=c, ms=6 if mk == "*" else 4,
                    mfc="none" if lab.endswith("unadjusted") else c, capsize=2, label=lab)
    ax.axhline(1, color="0.5", lw=0.8)
    ax.set_xticks(np.arange(len(cats)), ["All basins"] + BASINS)
    ax.set(ylabel="Upper-ocean rate ratio (γ$^n$ ≤ 27.0)", title="TRACE shortfall: observations and two Earth system models")
    ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.12), frameon=False)
    save(fig, "fig13_models")


def main():
    cfg = load_config()
    try:
        fig1_map_coast(cfg)
    except ImportError:
        print("cartopy unavailable: fig1 left as is")
    fig7_differential()
    fig8_scale()
    fig9_trace_sensitivity()
    fig10_emlr_products()
    fig11_cstar()
    fig12_perfect_model()
    fig13_models()
    print("round-2 figures written to", OUT)


if __name__ == "__main__":
    main()
