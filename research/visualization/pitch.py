"""Build a short slide-deck PDF pitch for an advisor meeting.

Usage: python -m research.visualization.pitch
Output: paper/pitch_for_professor.pdf (+ PNG previews in research/figures/pitch/)
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch

from research.data_io import ROOT

FIG = ROOT / "research" / "figures"
OUT = ROOT / "paper" / "pitch_for_professor.pdf"
PREV = FIG / "pitch"
W, H = 13.333, 7.5
INK, MUTED, ACCENT, GOOD, WARN = "#1b2a3a", "#5b6b7b", "#1f6f9f", "#2e7d32", "#b5651d"
plt.rcParams.update({"font.family": "DejaVu Sans"})


def slide(title: str | None = None, kicker: str | None = None):
    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor("white")
    if kicker:
        fig.text(0.05, 0.925, kicker.upper(), fontsize=11, color=ACCENT, weight="bold")
    if title:
        fig.text(0.05, 0.86, title, fontsize=25, color=INK, weight="bold", va="center")
    return fig


def bullets(fig, items, x=0.05, y=0.74, size=16, gap=0.075, width=None):
    for it in items:
        sub = it.startswith("  ")
        txt = it.strip()
        fig.text(x + (0.03 if sub else 0), y, ("–  " if sub else "•  ") + txt,
                 fontsize=size - (2 if sub else 0), color=MUTED if sub else INK, va="top", wrap=True)
        y -= gap * (0.8 if sub else 1)
    return y


def box(fig, x, y, w, h, text, color=ACCENT, size=15, title=None):
    ax = fig.add_axes([x, y, w, h])
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0.0,rounding_size=0.04",
                                transform=ax.transAxes, facecolor=color, alpha=0.09, edgecolor=color, lw=1.5))
    if title:
        ax.text(0.04, 0.88, title, fontsize=size, weight="bold", color=color, va="top", transform=ax.transAxes)
    ax.text(0.04, 0.88 - (0.2 if title else 0), text, fontsize=size - 1, color=INK, va="top",
            transform=ax.transAxes, wrap=True, linespacing=1.45)


def image(fig, path, rect):
    ax = fig.add_axes(rect)
    ax.imshow(mpimg.imread(path))
    ax.axis("off")


def footer(fig, n):
    fig.text(0.05, 0.035, "GLODAPv2.2023 · TRACE · eMLR(C*) — draft results, not yet peer-reviewed",
             fontsize=9, color=MUTED)
    fig.text(0.95, 0.035, str(n), fontsize=10, color=MUTED, ha="right")


def build():
    slides = []

    # 1 Title
    f = slide()
    f.text(0.05, 0.70, "Is the ocean storing human CO₂\nas fast as tracer ages predict?", fontsize=36,
           weight="bold", color=INK, va="top", linespacing=1.2)
    f.text(0.05, 0.43, "A sample-matched test of 50 years of GLODAP carbon data\nagainst tracer-based anthropogenic-carbon estimates",
           fontsize=18, color=MUTED, va="top", linespacing=1.4)
    f.text(0.05, 0.20, "Ishaan Sodhi   ·   14 November 2026 (updated with Rounds 2–3)", fontsize=14, color=INK)
    f.add_artist(plt.Line2D([0.05, 0.35], [0.27, 0.27], color=ACCENT, lw=3))
    slides.append(f)

    # 2 Where it started
    f = slide("It started as an audit of our TCO₂ regression project", "Background")
    bullets(f, [
        "Original claim: DIC ≈ a 5-parameter 'rank-1 quadratic' in S, T, AOU — 'low-dimensional'",
        "Re-analysis found every model under-predicts later years: +4 to +11 µmol kg⁻¹ bias",
        "  The bias grows with time and is largest near the surface — the shape of human CO₂ uptake",
        "Adding a simple time term removes 77–83 % of that bias (median)",
        "Which model 'wins' flips with the test design; rank-1 vs rank-2 vs quadratic",
        "  are statistically equivalent within ±4 µmol kg⁻¹ (DIC measurement precision)",
        "15 % of the data used were calculated or non-QC'd DIC (GLODAP flags ignored)",
    ], y=0.76, size=16, gap=0.085)
    box(f, 0.05, 0.08, 0.9, 0.12, "So the 'error' was a signal. Instead of hiding it, we made it the research question.",
        color=ACCENT, size=16)
    slides.append(f)

    # 3 Question
    f = slide("Research question", "The idea")
    f.text(0.05, 0.70, "Within the same places and at the same water properties (θ, S, AOU),\n"
                       "does ocean DIC rise at the rate that tracer-based models of\n"
                       "anthropogenic carbon (TRACE) predict — layer by layer, basin by basin?",
           fontsize=21, color=INK, va="top", linespacing=1.5, style="italic")
    box(f, 0.05, 0.10, 0.42, 0.32,
        "• DIC-based methods (eMLR) and tracer-based\n  methods (TTD/TRACE) both estimate human CO₂\n"
        "• They are usually compared as maps or inventories\n"
        "• Never with one identical estimator at the\n  same samples over the whole record (that we found)",
        title="Why it matters", size=14)
    box(f, 0.53, 0.10, 0.42, 0.32,
        "• Agreement → validates TRACE's rates\n  from carbon data directly\n"
        "• Disagreement → either TRACE is biased,\n  or natural carbon / circulation is changing\n"
        "• Both outcomes are publishable",
        title="What we learn", color=GOOD, size=14)
    slides.append(f)

    # 4 Method
    f = slide("Design: one estimator, two inputs, same samples", "Method")
    bullets(f, [
        "Data: 342,277 strict-QC GLODAPv2.2023 samples, 515 cruises, 1972–2021",
        "Grid: 4°×4° cells revisited over ≥ 10 years; 9 neutral-density layers; 4 basins",
        "Model per basin × layer:  y = cell effect + f(θ, S, AOU) + β·year",
        "  Run with y = observed DIC  →  β_DIC;   with y = TRACE C_ant  →  β_Cant",
        "New: drift-immune version — upper layer vs deep layers within the same cruises",
        "  Cruise fixed effects cancel any per-cruise calibration offset, even one drifting in time",
        "Uncertainty: cruise-level bootstrap shared across all groups (valid pooled CIs)",
        "Pre-registered pass/fail rules; validated on synthetic data with drift (bias ≤ 0.002)",
    ], y=0.76, size=16, gap=0.08)
    slides.append(f)

    # 5 Result 1
    f = slide("Result 1 — a clear human-CO₂ signature", "Results")
    image(f, FIG / "fig2_trend_profiles.png", [0.03, 0.22, 0.94, 0.58])
    f.text(0.05, 0.13, "≈ 1 µmol kg⁻¹ yr⁻¹ at the surface in every basin (matches CO₂ equilibrium theory), decaying with depth;\n"
                       "reaches North Atlantic Deep Water but not old Indian deep water.",
           fontsize=13, color=INK, va="top", linespacing=1.5)
    slides.append(f)

    # 6 Result 2
    f = slide("Result 2 — upper ocean 16–18 % faster than TRACE, not drift", "Results")
    image(f, FIG / "fig7_drift_immune_excess.png", [0.03, 0.22, 0.94, 0.58])
    f.text(0.05, 0.13, "Pooled upper ocean (γⁿ ≤ 27.0): observed / TRACE = 1.16 [1.11, 1.21].  Drift-immune (cruise offsets removed): 1.18 [1.11, 1.23].\n"
                       "Similar in all four basins (16–20 %); 10 of 12 upper layers survive the drift test; robust to the choice of deep reference (1.16–1.19).",
           fontsize=13, color=INK, va="top", linespacing=1.5)
    slides.append(f)

    # 7 Result 3
    f = slide("Result 3 — two eMLR(C*) products agree", "Results")
    image(f, FIG / "fig10_emlr_products.png", [0.03, 0.22, 0.94, 0.58])
    f.text(0.05, 0.13, "At the same samples: Gruber 2019 (1994–2007) 1.13;  Müller 2023 (1994–2004) 1.20 and (2004–2014) 1.13 × TRACE.\n"
                       "Even their lowest sensitivity cases exceed TRACE.  Caveat: both use GLODAP DIC, so not independent.",
           fontsize=13, color=INK, va="top", linespacing=1.5)
    slides.append(f)

    # 8 Result 4 (model)
    f = slide("Result 4 — an independent check in an Earth system model", "Results")
    image(f, FIG / "fig12_perfect_model.png", [0.03, 0.22, 0.94, 0.58])
    f.text(0.05, 0.13, "MRI-ESM2-0 sampled at every GLODAP sample (true C_ant known exactly; no GLODAP DIC involved):\n"
                       "our estimator is unbiased (factor −0.01), and TRACE falls 12 % short of the model's true upper-ocean uptake.",
           fontsize=13, color=INK, va="top", linespacing=1.5)
    slides.append(f)

    # 9 Mechanism
    f = slide("Why? Most of it is TRACE's surface CO₂ boundary condition", "Mechanism")
    bullets(f, [
        "TRACE's built-in surface-ocean CO₂ history rises at only 0.88 × the atmospheric rate",
        "Swap in the real atmospheric record:",
        "  Observations: excess falls from 16 % to 6 % [1, 10] → explains ≈ 60 %",
        "  Earth system model: shortfall falls from 12 % to 2 % → explains ≈ 85–90 %",
        "Time series (BATS, HOT, ESTOC) show surface pCO₂ keeping pace with the atmosphere",
        "Other TRACE knobs (age-distribution shape, ages ±20 %, θ input) matter little",
        "Deep Atlantic: 25–35 % below TRACE in data and eMLR — slowing deep ventilation",
    ], y=0.76, size=16, gap=0.085)
    box(f, 0.05, 0.08, 0.9, 0.12, "The mechanism was found after seeing the data (labelled post hoc); the model test was logged before running.",
        color=WARN, size=14)
    slides.append(f)

    # 10 Confidence
    f = slide("How sure are we?", "Honest assessment")
    box(f, 0.05, 0.14, 0.42, 0.62,
        "✓ Survives cruise-offset/drift removal\n"
        "✓ Same in all four basins; pooled CIs\n"
        "✓ 7 TRACE structural variants\n"
        "✓ C* (phosphate + alkalinity): 9/12 layers\n"
        "✓ Two eMLR(C*) products agree\n"
        "✓ ESM: estimator unbiased; TRACE short\n"
        "✓ Robust to deep reference choice\n"
        "✓ Pre-registered; all runs logged",
        title="What holds up", color=GOOD, size=15)
    box(f, 0.53, 0.14, 0.42, 0.62,
        "! All observational evidence uses GLODAP DIC\n"
        "! One Earth system model so far; its\n   basin pattern differs from the data\n"
        "! Within-cruise, depth-dependent\n   artefacts not excluded\n"
        "! Pacific thermocline depends on AOU vs\n   phosphate correction\n"
        "! Not yet reviewed by a carbon expert",
        title="What could break it", color=WARN, size=15)
    slides.append(f)

    # 11 Claim
    f = slide("The claim we think we can defend", "Takeaway")
    box(f, 0.05, 0.50, 0.9, 0.26,
        "Upper-ocean anthropogenic carbon accumulates 16–18 % faster than TRACE predicts at the\n"
        "same samples. The excess is not caused by cruise-level calibration drift, is matched by two\n"
        "eMLR(C*) products, and is reproduced (12 %) in an Earth system model. Most of it traces to\n"
        "TRACE's damped surface CO₂ boundary condition: an observational constraint on C_ant rates.",
        color=ACCENT, size=16)
    bullets(f, [
        "Side results: a reusable drift-immune diagnostic for any interior C_ant product;",
        "  alkalinity drift at depth (matters for C* methods); static DIC regressions must model time",
        "What we are NOT claiming: that the TTD approach is wrong, or that natural carbon is increasing",
    ], y=0.40, size=15, gap=0.08)
    slides.append(f)

    # 12 Next steps / asks
    f = slide("Next steps — and where we'd value your input", "Plan")
    box(f, 0.05, 0.14, 0.42, 0.62,
        "1. Second Earth system model (NorESM2-LM):\n    code ready, needs ~45 GB download\n"
        "2. Contact TRACE / eMLR(C*) developers\n"
        "3. Convert to journal template; GBC or\n    Biogeosciences (GRL if model 2 agrees)\n"
        "4. OSM 2028 / EGU 2027 abstract now\n"
        "5. Zenodo release of code + results",
        title="Next 2–4 weeks", size=15)
    box(f, 0.53, 0.14, 0.42, 0.62,
        "• Is the boundary-condition explanation\n  plausible, or a known TRACE issue?\n"
        "• Best way to rule out within-cruise,\n  depth-dependent DIC artefacts?\n"
        "• Who should we show this to first?\n"
        "• Should the original rank-1 paper be\n  retired or reframed?\n"
        "• Venue and timeline you'd recommend",
        title="Questions for you", color=GOOD, size=15)
    slides.append(f)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    PREV.mkdir(parents=True, exist_ok=True)
    with PdfPages(OUT) as pdf:
        for i, f in enumerate(slides, 1):
            if i > 1:
                footer(f, i)
            pdf.savefig(f)
            f.savefig(PREV / f"slide{i:02d}.png", dpi=60)
            plt.close(f)
    print(OUT, len(slides), "slides")


if __name__ == "__main__":
    build()
