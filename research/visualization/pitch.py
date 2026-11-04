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
    f.text(0.05, 0.20, "Ishaan Sodhi   ·   5 November 2026", fontsize=14, color=INK)
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
        "  Cell effects → compare the same place over time (removes where-ships-went bias)",
        "  Run with y = observed DIC  →  β_DIC;   with y = TRACE C_ant  →  β_Cant",
        "Uncertainty: cruise-level bootstrap (fixed the repo's broken bootstrap)",
        "Pre-registered hypotheses and pass/fail rules before running anything",
        "Validated: synthetic data with known trends → bias ≤ 0.012 µmol kg⁻¹ yr⁻¹",
    ], y=0.76, size=16, gap=0.08)
    slides.append(f)

    # 5 Result 1
    f = slide("Result 1 — a clear human-CO₂ signature", "Results")
    image(f, FIG / "fig2_trend_profiles.png", [0.03, 0.22, 0.94, 0.58])
    f.text(0.05, 0.13, "≈ 1 µmol kg⁻¹ yr⁻¹ at the surface in every basin (matches CO₂ equilibrium theory), decaying with depth;\n"
                       "reaches North Atlantic Deep Water (0.29) but not old Indian deep water (≈0).  TRACE within ±15 % in 18 of 20 upper layers.",
           fontsize=13, color=INK, va="top", linespacing=1.5)
    slides.append(f)

    # 6 Result 2
    f = slide("Result 2 — upper ocean: 10–25 % faster than TRACE", "Results")
    image(f, FIG / "fig4_delta_robustness.png", [0.03, 0.22, 0.94, 0.58])
    f.text(0.05, 0.13, "δ = observed − TRACE ≈ +0.10 to +0.17 µmol kg⁻¹ yr⁻¹ in light layers (all basins except Atlantic surface).\n"
                       "Same sign in all 11 sensitivity variants and every leave-one-cruise-out test; unchanged using post-1990 data only.",
           fontsize=13, color=INK, va="top", linespacing=1.5)
    slides.append(f)

    # 7 Result 3
    f = slide("Result 3 — eMLR(C*) sides with the data", "Results")
    image(f, FIG / "fig6_emlr_comparison.png", [0.03, 0.22, 0.94, 0.58])
    f.text(0.05, 0.13, "Gruber et al. 2019 eMLR(C*) also exceeds TRACE by 10–35 % in upper layers → 13 of 15 robust layers point to 'TRACE too low'.\n"
                       "Southern Ocean intermediate waters show the opposite sign.  Caveat: eMLR uses the same GLODAP DIC (not fully independent).",
           fontsize=13, color=INK, va="top", linespacing=1.5)
    slides.append(f)

    # 8 Confidence
    f = slide("How sure are we?", "Honest assessment")
    box(f, 0.05, 0.18, 0.42, 0.58,
        "✓ Estimator unbiased on synthetic data\n"
        "✓ 11 sensitivity variants: same sign\n"
        "✓ Leave-one-cruise-out: no sign flips\n"
        "✓ Post-1990 only: excess unchanged\n"
        "✓ Phosphate instead of AOU: same sign\n"
        "✓ QC choice barely matters (0/14 layers)\n"
        "✓ Pre-registered; all runs logged",
        title="What holds up", color=GOOD, size=15)
    box(f, 0.53, 0.18, 0.42, 0.58,
        "! Excess (2.5–4 µmol kg⁻¹ over the record) is\n   at the level of possible calibration drift\n"
        "! eMLR(C*) is not independent of the data\n"
        "! TRACE run only with default settings\n"
        "! One linear trend over 1972–2021\n"
        "! CIs cover ≈ 91 % (not 95 %) under\n   large cruise offsets\n"
        "! Literature search not yet systematic",
        title="What could break it", color=WARN, size=15)
    slides.append(f)

    # 9 Claim
    f = slide("The claim we think we can defend", "Takeaway")
    box(f, 0.05, 0.50, 0.9, 0.26,
        "Observed interior-ocean carbon and tracer-based anthropogenic-carbon rates agree to first order,\n"
        "but DIC-based estimates accumulate systematically faster than TRACE in the upper ocean outside the\n"
        "Southern Ocean — a robust DIC-versus-tracer discrepancy that is not yet distinguishable from\n"
        "measurement drift or TRACE parameter uncertainty.",
        color=ACCENT, size=16)
    bullets(f, [
        "Side result (methods note): static DIC regressions must model time / C_ant;",
        "  temporal 'validation' of empirical carbon models is otherwise confounded",
        "What we are NOT claiming: that TRACE is wrong, or that natural carbon is increasing",
    ], y=0.40, size=15, gap=0.08)
    slides.append(f)

    # 10 Next steps / asks
    f = slide("Next steps — and where we'd value your input", "Plan")
    box(f, 0.05, 0.14, 0.42, 0.62,
        "1. TRACE sensitivity: disequilibrium and\n    TTD-shape parameters (days)\n"
        "2. Independent comparator: Müller 2023\n    field + ocean-model C_ant (1–2 weeks)\n"
        "3. Drift tests by era / CRM use (days)\n"
        "4. Systematic literature search\n"
        "5. Target: Biogeosciences / Ocean Science\n    (GRL if the excess survives 1–3);\n"
        "    OSM / EGU abstract now",
        title="Next 2–4 weeks", size=15)
    box(f, 0.53, 0.14, 0.42, 0.62,
        "• Is the upper-ocean excess plausible to you,\n  or a known TTD/TRACE issue?\n"
        "• Best way to rule out DIC calibration drift?\n"
        "• Access to ocean-model C_ant output\n  (e.g., RECCAP2) for an independent test?\n"
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
