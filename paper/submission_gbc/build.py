"""Build the journal-format (AGU / Global Biogeochemical Cycles layout) submission draft.

Reads paper/MANUSCRIPT.md, paper/FIGURE_CAPTIONS.md, paper/REFERENCES.bib and the front matter in
this folder; nothing in paper/ is modified. Steps:
  * in-text author-year citations -> natbib keys (\\citep / \\citet via pandoc --natbib);
  * figures renumbered by order of first mention in the main text; Round-1 Figs. 3, 4 and 6 move
    to the supplement as Figs. S2-S4 (FINAL_REVIEW §5 item 6);
  * each main-text figure is placed after the paragraph that first cites it;
  * Data/Code availability become the AGU "Open Research" section.
Output: main.tex, figures/, figure_map.csv; compile with `tectonic main.tex`.

Usage: python paper/submission_gbc/build.py
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
FIGS = PAPER.parent / "research" / "figures"
SUPP_OLD = ["S1", "3", "4", "6"]  # Round-1 figures moved to the supplement, in this order -> S1..S4

CITE = {  # (surname as written, year) -> bib key
    ("Sabine", "2004"): "Sabine2004", ("Sabine", "2008"): "Sabine2008", ("Gruber", "2019"): "Gruber2019",
    ("Gruber", "1996"): "Gruber1996", ("Müller", "2023"): "Mueller2023", ("Friis", "2005"): "Friis2005",
    ("Clement & Gruber", "2018"): "ClementGruber2018", ("Carter", "2017"): "Carter2017",
    ("Carter", "2025"): "Carter2025", ("Steinfeldt", "2024"): "Steinfeldt2024", ("Sandborn", "2026"): "Sandborn2026",
    ("Olsen", "2016"): "Olsen2016", ("Lauvset", "2024"): "Lauvset2024", ("Bates", "2014"): "Bates2014",
    ("Eyring", "2016"): "Eyring2016", ("Yukimoto", "2019"): "Yukimoto2019", ("Seland", "2020"): "Seland2020",
}
NAME = r"(Clement & Gruber|[A-ZÄÖÜ][a-zäöü]+)(?: et al\.)?"


def key(name, year):
    return "@" + CITE[(name, year)]


def convert_citations(text: str) -> str:
    # narrative: "Müller et al. (2023)" -> @Mueller2023
    text = re.sub(NAME + r" \((\d{4})\)", lambda m: key(m[1], m[2]), text)

    # parenthetical groups: "(e.g., A et al., 2004; B, 2019, Table 2)" -> "[e.g., @A; @B, Table 2]"
    def group(m):
        inner = re.sub(NAME + r", (\d{4})", lambda n: key(n[1], n[2]), m[1])
        return "[" + inner + "]" if "@" in inner else m[0]
    return re.sub(r"\(([^()]*\d{4}[^()]*)\)", group, text)


def outside_code(text: str, fn) -> str:
    parts = text.split("`")
    return "`".join(fn(p) if i % 2 == 0 else p for i, p in enumerate(parts))


SUP = dict(zip("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789−+"))
SUB = dict(zip("₀₁₂₃₄₅₆₇₈₉", "0123456789"))


def symbols(t: str) -> str:
    t = t.replace("C*", r"C\*")
    t = re.sub("[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+", lambda m: "^" + "".join(SUP[c] for c in m[0]) + "^", t)  # font lacks these
    t = re.sub("[₀₁₂₃₄₅₆₇₈₉]+", lambda m: "~" + "".join(SUB[c] for c in m[0]) + "~", t)
    t = re.sub(r"(?<=[A-Za-zβδγΔ])_([A-Za-z0-9]+)", r"~\1~", t)  # C_ant -> C~ant~
    return re.sub(r"(?<=~)\^([A-Za-z]+)", r"^\1^", t)  # C~ant~^true -> superscript


def captions() -> dict[str, tuple[str, str]]:
    """old number -> (file stem, one-paragraph caption)."""
    out = {}
    blocks = re.split(r"\n(?=\*\*Figure )", (PAPER / "FIGURE_CAPTIONS.md").read_text())
    for b in blocks[1:]:
        m = re.match(r"\*\*Figure (S?\d+) \(`([^`]+)`\)\.\*\*\s*(.*)", b, re.S)
        lines = [l.strip() for l in m[3].strip().splitlines() if l.strip() and not l.startswith("*Purpose")]
        lines = [l[2:] if l.startswith("- ") else l for l in lines]
        cap = " ".join(l if l[-1] in ".:;)" else l + "." for l in lines)
        out[m[1]] = (m[2], cap)
    return out


def main() -> None:
    md = (PAPER / "MANUSCRIPT.md").read_text()
    body = md[md.index("## 1. Introduction"):md.index("## References")]
    body = re.sub(r"## Data availability\n", "## Open Research {-}\n\n**Data.**\n\n", body)
    body = body.replace("## Code availability\n", "**Software.** ")
    body = re.sub(r"^(#+) \d+(\.\d+)*\.? ", r"\1 ", body, flags=re.M)  # LaTeX numbers sections
    lines = body.split("\n")  # pandoc needs a blank line before a list that follows a paragraph
    body = "\n".join(("\n" + l) if re.match(r"\s*(- |\d+\. )", l) and i and lines[i - 1].strip()
                      and not re.match(r"\s*(- |\d+\. )", lines[i - 1]) and not lines[i - 1].startswith("|")
                      else l for i, l in enumerate(lines))

    caps = captions()
    caps = {k: v for k, v in caps.items() if (FIGS / f"{v[0]}.pdf").exists()}
    # numbering: main-text figures by first mention, Round-1 3/4/6 to the supplement
    mention = []
    for m in re.finditer(r"Figs?\. (S?\d+)", body):
        if m[1] not in mention:
            mention.append(m[1])
    main_old = [n for n in mention if n not in SUPP_OLD and n in caps]
    main_old += [n for n in caps if n not in main_old and n not in SUPP_OLD]  # uncited: appended
    new = {o: str(i + 1) for i, o in enumerate(main_old)}
    new.update({o: f"S{i + 1}" for i, o in enumerate(SUPP_OLD)})
    body = re.sub(r"(Figs?\. )(S?\d+)", lambda m: m[1] + "\0" + new.get(m[2], m[2]), body).replace("\0", "")

    # place each main figure after the paragraph that first cites it
    (HERE / "figures").mkdir(exist_ok=True)
    paras = body.split("\n\n")
    placed = set()
    out = []
    for p in paras:
        out.append(p)
        for o in main_old:
            if o not in placed and re.search(rf"Figs?\. {new[o]}\b", p):
                stem, cap = caps[o]
                shutil.copy(FIGS / f"{stem}.pdf", HERE / "figures" / f"{stem}.pdf")
                out.append(f"![{cap}](figures/{stem}.pdf){{width=100%}}")
                placed.add(o)
    for o in main_old:  # uncited figures go at the end of the results
        if o not in placed:
            stem, cap = caps[o]
            shutil.copy(FIGS / f"{stem}.pdf", HERE / "figures" / f"{stem}.pdf")
            out.insert(len(out) - 1, f"![{cap}](figures/{stem}.pdf){{width=100%}}")
    body = "\n\n".join(out)

    paths = re.findall(r"\]\(figures/[^)]*\)", body)  # protect image paths from symbols()
    for i, pth in enumerate(paths):
        body = body.replace(pth, f"](FIGPATH{i})", 1)
    body = outside_code(convert_citations(body), symbols)
    for i, pth in enumerate(paths):
        body = body.replace(f"](FIGPATH{i})", pth, 1)
    # captions carry citations and symbols too (already inside the image text)
    front = (HERE / "frontmatter.md").read_text()
    src = front + "\n\n" + body
    (HERE / "_build.md").write_text(src)
    subprocess.run(["pandoc", str(HERE / "_build.md"), "-f", "markdown+subscript+superscript",
                    "-t", "latex", "--natbib", "--shift-heading-level-by=-1", "-o", str(HERE / "_body.tex")], check=True)
    tex = (HERE / "_body.tex").read_text()
    tex = tex.replace(r"\pandocbounded{", "{")
    template = (HERE / "template.tex").read_text()
    (HERE / "main.tex").write_text(template.replace("%%BODY%%", tex))
    bib = re.sub(r"(@\w+\{[^,]+,)\s*%[^\n]*", r"\1", (PAPER / "REFERENCES.bib").read_text())  # BibTeX: no comments in entries
    bib = re.sub(r",?\s*note = \{\[V[^\n]*\]\}(?=\})", "", bib)  # internal verification tags
    bib = re.sub(r",?\s*note = \{File [^}]*\}(?=\})", "", bib)  # checksums belong in the repo, not the paper
    (HERE / "references.bib").write_text(bib)
    rows = ["old,new,file"] + [f"{o},{new[o]},{caps[o][0]}" for o in sorted(new, key=lambda o: new[o]) if o in caps]
    (HERE / "figure_map.csv").write_text("\n".join(rows) + "\n")
    for f in ("_build.md", "_body.tex"):
        (HERE / f).unlink()
    print("figure map:", {o: new[o] for o in new if o in caps})


if __name__ == "__main__":
    main()
