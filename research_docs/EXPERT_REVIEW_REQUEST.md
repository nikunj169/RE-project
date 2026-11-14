# Expert review request (draft, not sent)

*Prepared 2026-11-14 for FINAL_REVIEW §5 item 5. The author decides whom to contact and sends it.*

## Suggested reviewers (by role)
* **A TRACE developer** (corresponding authors of Carter et al., 2025, or Sandborn et al., 2026). They can judge the boundary-condition finding (E12, E13d) and whether the "adjusted" CO₂ trajectory is intended to be damped.
* **An eMLR(C*) developer** (Gruber et al., 2019, or Müller et al., 2023). They can judge the comparison at the same samples (E10) and the deep Atlantic deficit.
* **A GLODAP data-quality specialist.** They can judge the cruise-offset and drift arguments (E06, E16) and the TAlk drift finding (E08).

## Draft e-mail

> Subject: Request for informal comments: sample-matched comparison of GLODAP DIC trends with TRACE C_ant
>
> Dear Dr. [name],
>
> I am writing to ask whether you would be willing to give informal comments on a manuscript before journal submission. We apply one within-location trend estimator to observed GLODAPv2.2023 DIC and to TRACE C_ant at the same 342,277 samples (1972–2021). To remove calibration drift, a differential version with cruise fixed effects is immune to any additive per-cruise offset.
>
> The main result is that upper-ocean DIC accumulates 16–18 % faster than TRACE predicts:
> * the excess is similar in all four basins;
> * two eMLR(C*) products at the same samples show a similar excess;
> * an Earth system model sampled like GLODAP shows the same shortfall.
>
> Replacing TRACE's adjusted surface CO₂ trajectory with the atmospheric record removes about 60 % of the excess in observations and most of it in the model. Your work on [TRACE / eMLR(C*) / GLODAP] makes your view on [the boundary condition / the comparison / the drift controls] especially valuable.
>
> The manuscript (about 18 pages), supplement and code are attached or linked. Even brief comments would be much appreciated.
>
> With best regards,
> Ishaan Sodhi

## Attach
* `paper/submission_gbc/main.pdf` (journal-format draft)
* `paper/SUPPLEMENTARY_MATERIAL.md` (or a PDF of it)
* The repository link, or a Zenodo DOI once released

## Remaining steps for the author before submission
1. **Expert review.** Send the request above.
2. **Front matter.** Add affiliation and e-mail in `paper/submission_gbc/template.tex`.
3. **Licence and Zenodo.**
   * Choose a licence (none is set) and add `LICENSE`.
   * Enable the repository in Zenodo's GitHub integration and create a GitHub release. `.zenodo.json` and `CITATION.cff` supply the metadata.
   * Add the DOI to the Open Research section.
4. **AGU class.** Move `main.tex` content into the AGU class from AGU's Overleaf template. This is not on CTAN, so the draft uses `article` with the same sections.
