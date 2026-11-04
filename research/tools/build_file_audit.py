"""Build research_docs/FILE_AUDIT.csv: one row per tracked repository file.

Annotations are hand-written from direct inspection (see REPOSITORY_INTELLIGENCE.md).
Rules are matched in order; the first matching glob wins. ``inspection_status`` values:
  read_full     - entire file read
  read_partial  - docstring, structure and I/O read; logic not fully audited
  profiled      - data file: schema/size/summary statistics inspected
  listed_only   - generated artifact (figure/PDF/TeX table/binary); existence recorded only
"""

from __future__ import annotations

import csv
import fnmatch
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "research_docs" / "FILE_AUDIT.csv"

# (glob, file_type, purpose, inputs, outputs, dependencies, scientific_relevance, status, issues)
RULES: list[tuple[str, ...]] = [
    ("README.md", "doc", "G3 study overview and quick-start", "-", "-", "-", "high: states claims", "read_full",
     "Lists src/visualization modules that do not exist; claims no hard-coded results (false, see scripts/10)"),
    (".gitignore", "config", "git ignore rules", "-", "-", "-", "none", "read_full",
     "Ignores *.csv/*.png/data/ yet these were force-committed"),
    ("Screenshot*.png", "image", "unreferenced screenshot", "-", "-", "-", "none", "listed_only", "Orphan file"),
    ("config/global_config.yaml", "config", "paths, variable names, QC thresholds, basins, splits, seeds, bootstrap, TOST margins", "-", "-", "-", "high", "read_full",
     "No GLODAP QC-flag filtering configured"),
    ("config/model_families.yaml", "config", "model family metadata", "-", "-", "-", "medium", "read_full", "Descriptive only; not loaded by scripts"),
    ("config/pysr_config.yaml", "config", "PySR operator sets A-D and common settings", "-", "-", "-", "high", "read_full",
     "Settings silently capped in experiments/run_symbolic_search.py; subsample key mismatch (subsample.max_train_samples vs common.subsample_max)"),
    ("config/selected_symbolic_candidate.yaml", "config", "frozen per-basin SR equations (written by scripts/09)", "symbolic_candidates_all.csv", "-", "-", "high", "read_full",
     "Selected by min 2015-17 RMSE over 599 candidates -> optimistic validation score"),
    ("config/watermass.yaml", "config", "threshold water-mass definitions", "-", "-", "-", "medium", "read_full",
     "In-situ T, depth-only; NADW not restricted to Atlantic; duplicated in src/analysis/watermass.py"),
    ("data/processed/glodap_qc_filtered.csv", "data", "QC-filtered GLODAP subset used by every G3 script", "GLODAPv2.2023 .mat (absent)", "-", "-", "critical", "profiled",
     "512,384 rows; only S,T(in-situ),AOU,TCO2 + metadata; no flags/pressure/theta/O2/nutrients/TAlk"),
    ("data/processed/basin_counts.csv", "data", "per-basin/per-split counts", "scripts/01", "-", "-", "low", "read_full", ""),
    ("data/processed/qc_stats.json", "data", "QC filter funnel counts", "scripts/01", "-", "-", "low", "read_full", ""),
    ("environment/environment.yml", "config", "conda env spec", "-", "-", "-", "low", "read_full", "Pins python<3.14 but results produced on 3.14.3"),
    ("environment/requirements.txt", "config", "pip requirements", "-", "-", "-", "low", "read_full", "Unpinned lower bounds only"),
    ("environment/julia/VERSIONS.md", "doc", "verified tool versions", "-", "-", "-", "low", "read_full", "PySR 2.3.0 used; 2.7.2 now installed"),
    ("src/__init__.py", "python", "package marker", "-", "-", "-", "none", "read_full", ""),
    ("src/*/__init__.py", "python", "package marker", "-", "-", "-", "none", "read_full", ""),
    ("src/config_loader.py", "python", "YAML config loading and path resolution", "config/*.yaml", "dict/Path", "pyyaml", "low", "read_full", ""),
    ("src/data_loader.py", "python", "load GLODAP .mat (scipy/h5py), range QC, build DataFrame", ".mat", "DataFrame, QC stats", "scipy, h5py, config_loader", "critical", "read_full",
     "No QC flags; open-interval range filters only; pressure/theta/O2/nutrients/TAlk not extracted"),
    ("src/basins.py", "python", "basin assignment by GLODAP region code + 35S cut", "DataFrame", "Series", "numpy, pandas", "high", "read_full",
     "Arctic (region 4) -> Unclassified and dropped"),
    ("src/splits.py", "python", "temporal split, spatial blocks, cruise groups, leakage assertions", "DataFrame", "masks/groups", "sklearn", "high", "read_full", "get_grouped_kfold_splits unused"),
    ("src/models/base.py", "python", "abstract model interface fit(S,T,A,Y)/predict", "-", "-", "abc", "medium", "read_full", ""),
    ("src/models/mean_baseline.py", "python", "Model 0: training mean", "S,T,A,Y", "pred", "base", "low", "read_full", ""),
    ("src/models/linear.py", "python", "Model 1: OLS on S,T,AOU/100", "S,T,A,Y", "pred", "sklearn", "medium", "read_full", ""),
    ("src/models/rank1_quadratic.py", "python", "Model 2: (aS+bT+gA/100+d)^2+e via curve_fit", "S,T,A,Y", "pred, Hessian eigenvalue", "scipy", "critical", "read_full",
     "Special case of rank-2 and full quadratic (nested) - not acknowledged in inference"),
    ("src/models/rank2_quadratic.py", "python", "Model 3: z1^2+z2^2+e via Nelder-Mead (5 restarts)", "S,T,A,Y", "pred", "scipy, rank1", "high", "read_full",
     "Nelder-Mead on 9 params is weak; training RMSE checked <= rank-1 so convergence adequate"),
    ("src/models/full_quadratic.py", "python", "Model 4: OLS full 2nd-order polynomial (raw AOU)", "S,T,A,Y", "pred", "sklearn", "high", "read_full",
     "Raw AOU^2 ~1e5 -> ill-conditioned; small numeric discrepancies vs scaled fit"),
    ("src/models/cubic.py", "python", "Model 5: OLS full cubic (raw AOU)", "S,T,A,Y", "pred", "sklearn", "medium", "read_full", "Ill-conditioned (AOU^3 ~3e7)"),
    ("src/models/symbolic.py", "python", "PySR wrapper; extracts hall-of-fame", "S,T,A,Y", "candidate DataFrame", "pysr", "high", "read_full", ""),
    ("src/models/frozen_symbolic.py", "python", "evaluate frozen sympy expression without refit", "equation", "pred", "sympy", "medium", "read_full",
     "Used in blocked CV with constants fit on test blocks (leakage)"),
    ("src/analysis/metrics.py", "python", "RMSE/MAE/R2/bias/within-X metrics", "y, yhat", "dict", "sklearn", "medium", "read_full", ""),
    ("src/analysis/pareto.py", "python", "Pareto frontier + legacy family classifier", "candidates", "frontier", "pandas", "medium", "read_full",
     "classify_equation_family looks for x0/x1/x2 but equations use sal/tmp/aou100 -> all 'constant'"),
    ("src/analysis/family_classification.py", "python", "string-heuristic SR family classifier", "equation str", "label", "re", "high", "read_full",
     "Not a test of algebraic rank; rational forms labelled linear; drives the 24/48 recurrence claim"),
    ("src/analysis/depth_analysis.py", "python", "metrics per depth bin", "y, yhat, depth", "DataFrame", "metrics", "medium", "read_full", ""),
    ("src/analysis/watermass.py", "python", "threshold water-mass labels + metrics", "DataFrame", "labels", "metrics", "medium", "read_full", "See config/watermass.yaml issues"),
    ("src/analysis/derivatives.py", "python", "analytic partial derivatives of rank-1", "coeffs, S,T,A", "summary", "numpy", "low", "read_full", ""),
    ("src/analysis/southern_ocean.py", "python", "SO abyss global vs local rank-1 refit", "SO arrays", "dict", "scipy", "medium", "read_full", ""),
    ("src/validation/bootstrap.py", "python", "observation and cluster bootstrap prediction intervals", "model, data, groups", "PI bounds", "numpy", "high", "read_full",
     "BUG: np.isin(groups, sampled_groups) drops duplicate cluster draws -> not a bootstrap"),
    ("src/validation/equivalence.py", "python", "per-observation TOST; bootstrap CI of paired AE difference", "y, preds, groups", "dict", "scipy", "high", "read_full",
     "Same np.isin bug; naive TOST ignores dependence"),
    ("src/validation/spatial.py", "python", "GroupKFold spatial/cruise CV helpers", "DataFrame", "fold metrics", "sklearn", "medium", "read_full", "Superseded by scripts/09"),
    ("scripts/01_preprocess.py", "python", "raw .mat -> processed CSV + basin counts + QC stats", ".mat", "data/processed/*", "src", "critical", "read_full", "Cannot run: raw data absent"),
    ("scripts/02_baseline_reproduction.py", "python", "rank-1 and full-quad per basin vs old paper values", "processed CSV", "baseline_*.csv/md", "src", "high", "read_full",
     "Old-paper reference values hard-coded (for discrepancy logging only)"),
    ("scripts/03_model_families.py", "python", "fit model hierarchy 0-5 per basin; val/ext metrics", "processed CSV", "model_family_*.csv", "src", "critical", "read_full", ""),
    ("scripts/04_pysr_pareto_search.py", "python", "older single-seed PySR search", "processed CSV", "symbolic_candidates_all.csv", "src, pysr", "low", "read_full",
     "Superseded; broken family classifier; would overwrite results"),
    ("scripts/05_equivalence_analysis.py", "python", "TOST + cluster bootstrap rank-1 vs full-quad", "processed CSV", "equivalence_analysis.csv", "src", "high", "read_full", "Inherits bootstrap bug"),
    ("scripts/06_depth_watermass_southern_ocean.py", "python", "depth, Atlantic water-mass, derivative, SO abyss analyses", "processed CSV", "several csv", "src", "high", "read_full",
     "Derivatives computed on all years incl. external holdout"),
    ("scripts/07_generate_figures.py", "python", "matplotlib figures from results CSVs", "results/processed", "results/figures/*.png", "matplotlib", "low", "read_partial", ""),
    ("scripts/08_generate_tables.py", "python", "LaTeX tables from results CSVs", "results/processed", "results/tables/*.tex", "pandas", "low", "read_partial", ""),
    ("scripts/09_spatial_cruise_validation.py", "python", "select frozen SR candidate; GroupKFold spatial & cruise CV pre-2018", "processed CSV, SR candidates", "spatial_cruise_*.csv, selected yaml", "src", "critical", "read_full",
     "Frozen SR not refit per fold (leakage); overwrites config/selected_symbolic_candidate.yaml"),
    ("scripts/10_finalize_analysis_outputs.py", "python", "recurrence table, Pareto summary, provenance", "results", "tables, REPRODUCIBILITY_FINAL", "src", "medium", "read_full",
     "Hard-coded rank-1 RMSE dict"),
    ("scripts/11_final_model_hierarchy.py", "python", "definitive hierarchy incl. frozen SR, AIC/BIC on train", "processed CSV", "final_model_hierarchy.csv", "src", "critical", "read_full", ""),
    ("scripts/12_final_numeric_audit.py", "python", "assert manuscript numbers vs result files", "results", "FINAL_MANUSCRIPT_NUMERIC_AUDIT.md", "pandas", "medium", "read_full",
     "Asserts exact counts (599, 86) - brittle"),
    ("scripts/13_build_publication_assets.py", "python", "publication figures/tables + empirical PI coverage", "results, processed CSV", "final_latex/*", "matplotlib", "medium", "read_partial",
     "PI from training residual quantiles; no time term"),
    ("scripts/run_all.py", "python", "pipeline orchestrator", "-", "-", "os.system", "low", "read_full", "Stale: runs step 04, omits 09-13, uses bare 'python'"),
    ("experiments/run_symbolic_search.py", "python", "the PySR experiment actually used (48 runs)", "processed CSV", "results/raw/symbolic_regression/*", "src, pysr", "critical", "read_full",
     "Silently caps YAML settings (pop<=20, iters<=50, maxsize<=25)"),
    ("experiments/analyze_pareto.py", "python", "Pareto plots vs reference models", "SR candidates", "figures, frequency csv", "src", "low", "read_partial", ""),
    ("tests/__init__.py", "python", "package marker", "-", "-", "-", "none", "read_full", ""),
    ("tests/test_basins.py", "python", "22 unit tests: basins, splits, models, metrics, TOST, Pareto", "synthetic data", "-", "pytest, src", "medium", "read_partial",
     "No test for bootstrap correctness or family classifier; not executed (sandbox)"),
    ("results/raw/symbolic_regression/*.csv", "data", "PySR candidate libraries / Pareto / run log", "experiments/run_symbolic_search", "-", "-", "high", "profiled", "599 candidates; family labels heuristic"),
    ("results/raw/symbolic_regression/*.yaml", "data", "SR experiment manifest", "-", "-", "-", "medium", "read_full", ""),
    ("results/processed/*.csv", "data", "processed analysis outputs", "scripts 02-13", "-", "-", "high", "profiled", "Temporal metrics include unmodelled time drift (see report 8.1)"),
    ("results/processed/*.md", "doc", "markdown renderings of result tables", "scripts", "-", "-", "medium", "read_partial", ""),
    ("results/processed/*.yaml", "data", "blocked-validation manifest", "scripts/09", "-", "-", "low", "read_full", ""),
    ("results/processed/*.pdf", "pdf", "equations summary PDF", "-", "-", "-", "low", "listed_only", ""),
    ("results/tables/*.csv", "data", "recurrence and Pareto summary tables", "scripts/10", "-", "-", "medium", "profiled", ""),
    ("results/tables/*.tex", "latex", "generated LaTeX tables", "scripts 08/13", "-", "-", "low", "listed_only", ""),
    ("results/figures/*", "image", "generated figures", "scripts 07/13", "-", "-", "low", "listed_only", ""),
    ("results/*.md", "doc", "author audit / decision / risk reports", "-", "-", "-", "high", "read_partial",
     "SCIENTIFIC_DECISION_REPORT and FINAL_RESULTS_SUMMARY contradict each other on Pareto optimality"),
    ("results/*.yaml", "data", "provenance record", "scripts/10", "-", "-", "medium", "read_full", "raw_data_sha256 null (raw absent at time?)"),
    ("logs/*", "log", "PySR run log", "-", "-", "-", "low", "read_partial", ""),
    ("outputs/*", "binary", "PySR checkpoint / hall of fame from one run", "-", "-", "-", "low", "listed_only", ""),
    ("final_latex/sections/*.tex", "latex", "current manuscript sections", "-", "-", "-", "high", "read_partial", "Claims 'robust temporal transferability' (contradicted by drift)"),
    ("final_latex/main.tex", "latex", "manuscript root + abstract", "-", "-", "-", "high", "read_full", ""),
    ("final_latex/*", "latex/asset", "manuscript assets", "-", "-", "-", "low", "listed_only", ""),
    ("submission/manuscript_source/*", "latex/asset", "submission manuscript copy", "-", "-", "-", "medium", "listed_only", "Divergent copy"),
    ("submission/README.md", "doc", "submission package description", "-", "-", "-", "medium", "read_partial", ""),
    ("submission/verification_report.md", "doc", "author verification of manuscript numbers", "-", "-", "-", "medium", "read_partial", "Status B: revision required"),
    ("submission/*", "latex/asset", "submission package", "-", "-", "-", "low", "listed_only", ""),
    ("paper/*", "latex/pdf", "older manuscript copy", "-", "-", "-", "low", "listed_only", "Superseded"),
    ("paper_build/*", "build", "LaTeX build artifacts", "-", "-", "-", "none", "listed_only", ""),
    ("paper_final/src/phase2/*", "python", "G1 Southern Ocean CO2(aq) RF/XGBoost/SHAP/PySR pipeline", "southern_ocean_training.csv", "phase2 results", "sklearn, xgboost, shap, pysr", "low (obsolete target)", "read_partial", "Obsolete"),
    ("paper_final/src/phase3/*", "python", "G1 equation discovery/validation for CO2(aq)", "csv", "figures", "pysr", "low", "read_partial", "Obsolete"),
    ("paper_final/src/*.py", "python", "G1 preprocessing (gsw, CO2(aq))", ".mat", "csv", "gsw", "low", "read_partial", "config.py empty"),
    ("paper_final/*.py", "python", "G2 validation / fix / CANYON-B / Revelle scripts", "GLODAP .mat (absent)", "csv/png in cwd", "scipy, sklearn", "low-medium", "read_partial",
     "Cannot run without raw data; superseded by G3; several self-documented errors"),
    ("paper_final/notebooks/*", "notebook", "exploration notebook", "-", "-", "-", "none", "read_full", "Empty file (0 bytes)"),
    ("paper_final/models/*", "binary", "joblib RF/XGB models (G1)", "-", "-", "-", "none", "listed_only", ""),
    ("paper_final/data/processed/*", "data", "G1 Southern Ocean training set", "-", "-", "-", "low", "profiled", "13,191 rows"),
    ("paper_final/*.md", "doc", "G2 manuscript draft", "-", "-", "-", "low", "listed_only", "Superseded"),
    ("paper_final/*.txt", "doc", "G2 run summaries", "-", "-", "-", "low", "listed_only", ""),
    ("paper_final/*", "artifact", "G1/G2 outputs (csv/png/pdf/latex)", "-", "-", "-", "low", "listed_only", "Superseded"),
    ("other/glodap_canyonb_results.csv", "data", "CANYON-B DIC predictions with pres/oxy", "other/canyonb2.py", "-", "-", "medium", "profiled",
     "Generated with fixed date 2020-01-01 for every row; no year/cruise key"),
    ("other/compressed_data.csv", "data", "unknown", "-", "-", "-", "unknown", "profiled", "Non-UTF8 content; unreferenced"),
    ("other/*.py", "python", "G2 scratch scripts", ".mat", "-", "scipy", "low", "read_partial", "main.py empty; 1.py references undefined variables"),
    ("other/*", "artifact", "scratch outputs", "-", "-", "-", "none", "listed_only", ""),
    ("research_docs/*", "doc", "THIS research program's documentation", "-", "-", "-", "-", "authored", ""),
    ("research/*", "python", "THIS research program's code", "-", "-", "-", "-", "authored", ""),
]


def classify(path: str) -> tuple[str, ...]:
    for rule in RULES:
        if fnmatch.fnmatch(path, rule[0]):
            return rule[1:]
    return ("unknown", "", "", "", "", "", "not_inspected", "no rule matched")


def main() -> None:
    files = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT, text=True).split("\0")
    files += [p for p in subprocess.check_output(
        ["git", "ls-files", "-z", "--others", "--exclude-standard"], cwd=ROOT, text=True).split("\0")]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file_path", "file_type", "purpose", "inputs", "outputs", "dependencies",
                    "scientific_relevance", "inspection_status", "issues"])
        for f in sorted(set(files) - {""}):
            w.writerow([f, *classify(f)])
    print(f"wrote {OUT} ({len(set(files))} files)")


if __name__ == "__main__":
    main()
