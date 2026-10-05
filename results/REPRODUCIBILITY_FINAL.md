# Reproducibility Final Record

This record freezes the software, data, split, and holdout provenance for the final analysis.

project_root: /Users/nikunjmahajan/Desktop/RE/revised_symbolic_tco2
raw_data: /Users/nikunjmahajan/Desktop/RE/revised_symbolic_tco2/data/raw/GLODAPv2.2023_Merged_Master_File.mat
raw_data_sha256: 002881fa71923d15c1bd5aca3b99bc220ad5d1c77a90f130ba1ff0d3910d8766
python: 3.14.3 (main, Feb  3 2026, 15:32:20) [Clang 17.0.0 (clang-1700.6.3.2)]
platform: macOS-15.6.1-arm64-arm-64bit-Mach-O
machine: arm64
numpy: 2.5.2
pandas: 3.0.5
julia: julia version 1.12.7
pysr: 2.3.0
external_holdout_start_year: 2018
selection_data: pre-2018 only; validation/model selection is 2015-2017
external_holdout_locked: true
symbolic_manifest_sha256: 76472797c6e298f8a01c95d4b2ed98eeba9dd2e33b0cfe2ff4c4b19d461e4f87

Exact final-output command sequence:

```bash
python scripts/01_preprocess.py
python scripts/02_baseline_reproduction.py
python scripts/03_model_families.py
python experiments/run_symbolic_search.py
python scripts/09_spatial_cruise_validation.py
python scripts/10_finalize_analysis_outputs.py
python scripts/07_generate_figures.py
python scripts/08_generate_tables.py
```
