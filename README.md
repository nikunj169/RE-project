# How Low-Dimensional Is the S–T–AOU Dependence of Ocean Total Dissolved Inorganic Carbon?

## A Symbolic Regression Analysis (Revised Study)

This repository contains the complete analysis pipeline for the revised study investigating the dimensional structure of the empirical relationship between ocean total dissolved inorganic carbon (TCO2) and salinity (S), temperature (T), and apparent oxygen utilisation (AOU).

## Central Scientific Question

**How low-dimensional is the empirical S–T–AOU dependence of ocean TCO2, and how much predictive accuracy can be retained by imposing increasingly strong structural constraints?**

## Key Differences from the Original Study

1. **Model families, not single equations**: Explicit comparison of linear, rank-1 quadratic, rank-2 quadratic, full quadratic, cubic, and symbolic regression candidates
2. **Pareto frontier analysis**: Symbolic regression candidates evaluated on accuracy-complexity tradeoff
3. **Multiple operator spaces**: Systematic investigation across different PySR operator sets
4. **Repeated symbolic regression**: Stability across 8-10 random seeds per configuration
5. **Spatial/cruise-blocked validation**: Not just temporal holdout
6. **Dependence-aware inference**: Cluster bootstrap, not naive per-observation TOST
7. **Equivalence margin sensitivity**: Multiple margins (1, 2, 5, 10 μmol/kg)
8. **Enhanced external holdout analysis**: Detailed basin/depth/water-mass breakdown
9. **No hard-coded results**: All numbers generated from code

## Quick Start

```bash
# 1. Install environment
cd environment
conda env create -f environment.yml
conda activate revised_tco2

# 2. Run the final analysis/output sequence
python scripts/01_preprocess.py
python scripts/02_baseline_reproduction.py
python scripts/03_model_families.py
python experiments/run_symbolic_search.py
python scripts/09_spatial_cruise_validation.py
python scripts/11_final_model_hierarchy.py
python scripts/10_finalize_analysis_outputs.py
python scripts/07_generate_figures.py
python scripts/08_generate_tables.py
python scripts/12_final_numeric_audit.py

# 4. Compile the manuscript
mkdir -p paper_build
tectonic -X compile --keep-logs --keep-intermediates --outdir paper_build paper/latex/main.tex

# 5. Run tests
python3 -m pytest tests/ -q
```

## Directory Structure

```
revised_symbolic_tco2/
├── README.md                      # This file
├── environment/                   # Conda/pip environment files
├── config/                        # YAML configuration files
│   ├── global_config.yaml         # Paths, seeds, splits
│   ├── model_families.yaml        # Model family definitions
│   ├── pysr_config.yaml           # PySR search configurations
│   └── basins.yaml                # Basin definitions
├── data/
│   ├── raw/                       # Raw GLODAP .mat file (symlink)
│   └── processed/                 # QC-filtered, split datasets
├── src/
│   ├── __init__.py
│   ├── data_loader.py             # GLODAP loading and QC
│   ├── basins.py                  # Basin definitions and masking
│   ├── splits.py                  # Temporal/spatial splitting
│   ├── models/
│   │   ├── __init__.py
│   │   ├── base.py                # Abstract model interface
│   │   ├── mean_baseline.py       # Model 0: training mean
│   │   ├── linear.py              # Model 1: linear regression
│   │   ├── rank1_quadratic.py     # Model 2: rank-1 quadratic
│   │   ├── rank2_quadratic.py     # Model 3: rank-2 quadratic
│   │   ├── full_quadratic.py      # Model 4: full 2nd-order polynomial
│   │   ├── cubic.py               # Model 5: cubic benchmark
│   │   └── symbolic.py            # Model 6+: PySR symbolic candidates
│   ├── analysis/
│   │   ├── __init__.py
│   │   ├── metrics.py             # RMSE, MAE, R2, bias, etc.
│   │   ├── pareto.py              # Pareto frontier computation
│   │   ├── depth_analysis.py      # Depth-stratified analysis
│   │   ├── watermass.py           # Water-mass classification
│   │   ├── derivatives.py         # Derivative/sensitivity analysis
│   │   ├── family_classification.py # Equation family classification
│   │   └── southern_ocean.py      # Southern Ocean failure analysis
│   ├── validation/
│   │   ├── __init__.py
│   │   ├── temporal.py            # Temporal holdout validation
│   │   ├── spatial.py             # Spatial/cruise-blocked validation
│   │   ├── equivalence.py         # TOST and equivalence testing
│   │   ├── bootstrap.py           # Dependence-aware bootstrap
│   │   └── prediction_intervals.py # PI construction and calibration
│   └── visualization/
│       ├── __init__.py
│       ├── pareto_plots.py        # Pareto frontier figures
│       ├── comparison_plots.py    # Model family comparison
│       ├── scatter_plots.py       # Observed vs predicted
│       ├── spatial_plots.py       # Spatial validation maps
│       ├── depth_plots.py         # Depth-stratified figures
│       ├── watermass_plots.py     # Water-mass analysis figures
│       └── stability_plots.py     # Stability/robustness figures
├── scripts/                       # End-to-end analysis scripts
├── experiments/                   # Experiment configurations/logs
├── results/
│   ├── raw/                       # Unprocessed experiment outputs
│   ├── processed/                 # Cleaned, aggregated results
│   ├── tables/                    # Publication-ready tables
│   └── figures/                   # Publication-ready figures
├── models/                        # Saved model objects
├── paper/
│   ├── latex/                     # LaTeX manuscript
│   ├── figures/                   # Paper figure copies
│   └── tables/                    # Paper table copies
├── supplementary/                 # Supplementary material
├── tests/                         # Automated tests
└── logs/                          # Experiment logs
```

## Data

The analysis uses GLODAP v2.2023 (Lauvset et al., 2023). The raw `.mat` file should be placed at:
```
data/raw/GLODAPv2.2023_Merged_Master_File.mat
```
Or create a symlink to the existing file in the old project.

## Configuration

All analysis parameters are defined in YAML files under `config/`. Key configuration:

- **Quality Control**: S∈(25,42), T∈(-2.5,35°C), TCO2∈(1700,2600), AOU>-50
- **Basins**: Atlantic (region 1, lat>-35°), Indian (region 16, lat>-35°), Pacific (region 8, lat>-35°), Southern Ocean (lat<-35°)
- **Temporal Split**: Train<2015, Validation 2015-2018, External≥2018
- **Random Seeds**: 42, 123, 456, 789, 1024, 2048, 4096, 8192

## Dependencies

- Python ≥3.11
- Julia ≥1.9 (for PySR)
- PySR ≥1.0
- scipy, numpy, pandas, scikit-learn, matplotlib, seaborn
- gsw (Gibbs SeaWater)
- xarray (optional, for spatial operations)
- See `environment/environment.yml` for full specification

## Reproducibility

All random seeds are fixed and documented in `config/global_config.yaml`.
Every numerical result in the paper is generated by code and saved to `results/`.
No results are manually hard-coded.

## License

[To be determined]
