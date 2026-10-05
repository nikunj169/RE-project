# Julia Environment Versions

## Verified Versions (2026-09-10)

- **Python**: 3.14.3 (Homebrew, Apple Silicon)
- **Julia**: 1.12.7 (Homebrew)
- **PySR**: 2.3.0
- **JuliaCall**: 0.9.35
- **SymbolicRegression.jl**: ~2.4 (managed by PySR)
- **scipy**: 1.18.1
- **numpy**: 2.5.2
- **scikit-learn**: 1.9.0
- **pandas**: 3.0.5

## Installation Commands

```bash
# Julia
brew install julia

# Python packages
pip install pysr scipy numpy pandas scikit-learn matplotlib seaborn pyyaml tabulate jinja2

# PySR Julia backend
python -m pysr install
```

## Notes

- PySR 2.3.0 requires `tournament_selection_n < population_size`
- For reproducibility, set `deterministic=True` and `parallelism='serial'` in PySRRegressor
- The numpy directory in the parent RE/ folder can shadow the system numpy; always run from the project directory
