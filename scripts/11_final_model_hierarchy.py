"""Create the definitive temporal/external model hierarchy and information criteria."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.analysis.metrics import compute_metrics
from src.basins import assign_basin
from src.config_loader import get_path, load_config
from src.models.cubic import CubicPolynomial
from src.models.frozen_symbolic import FrozenSymbolicModel
from src.models.full_quadratic import FullQuadratic
from src.models.linear import LinearModel
from src.models.mean_baseline import MeanBaseline
from src.models.rank1_quadratic import Rank1Quadratic
from src.models.rank2_quadratic import Rank2Quadratic

BASINS = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]
FACTORIES = {
    "mean": MeanBaseline,
    "linear": LinearModel,
    "rank1": Rank1Quadratic,
    "rank2": Rank2Quadratic,
    "full_quadratic": FullQuadratic,
    "cubic": CubicPolynomial,
}
PARAMETERS = {
    "mean": 1,
    "linear": 4,
    "rank1": 5,
    "rank2": 9,
    "full_quadratic": 10,
    "cubic": 20,
}
COMPLEXITIES = {
    "mean": 1,
    "linear": 4,
    "rank1": 8,
    "rank2": 16,
    "full_quadratic": 20,
    "cubic": 30,
}


def info_criteria(y_true, y_pred, parameter_count):
    residual = np.asarray(y_true, dtype=float) - np.asarray(y_pred, dtype=float)
    n = len(residual)
    rss = float(np.sum(residual**2))
    sigma2 = max(rss / n, np.finfo(float).tiny)
    loglik = -0.5 * n * (np.log(2 * np.pi * sigma2) + 1)
    return {
        "RSS_train": rss,
        "loglik_train": loglik,
        "AIC_train": -2 * loglik + 2 * parameter_count,
        "BIC_train": -2 * loglik + parameter_count * np.log(n),
    }


def load_selected():
    path = Path(__file__).resolve().parents[1] / "config" / "selected_symbolic_candidate.yaml"
    with open(path) as handle:
        payload = yaml.safe_load(handle)
    return {row["basin"]: row for row in payload["candidates"]}


def main():
    root = Path(__file__).resolve().parents[1]
    config = load_config()
    data_path = get_path("processed_data") / "glodap_qc_filtered.csv"
    df = pd.read_csv(data_path)
    df["basin"] = assign_basin(df, config["basins"])
    selected = load_selected()
    rows = []

    for basin in BASINS:
        sub = df[df["basin"] == basin].copy()
        train = sub[sub.year < 2015]
        validation = sub[(sub.year >= 2015) & (sub.year < 2018)]
        external = sub[sub.year >= 2018]
        for model_name, factory in FACTORIES.items():
            model = factory()
            model.fit(train.salinity, train.temperature, train.aou, train.tco2)
            train_pred = model.predict(train.salinity, train.temperature, train.aou)
            val_pred = model.predict(validation.salinity, validation.temperature, validation.aou)
            ext_pred = model.predict(external.salinity, external.temperature, external.aou)
            train_metrics = compute_metrics(train.tco2, train_pred)
            val_metrics = compute_metrics(validation.tco2, val_pred)
            ext_metrics = compute_metrics(external.tco2, ext_pred)
            row = {
                "basin": basin,
                "model": model_name,
                "parameter_count": PARAMETERS[model_name],
                "expression_complexity": COMPLEXITIES[model_name],
                "n_train": len(train),
                "n_validation": len(validation),
                "n_external": len(external),
                "train_RMSE": train_metrics["RMSE"],
                "train_MAE": train_metrics["MAE"],
                "train_R2": train_metrics["R2"],
                "train_Bias": train_metrics["Bias"],
                "validation_RMSE": val_metrics["RMSE"],
                "validation_MAE": val_metrics["MAE"],
                "validation_R2": val_metrics["R2"],
                "validation_Bias": val_metrics["Bias"],
                "external_RMSE": ext_metrics["RMSE"],
                "external_MAE": ext_metrics["MAE"],
                "external_R2": ext_metrics["R2"],
                "external_Bias": ext_metrics["Bias"],
                **info_criteria(train.tco2, train_pred, PARAMETERS[model_name]),
                "aic_bic_comparable": True,
                "selection_source": "parametric_family_predefined",
            }
            rows.append(row)

        candidate = selected[basin]
        frozen = FrozenSymbolicModel(candidate["equation"], candidate["sympy_format"], candidate["complexity"])
        frozen.fit(train.salinity, train.temperature, train.aou, train.tco2)
        train_pred = frozen.predict(train.salinity, train.temperature, train.aou)
        val_pred = frozen.predict(validation.salinity, validation.temperature, validation.aou)
        ext_pred = frozen.predict(external.salinity, external.temperature, external.aou)
        train_metrics = compute_metrics(train.tco2, train_pred)
        val_metrics = compute_metrics(validation.tco2, val_pred)
        ext_metrics = compute_metrics(external.tco2, ext_pred)
        rows.append(
            {
                "basin": basin,
                "model": "frozen_symbolic",
                "parameter_count": np.nan,
                "expression_complexity": candidate["complexity"],
                "n_train": len(train),
                "n_validation": len(validation),
                "n_external": len(external),
                "train_RMSE": train_metrics["RMSE"],
                "train_MAE": train_metrics["MAE"],
                "train_R2": train_metrics["R2"],
                "train_Bias": train_metrics["Bias"],
                "validation_RMSE": val_metrics["RMSE"],
                "validation_MAE": val_metrics["MAE"],
                "validation_R2": val_metrics["R2"],
                "validation_Bias": val_metrics["Bias"],
                "external_RMSE": ext_metrics["RMSE"],
                "external_MAE": ext_metrics["MAE"],
                "external_R2": ext_metrics["R2"],
                "external_Bias": ext_metrics["Bias"],
                "RSS_train": np.nan,
                "loglik_train": np.nan,
                "AIC_train": np.nan,
                "BIC_train": np.nan,
                "aic_bic_comparable": False,
                "selection_source": "temporal_validation_only_locked_external",
                "candidate_id": candidate["candidate_id"],
                "equation": candidate["equation"],
            }
        )

    result = pd.DataFrame(rows)
    rank1 = result[result.model == "rank1"].set_index("basin")["validation_RMSE"]
    result["validation_delta_vs_rank1"] = result.apply(lambda r: r.validation_RMSE - rank1[r.basin], axis=1)
    result["validation_pct_delta_vs_rank1"] = result["validation_delta_vs_rank1"] / result.apply(lambda r: rank1[r.basin], axis=1) * 100
    result["external_delta_vs_rank1"] = result.apply(lambda r: r.external_RMSE - result[(result.basin == r.basin) & (result.model == "rank1")].external_RMSE.iloc[0], axis=1)
    result["external_pct_delta_vs_rank1"] = result["external_delta_vs_rank1"] / result.apply(lambda r: result[(result.basin == r.basin) & (result.model == "rank1")].external_RMSE.iloc[0], axis=1) * 100
    output = get_path("results_processed") / "final_model_hierarchy.csv"
    result.to_csv(output, index=False)
    result.to_markdown(get_path("results_processed") / "final_model_hierarchy.md", index=False)
    print(output)


if __name__ == "__main__":
    main()
