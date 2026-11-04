"""Shared loading of config and the research dataset."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def load_config() -> dict:
    return yaml.safe_load((ROOT / "research" / "configs" / "main.yaml").read_text())


def load_data(cfg: dict, strict: bool = True) -> pd.DataFrame:
    """Load the research dataset, restricted to study basins and finite gamma.

    Adds ``layer`` (neutral-density interval label) and ``cell`` (cfg cell_deg grid).
    """
    from research.models.within_cell_trend import add_cells

    df = pd.read_csv(ROOT / cfg["data"]["path"])
    df = df[df["basin"].isin(cfg["basins"]) & df["gamma"].notna()]
    if strict:
        df = df[df["qc_strict"]]
    df = df.copy()
    df["layer"] = pd.cut(df["gamma"], cfg["gamma_edges"]).astype(str)
    df["cell"] = add_cells(df, cfg["cell_deg"])
    return df


def groups(df: pd.DataFrame, cfg: dict, min_span: float | None = None):
    """Yield (basin, layer, subset) for groups meeting the inclusion rules."""
    from research.models.within_cell_trend import select_cells

    min_span = cfg["min_span_years"] if min_span is None else min_span
    for (basin, layer), g in df.groupby(["basin", "layer"], sort=False):
        g = select_cells(g, min_span)
        if len(g) >= cfg["min_samples"] and g["cell"].nunique() >= cfg["min_cells"]:
            yield basin, layer, g
