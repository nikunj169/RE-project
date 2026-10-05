"""
Basin definitions and masking.

Implements the mutually exclusive basin definitions:
  Atlantic:       region code 1, latitude > -35°
  Indian:         region code 16, latitude > -35°
  Pacific:        region code 8, latitude > -35°
  Southern Ocean: latitude < -35°, all region codes
"""

import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)

# Basin definitions (matches config/basins.yaml)
BASIN_DEFS = {
    "Atlantic": {"codes": [1], "lat_min": -35.0},
    "Indian": {"codes": [16], "lat_min": -35.0},
    "Pacific": {"codes": [8], "lat_min": -35.0},
    "Southern Ocean": {"codes": None, "lat_max": -35.0},
}


def assign_basin(df, config=None):
    """
    Assign each row to a basin based on region code and latitude.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'region' and 'latitude' columns.
    config : dict, optional
        Basin config dict. If None, uses built-in definitions.

    Returns
    -------
    pd.Series
        Basin name for each row.
    """
    if config is None:
        defs = BASIN_DEFS
    else:
        defs = {}
        for name, cfg in config.items():
            # Normalize config keys to match internal format
            normalized = dict(cfg)
            if "region_codes" in normalized:
                normalized["codes"] = normalized.pop("region_codes")
            # Normalize basin name (SouthernOcean -> Southern Ocean)
            display_name = name.replace("_", " ")
            if display_name == "SouthernOcean":
                display_name = "Southern Ocean"
            defs[display_name] = normalized

    basin = pd.Series("Unclassified", index=df.index)

    # Southern Ocean first (lat-based, captures everything south)
    if "Southern Ocean" in defs:
        so_cfg = defs["Southern Ocean"]
        mask = df["latitude"] < so_cfg["lat_max"]
        basin[mask] = "Southern Ocean"

    # Then region-based basins (override SO if lat > threshold)
    for name in ["Atlantic", "Indian", "Pacific"]:
        if name in defs:
            cfg = defs[name]
            codes = cfg.get("codes", [])
            lat_min = cfg.get("lat_min", -90)
            mask = np.isin(df["region"], codes) & (df["latitude"] >= lat_min)
            basin[mask] = name

    return basin


def get_basin_mask(df, basin_name, config=None):
    """
    Get boolean mask for a specific basin.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'region' and 'latitude'.
    basin_name : str
        One of: 'Atlantic', 'Indian', 'Pacific', 'Southern Ocean'.
    config : dict, optional
        Basin config dict.

    Returns
    -------
    np.ndarray of bool
    """
    basins = assign_basin(df, config)
    return (basins == basin_name).values


def verify_mutual_exclusivity(df, config=None):
    """
    Verify that basin definitions are mutually exclusive.

    Returns
    -------
    dict
        Overlap statistics.
    """
    basins = assign_basin(df, config)
    basin_names = ["Atlantic", "Indian", "Pacific", "Southern Ocean"]

    overlaps = {}
    for i, b1 in enumerate(basin_names):
        for j, b2 in enumerate(basin_names):
            if j <= i:
                continue
            m1 = (basins == b1).values
            m2 = (basins == b2).values
            overlap = (m1 & m2).sum()
            overlaps[f"{b1} ∩ {b2}"] = int(overlap)

    n_total = (basins != "Unclassified").sum()
    n_sum = sum((basins == b).sum() for b in basin_names)

    result = {
        "overlaps": overlaps,
        "total_unique": int(n_total),
        "sum_of_basins": int(n_sum),
        "mutually_exclusive": n_total == n_sum,
    }

    if result["mutually_exclusive"]:
        logger.info("Basin definitions are mutually exclusive ✓")
    else:
        logger.warning("Basin definitions have OVERLAP — check definitions!")

    return result


def get_basin_counts(df):
    """Get observation counts per basin and per temporal split."""
    from .splits import split_temporal

    basins = assign_basin(df)
    splits = split_temporal(df)

    counts = []
    for basin_name in ["Atlantic", "Indian", "Pacific", "Southern Ocean"]:
        m = (basins == basin_name).values
        sub = df[m]
        sub_splits = {k: v[m] for k, v in splits.items()}
        counts.append({
            "Basin": basin_name,
            "N_total": int(m.sum()),
            "N_train": int(sub_splits["train"].sum()),
            "N_validation": int(sub_splits["validation"].sum()),
            "N_external": int(sub_splits["external"].sum()),
        })

    total = {k: sum(c[k] for c in counts) for k in ["N_total", "N_train", "N_validation", "N_external"]}
    total["Basin"] = "TOTAL"
    counts.append(total)

    return pd.DataFrame(counts)
