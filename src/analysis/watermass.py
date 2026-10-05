"""
Water-mass classification and analysis.

Implements threshold-based water-mass classification following
oceanographic conventions (Tomczak, 2002).
"""

import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)

# Default thresholds from config/watermass.yaml
WATER_MASS_DEFS = {
    "Surface": {"depth_max": 200},
    "Thermocline": {"depth_min": 200, "depth_max": 800, "temperature_min": 8.0},
    "AAIW": {
        "temperature_min": 1.0, "temperature_max": 8.0,
        "salinity_min": 33.5, "salinity_max": 34.6,
        "depth_min": 300, "depth_max": 1500,
        "latitude_max": -15.0,
    },
    "NADW": {
        "temperature_min": 1.0, "temperature_max": 5.0,
        "salinity_min": 34.8, "salinity_max": 35.1,
        "depth_min": 1000, "depth_max": 4000,
    },
    "AABW": {
        "temperature_max": 1.0,
        "salinity_min": 34.6, "salinity_max": 34.75,
        "depth_min": 3000,
    },
}


def classify_water_mass(df, defs=None):
    """
    Assign water-mass labels to each observation.

    Priority: AABW > NADW > AAIW > Thermocline > Surface
    (more specific classifications override more general ones).

    Returns
    -------
    pd.Series
        Water-mass label for each row.
    """
    if defs is None:
        defs = WATER_MASS_DEFS

    wm = pd.Series("Unclassified", index=df.index)

    # Apply in reverse priority order (Surface last)
    for name in ["Surface", "Thermocline", "AAIW", "NADW", "AABW"]:
        if name not in defs:
            continue
        cfg = defs[name]
        mask = np.ones(len(df), dtype=bool)

        if "depth_min" in cfg:
            mask &= df["depth"] >= cfg["depth_min"]
        if "depth_max" in cfg:
            mask &= df["depth"] < cfg["depth_max"]
        if "temperature_min" in cfg:
            mask &= df["temperature"] >= cfg["temperature_min"]
        if "temperature_max" in cfg:
            mask &= df["temperature"] < cfg["temperature_max"]
        if "salinity_min" in cfg:
            mask &= df["salinity"] >= cfg["salinity_min"]
        if "salinity_max" in cfg:
            mask &= df["salinity"] < cfg["salinity_max"]
        if "latitude_max" in cfg:
            mask &= df["latitude"] < cfg["latitude_max"]
        if "latitude_min" in cfg:
            mask &= df["latitude"] >= cfg["latitude_min"]

        wm[mask] = name

    return wm


def watermass_performance(y_true, y_pred, df, defs=None):
    """
    Compute prediction metrics per water mass.

    Returns
    -------
    pd.DataFrame
        One row per water mass with metrics.
    """
    from .metrics import compute_metrics

    wm = classify_water_mass(df, defs)
    results = []

    for name in ["Surface", "Thermocline", "AAIW", "NADW", "AABW", "Unclassified"]:
        mask = (wm == name).values
        if mask.sum() < 10:
            continue
        m = compute_metrics(y_true[mask], y_pred[mask])
        m["water_mass"] = name
        results.append(m)

    return pd.DataFrame(results)
