"""
GLODAP v2.2023 data loading and quality control.

Loads the .mat file, extracts required variables, applies QC filters,
and produces a clean DataFrame for analysis.
"""

import numpy as np
import pandas as pd
import scipy.io as sio
import os
import logging

from .config_loader import load_config, get_data_path

logger = logging.getLogger(__name__)


def load_glodap_raw(mat_path=None):
    """
    Load raw GLODAP variables from the .mat file.

    Parameters
    ----------
    mat_path : str or Path, optional
        Path to the .mat file. If None, uses config default.

    Returns
    -------
    dict of str -> np.ndarray
        Variable name to array mapping.
    """
    if mat_path is None:
        mat_path = get_data_path()

    mat_path = str(mat_path)
    logger.info(f"Loading GLODAP data from {mat_path}")

    try:
        mat = sio.loadmat(mat_path, squeeze_me=True)
    except NotImplementedError:
        import h5py
        logger.info("Detected v7.3 MAT file, using h5py")
        mat = {}
        with h5py.File(mat_path, "r") as f:
            for k, v in f.items():
                if not k.startswith("__"):
                    mat[k] = np.array(v)

    config = load_config()
    var_config = config["variables"]

    data = {}
    for key, mat_name in var_config["predictors"].items():
        if mat_name in mat:
            data[key] = np.squeeze(mat[mat_name]).astype(float)
        else:
            logger.warning(f"Variable {mat_name} not found in .mat file")

    target_name = var_config["target"]
    if target_name in mat:
        data["tco2"] = np.squeeze(mat[target_name]).astype(float)
    else:
        logger.warning(f"Target variable {target_name} not found")

    for key, mat_name in var_config["metadata"].items():
        if mat_name in mat:
            data[key] = np.squeeze(mat[mat_name]).astype(float)
        else:
            logger.warning(f"Metadata variable {mat_name} not found")

    n = len(next(iter(data.values())))
    logger.info(f"Loaded {n:,} raw observations")
    return data


def apply_quality_control(data, config=None):
    """
    Apply quality control filters to raw GLODAP data.

    Parameters
    ----------
    data : dict
        Raw variable arrays from load_glodap_raw().
    config : dict, optional
        QC config. If None, loads from global config.

    Returns
    -------
    pd.DataFrame
        QC-filtered data with all variables.
    dict
        QC statistics (counts before/after, per-filter rejections).
    """
    if config is None:
        config = load_config()

    qc = config["quality_control"]

    # Build initial finite mask
    all_vars = ["salinity", "temperature", "aou", "tco2"]
    mask = np.ones(len(data["salinity"]), dtype=bool)

    if qc.get("require_finite", True):
        for var in all_vars:
            if var in data:
                mask &= np.isfinite(data[var])

    n_finite = mask.sum()

    # Apply range filters
    filter_stats = {"total_raw": len(data["salinity"]), "finite": int(n_finite)}

    for var, bounds in qc.items():
        if var in ("require_finite",):
            continue
        if isinstance(bounds, dict) and var in data:
            if "min" in bounds:
                mask &= data[var] > bounds["min"]
            if "max" in bounds:
                mask &= data[var] < bounds["max"]
            filter_stats[f"{var}_range"] = int(mask.sum())

    n_after_qc = mask.sum()
    filter_stats["after_qc"] = int(n_after_qc)

    # Build DataFrame from filtered data
    df_data = {}
    for var in all_vars:
        if var in data:
            df_data[var] = data[var][mask]

    # Add metadata
    for key in ["latitude", "longitude", "depth", "year", "month", "region"]:
        if key in data:
            df_data[key] = data[key][mask]

    # Add cruise/station if available (for spatial blocking)
    for key in ["cruise", "station"]:
        if key in data:
            df_data[key] = data[key][mask]

    df = pd.DataFrame(df_data)

    logger.info(f"QC: {n_after_qc:,} / {len(data['salinity']):,} observations retained "
                f"({100*n_after_qc/len(data['salinity']):.1f}%)")

    return df, filter_stats


def load_and_qc(mat_path=None):
    """
    Convenience function: load raw data and apply QC in one step.

    Returns
    -------
    pd.DataFrame
        QC-filtered data.
    dict
        QC statistics.
    """
    raw = load_glodap_raw(mat_path)
    df, stats = apply_quality_control(raw)
    return df, stats


def save_processed(df, output_path=None):
    """Save processed DataFrame to CSV."""
    if output_path is None:
        from .config_loader import get_path
        output_dir = get_path("processed_data")
        os.makedirs(output_dir, exist_ok=True)
        output_path = output_dir / "glodap_qc_filtered.csv"

    df.to_csv(output_path, index=False)
    logger.info(f"Saved processed data to {output_path}")
    return output_path
