"""
Temporal and spatial splitting utilities.

Implements:
- Chronological temporal split (train < 2015, validation 2015-2018, external ≥ 2018)
- Spatial block assignment for cruise/spatial-blocked validation
"""

import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def split_temporal(df, train_end=2015, val_end=2018):
    """
    Create boolean masks for temporal splits.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'year' column.
    train_end : int
        Training includes years < train_end.
    val_end : int
        Validation includes train_end <= year < val_end.
        External includes year >= val_end.

    Returns
    -------
    dict of str -> np.ndarray[bool]
        Keys: 'train', 'validation', 'external'
    """
    yr = df["year"].values
    return {
        "train": yr < train_end,
        "validation": (yr >= train_end) & (yr < val_end),
        "external": yr >= val_end,
    }


def assign_spatial_blocks(df, lat_bin_size=5.0, lon_bin_size=10.0):
    """
    Assign spatial block identifiers based on latitude × longitude bins.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'latitude' and 'longitude'.
    lat_bin_size : float
        Latitude bin size in degrees.
    lon_bin_size : float
        Longitude bin size in degrees.

    Returns
    -------
    pd.Series
        Block identifier string for each row.
    """
    lat = df["latitude"].values
    lon = df["longitude"].values

    # Handle longitude wrapping
    lon = np.where(lon > 180, lon - 360, lon)

    lat_bin = np.floor(lat / lat_bin_size) * lat_bin_size
    lon_bin = np.floor(lon / lon_bin_size) * lon_bin_size

    blocks = [f"{la:.0f}_{lo:.0f}" for la, lo in zip(lat_bin, lon_bin)]
    return pd.Series(blocks, index=df.index)


def non_external_mask(df, external_start_year=2018):
    """Return the immutable eligibility mask used before external evaluation."""
    if "year" not in df.columns:
        raise ValueError("Blocked validation requires a year column")
    year = pd.to_numeric(df["year"], errors="coerce")
    if year.isna().any():
        raise ValueError("Blocked validation cannot use rows with missing year")
    return (year < external_start_year).to_numpy()


def assign_cruise_blocks(df):
    """Assign cruise identifiers without silently falling back to spatial blocks."""
    if "cruise" not in df.columns:
        raise ValueError("Cruise validation requires a cruise identifier column")
    cruise = df["cruise"]
    if cruise.isna().any():
        raise ValueError("Cruise validation requires non-missing cruise identifiers")
    return cruise.astype(str)


def validate_group_disjointness(train_groups, test_groups):
    """Return diagnostics and assert that grouped folds do not overlap."""
    train_set = set(map(str, train_groups))
    test_set = set(map(str, test_groups))
    overlap = sorted(train_set.intersection(test_set))
    if overlap:
        raise AssertionError(f"Groups overlap between train and test: {overlap[:5]}")
    return {
        "n_train_groups": len(train_set),
        "n_test_groups": len(test_set),
        "n_overlapping_groups": 0,
    }


def validate_temporal_exclusion(df, external_start_year=2018):
    """Assert that an evaluation frame excludes the locked external holdout."""
    mask = non_external_mask(df, external_start_year)
    if not bool(mask.all()):
        raise AssertionError(
            f"External holdout leakage detected: {(~mask).sum()} rows with "
            f"year >= {external_start_year}"
        )
    return {
        "n_rows": int(len(df)),
        "min_year": int(df["year"].min()),
        "max_year": int(df["year"].max()),
        "external_start_year": int(external_start_year),
    }


def get_group_diagnostics(groups):
    """Summarize group counts and group-size distribution."""
    sizes = pd.Series(groups).astype(str).value_counts()
    return {
        "n_groups": int(sizes.size),
        "min_group_size": int(sizes.min()),
        "median_group_size": float(sizes.median()),
        "max_group_size": int(sizes.max()),
    }


def get_grouped_kfold_splits(df, n_splits=5, group_col="cruise", random_state=42):
    """
    Generate GroupKFold splits using cruise or spatial blocks.

    Parameters
    ----------
    df : pd.DataFrame
        Data to split.
    n_splits : int
        Number of folds.
    group_col : str
        Column to use for grouping. If not present, uses spatial blocks.
    random_state : int
        Random seed for shuffling.

    Returns
    -------
    list of (train_idx, test_idx) tuples
    """
    from sklearn.model_selection import GroupKFold

    if group_col not in df.columns:
        groups = assign_spatial_blocks(df).values
    else:
        groups = df[group_col].values

    X = np.zeros((len(df), 1))  # dummy features
    y = df["tco2"].values

    gkf = GroupKFold(n_splits=n_splits)
    splits = list(gkf.split(X, y, groups))
    return splits


def verify_no_overlap(splits):
    """
    Verify that temporal splits have no overlap.

    Parameters
    ----------
    splits : dict of str -> np.ndarray[bool]

    Returns
    -------
    dict
        Overlap statistics.
    """
    names = list(splits.keys())
    results = {}
    for i, n1 in enumerate(names):
        for j, n2 in enumerate(names):
            if j <= i:
                continue
            overlap = (splits[n1] & splits[n2]).sum()
            results[f"{n1} ∩ {n2}"] = int(overlap)

    all_ok = all(v == 0 for v in results.values())
    results["no_overlap"] = all_ok
    return results
