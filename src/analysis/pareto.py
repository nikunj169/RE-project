"""
Pareto frontier computation for symbolic regression candidates.

Given a set of candidate equations with complexity and validation metrics,
compute the Pareto-optimal frontier (best RMSE at each complexity level).
"""

import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def compute_pareto_frontier(df, complexity_col="complexity", metric_col="RMSE"):
    """
    Compute the Pareto-optimal frontier.

    A candidate is Pareto-optimal if no other candidate has both
    lower complexity AND lower (or equal) metric value.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain complexity_col and metric_col.
    complexity_col, metric_col : str

    Returns
    -------
    pd.DataFrame
        Pareto-optimal subset, sorted by complexity.
    """
    df = df.dropna(subset=[complexity_col, metric_col]).copy()
    df = df.sort_values(complexity_col)

    pareto_mask = np.zeros(len(df), dtype=bool)
    best_metric = np.inf

    for i, (_, row) in enumerate(df.iterrows()):
        if row[metric_col] < best_metric:
            pareto_mask[i] = True
            best_metric = row[metric_col]

    pareto = df[pareto_mask].copy()
    logger.info(f"Pareto frontier: {len(pareto)} / {len(df)} candidates")
    return pareto


def classify_equation_family(equation_str):
    """
    Classify a symbolic equation into a structural family.

    This is a heuristic text-based classifier. Algebraically equivalent
    equations may have different textual representations.

    Parameters
    ----------
    equation_str : str
        The equation string from PySR.

    Returns
    -------
    str
        Family classification.
    """
    eq = equation_str.lower().strip()

    # Check for specific structural patterns
    if "x0" not in eq and "x1" not in eq and "x2" not in eq:
        return "constant"

    # Count occurrences of squared patterns
    has_square = ("**2" in eq or "^2" in eq or
                  any(f"x{i} * x{i}" in eq for i in range(3)))

    has_cross = any(
        f"x{i} * x{j}" in eq or f"x{j} * x{i}" in eq
        for i in range(3) for j in range(i+1, 3)
    )

    has_exp = "exp" in eq
    has_log = "log" in eq
    has_sqrt = "sqrt" in eq
    has_div = "/" in eq

    # Check for linear
    if not has_square and not has_cross and not has_exp and not has_log:
        if "+" in eq or "-" in eq:
            return "linear"

    # Check for rank-1 quadratic pattern: (a*x0 + b*x1 + c*x2 + d)^2 + e
    if has_square and not has_cross:
        return "rank-1 quadratic"

    if has_square and has_cross:
        return "full quadratic"

    if has_exp:
        return "exponential"

    if has_log:
        return "logarithmic"

    if has_div:
        return "rational"

    if has_sqrt:
        return "root"

    return "other"


def compute_pareto_with_family(candidates_df, S_val, T_val, Y_val, predict_fn,
                                complexity_col="complexity"):
    """
    Compute Pareto frontier and classify equation families.

    Parameters
    ----------
    candidates_df : pd.DataFrame
        Candidate equations with complexity and equation columns.
    S_val, T_val, Y_val : np.ndarray
        Validation data for computing RMSE.
    predict_fn : callable
        Function(equation_str_or_idx, S, T, A) -> predictions.

    Returns
    -------
    pd.DataFrame
        Candidates with RMSE, family classification, and Pareto status.
    """
    results = []
    for i, row in candidates_df.iterrows():
        try:
            y_pred = predict_fn(i, S_val, T_val, 0)  # AOU already in features
            rmse = np.sqrt(np.mean((Y_val - y_pred) ** 2))
        except Exception:
            rmse = np.nan

        family = classify_equation_family(str(row.get("equation", "")))

        results.append({
            "candidate_idx": i,
            "complexity": row[complexity_col],
            "loss": row.get("loss", np.nan),
            "RMSE_val": rmse,
            "family": family,
            "equation": str(row.get("equation", "")),
        })

    results_df = pd.DataFrame(results)

    # Compute Pareto frontier
    pareto = compute_pareto_frontier(results_df, complexity_col, "RMSE_val")
    results_df["is_pareto"] = results_df.index.isin(pareto.index)

    return results_df
