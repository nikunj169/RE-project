"""
Equation family classification for symbolic regression candidates.

Classifies candidates by structural form rather than string equality,
accounting for algebraically equivalent representations.
"""

import re
import numpy as np
import logging

logger = logging.getLogger(__name__)


def classify_structure(equation_str):
    """
    Classify an equation into a structural family.

    Parameters
    ----------
    equation_str : str
        The equation string (from PySR).

    Returns
    -------
    str
        One of: 'constant', 'linear', 'rank-1 quadratic', 'full quadratic',
        'polynomial higher-order', 'exponential', 'logarithmic',
        'rational', 'separable nonlinear', 'other'
    """
    eq = equation_str.lower().strip()

    # Detect operators
    has_exp = "exp" in eq
    has_log = "log" in eq
    has_sqrt = "sqrt" in eq
    has_square = ("** 2" in eq or "square" in eq or
                  "sal * sal" in eq or "tmp * tmp" in eq or "aou100 * aou100" in eq)
    has_div = "/" in eq

    # Detect variable usage (support both x0/x1/x2 and sal/tmp/aou100)
    has_sal = "sal" in eq or "x0" in eq
    has_tmp = "tmp" in eq or "x1" in eq
    has_aou = "aou100" in eq or "aou" in eq or "x2" in eq
    n_vars = sum([has_sal, has_tmp, has_aou])

    # Check for constant (no variables)
    if n_vars == 0:
        return "constant"

    # Check for cross-terms between different variables
    has_cross = (
        ("sal * tmp" in eq or "tmp * sal" in eq or
         "sal * aou" in eq or "aou * sal" in eq or
         "tmp * aou" in eq or "aou * tmp" in eq) or
        # Also check x0/x1/x2 pattern
        any(f"x{i} * x{j}" in eq or f"x{j} * x{i}" in eq
            for i in range(3) for j in range(i + 1, 3))
    )

    # Check for linear (only +, -, * with constants, no squares/cross/ex/log)
    if not has_square and not has_exp and not has_log and not has_sqrt:
        if not has_cross:
            return "linear"

    # Check for rank-1 quadratic pattern
    if has_square and not has_cross:
        return "rank-1 quadratic"

    # Check for full quadratic
    if has_square and has_cross:
        return "full quadratic"

    # Higher-order polynomial detection
    if "** 3" in eq or "cube" in eq or "aou100 * aou100 * aou100" in eq:
        return "polynomial higher-order"

    if has_exp and has_log:
        return "separable nonlinear"
    if has_exp:
        return "exponential"
    if has_log:
        return "logarithmic"
    if has_div:
        return "rational"
    if has_sqrt:
        return "root"

    return "other"


def classify_candidates(candidates_df, equation_col="equation"):
    """
    Classify all candidates in a DataFrame.

    Parameters
    ----------
    candidates_df : pd.DataFrame
    equation_col : str

    Returns
    pd.Series
        Family classification for each candidate.
    """
    return candidates_df[equation_col].apply(classify_structure)


def compute_family_recurrence(candidates_list):
    """
    Compute how often each equation family appears across multiple SR runs.

    Parameters
    ----------
    candidates_list : list of pd.DataFrame
        Each DataFrame from a different SR run (same basin/search space).

    Returns
    -------
    pd.DataFrame
        Family recurrence statistics.
    """
    from collections import Counter

    all_families = []
    for df in candidates_list:
        if df is not None and "family" in df.columns:
            # Count families among Pareto-optimal candidates
            pareto = df.get("is_pareto", pd.Series(True, index=df.index))
            families = df.loc[pareto, "family"].tolist()
            all_families.append(Counter(families))

    if not all_families:
        return pd.DataFrame()

    # Aggregate across runs
    total_counts = Counter()
    for c in all_families:
        total_counts.update(c)

    n_runs = len(all_families)
    results = []
    for family, count in total_counts.items():
        freq_per_run = count / n_runs
        results.append({
            "family": family,
            "total_count": count,
            "mean_per_run": freq_per_run,
            "fraction": count / sum(total_counts.values()),
        })

    return pd.DataFrame(results).sort_values("total_count", ascending=False)
