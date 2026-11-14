"""Unit tests for the research/ trend estimators (Round 2)."""

import numpy as np
import pandas as pd

from research.models.differential_trend import exact_twoway, fit_differential, twoway_demean
from research.models.within_cell_trend import fit_trend


def test_exact_twoway_matches_alternating_projections():
    rng = np.random.default_rng(0)
    n = 4000
    a, b = rng.integers(0, 200, n), rng.integers(0, 30, n)
    w = rng.integers(1, 3, n).astype(float)
    X = rng.normal(size=(n, 3))
    Y = (X @ [1.0, -2.0, 0.5] + rng.normal(size=200)[a] + rng.normal(size=30)[b] + 0.1 * rng.normal(size=n))[:, None]
    c = exact_twoway(X, Y, a, b, w)
    M = twoway_demean(np.hstack([X, Y]), a, b, w, tol=1e-13, max_iter=5000)
    sw = np.sqrt(w)[:, None]
    c2 = np.linalg.lstsq(M[:, :3] * sw, M[:, 3:] * sw, rcond=None)[0]
    assert np.allclose(c, c2, atol=1e-8)


def _synthetic(drift: float, seed: int = 1) -> tuple[pd.DataFrame, np.ndarray]:
    rng = np.random.default_rng(seed)
    rows = []
    for cr in range(60):
        year = 1975 + 45 * rng.random()
        off = drift * (year - 2000) + rng.normal(0, 3)
        for cell in rng.choice(20, 5, replace=False):
            for target in (True, False):
                th, s, aou = rng.normal(10, 3), rng.normal(35, 0.5), rng.normal(100, 30)
                beta = 0.8 if target else 0.1
                y = 2000 + 3 * th - 5 * s + 0.7 * aou + beta * (year - 2000) + off + cell + rng.normal(0, 2)
                rows.append((cr, f"{cell}", year, th, s, aou, y, target))
    df = pd.DataFrame(rows, columns=["cruise", "cell", "decimal_year", "theta", "salinity", "aou", "y", "target"])
    return df, df.pop("target").to_numpy()


def test_differential_estimator_immune_to_cruise_drift():
    for drift in (0.0, 0.1):
        df, mask = _synthetic(drift)
        est = fit_differential(df, mask, ["y"], form="linear")["point"]["y"]
        assert abs(est - 0.7) < 0.05


def test_single_layer_estimator_absorbs_drift():
    df, mask = _synthetic(0.1)
    est = fit_trend(df[mask], ["y"], form="linear").beta["y"]
    assert abs(est - 0.9) < 0.06  # 0.8 true + 0.1 drift


def test_wet_sampler_nearest_wet_column_and_depth_interpolation():
    from research.experiments.e13_sample_model import WetSampler

    lat_g, lon_g, dep_g = np.array([0.5, 1.5]), np.array([0.5, 1.5]), np.array([0.0, 100.0, 200.0])
    wet = np.ones((3, 2, 2), bool)
    wet[2, 0, 0] = False  # column (0,0) is shallow
    field = np.arange(3)[:, None, None] * 10.0 + np.zeros((3, 2, 2))
    field[:, 1, 1] += 1000.0
    s = WetSampler(lat_g, lon_g, dep_g, wet, np.array([0.5, 0.5]), np.array([0.5, 0.5]), np.array([50.0, 150.0]))
    out = s(field, np.array([0, 1]))
    assert (s.i[0], s.j[0]) == (0, 0) and out[0] == 5.0  # shallow sample stays in its own column
    assert (s.i[1], s.j[1]) != (0, 0) and wet[2, s.i[1], s.j[1]]  # deep sample moves to a wet neighbour
    assert out[1] in (15.0, 1015.0)
