"""
Tests for basin definitions.
"""

import numpy as np
import pandas as pd
import pytest

import sys
sys.path.insert(0, ".")

from src.basins import assign_basin, verify_mutual_exclusivity


def _make_test_df(n=1000):
    """Create a synthetic DataFrame for testing."""
    rng = np.random.RandomState(42)
    return pd.DataFrame({
        "salinity": rng.uniform(30, 40, n),
        "temperature": rng.uniform(-2, 30, n),
        "aou": rng.uniform(0, 300, n),
        "tco2": rng.uniform(1800, 2500, n),
        "latitude": rng.uniform(-70, 70, n),
        "longitude": rng.uniform(-180, 180, n),
        "depth": rng.uniform(0, 5000, n),
        "year": rng.choice(range(1990, 2024), n),
        "region": rng.choice([1, 8, 16], n),
    })


class TestBasinAssignment:
    def test_assigns_all_rows(self):
        df = _make_test_df()
        basins = assign_basin(df)
        assert len(basins) == len(df)

    def test_mutual_exclusivity(self):
        df = _make_test_df()
        result = verify_mutual_exclusivity(df)
        assert result["mutually_exclusive"], \
            f"Overlap detected: {result['overlaps']}"

    def test_atlantic_region_code(self):
        df = pd.DataFrame({
            "latitude": [0.0, 0.0, 0.0],
            "region": [1.0, 8.0, 16.0],
        })
        basins = assign_basin(df)
        assert basins.iloc[0] == "Atlantic"
        assert basins.iloc[1] == "Pacific"
        assert basins.iloc[2] == "Indian"

    def test_southern_ocean_latitude(self):
        df = pd.DataFrame({
            "latitude": [-40.0, -30.0, -50.0],
            "region": [1.0, 1.0, 8.0],
        })
        basins = assign_basin(df)
        assert basins.iloc[0] == "Southern Ocean"  # lat -40, region 1 but south
        assert basins.iloc[1] == "Atlantic"  # lat -30, region 1, north of -35
        assert basins.iloc[2] == "Southern Ocean"  # lat -50, south

    def test_no_overlap_at_boundary(self):
        """Observations exactly at -35° should go to one basin only."""
        df = pd.DataFrame({
            "latitude": [-35.0, -34.99, -35.01],
            "region": [1.0, 1.0, 1.0],
        })
        basins = assign_basin(df)
        # At -35°: latitude < -35 is False, so goes to Atlantic (region 1, lat > -35)
        assert basins.iloc[0] == "Atlantic"
        assert basins.iloc[1] == "Atlantic"  # -34.99 > -35
        assert basins.iloc[2] == "Southern Ocean"  # -35.01 < -35


class TestTemporalSplits:
    def test_external_exclusion_and_group_disjointness(self):
        from src.splits import non_external_mask, validate_group_disjointness
        df = pd.DataFrame({"year": [2010, 2017, 2018, 2020]})
        mask = non_external_mask(df, 2018)
        assert mask.tolist() == [True, True, False, False]
        diagnostics = validate_group_disjointness(["a", "b"], ["c"])
        assert diagnostics["n_overlapping_groups"] == 0
        with pytest.raises(AssertionError):
            validate_group_disjointness(["a"], ["a"])

    def test_cruise_identifier_is_strict(self):
        from src.splits import assign_cruise_blocks
        with pytest.raises(ValueError):
            assign_cruise_blocks(pd.DataFrame({"latitude": [0], "longitude": [0]}))
        groups = assign_cruise_blocks(pd.DataFrame({"cruise": [1, 1, 2]}))
        assert groups.tolist() == ["1", "1", "2"]

    def test_no_overlap(self):
        from src.splits import split_temporal, verify_no_overlap
        df = pd.DataFrame({"year": [2010, 2015, 2018, 2020]})
        splits = split_temporal(df)
        overlap = verify_no_overlap(splits)
        assert overlap["no_overlap"] is True

    def test_correct_assignment(self):
        from src.splits import split_temporal
        df = pd.DataFrame({"year": [2010, 2014, 2015, 2017, 2018, 2020]})
        splits = split_temporal(df)
        assert splits["train"].tolist() == [True, True, False, False, False, False]
        assert splits["validation"].tolist() == [False, False, True, True, False, False]
        assert splits["external"].tolist() == [False, False, False, False, True, True]


class TestModels:
    def test_mean_baseline(self):
        from src.models.mean_baseline import MeanBaseline
        rng = np.random.RandomState(42)
        S = rng.uniform(30, 40, 100)
        T = rng.uniform(0, 20, 100)
        A = rng.uniform(0, 200, 100)
        Y = rng.uniform(2000, 2200, 100)
        model = MeanBaseline()
        model.fit(S, T, A, Y)
        pred = model.predict(S, T, A)
        assert np.allclose(pred, np.mean(Y))

    def test_linear_model(self):
        from src.models.linear import LinearModel
        rng = np.random.RandomState(42)
        S = rng.uniform(30, 40, 100)
        T = rng.uniform(0, 20, 100)
        A = rng.uniform(0, 200, 100)
        Y = 100 * S + 50 * T + 0.5 * A + rng.normal(0, 1, 100)
        model = LinearModel()
        model.fit(S, T, A, Y)
        pred = model.predict(S, T, A)
        assert len(pred) == 100
        assert model.is_fitted

    def test_rank1_quadratic(self):
        from src.models.rank1_quadratic import Rank1Quadratic
        rng = np.random.RandomState(42)
        n = 500
        S = rng.uniform(33, 37, n)
        T = rng.uniform(0, 25, n)
        A = rng.uniform(50, 250, n)
        # Generate data from known rank-1 model
        alpha, beta, gamma, delta, eps = 1.5, -0.3, 2.0, 10.0, 2000.0
        core = alpha * S + beta * T + gamma * (A / 100) + delta
        Y = core ** 2 + eps + rng.normal(0, 10, n)

        model = Rank1Quadratic()
        model.fit(S, T, A, Y)
        pred = model.predict(S, T, A)
        assert len(pred) == n
        assert model.is_fitted

        # Check that recovered coefficients are close to truth
        p = model.params_
        assert abs(p["alpha"] - alpha) < 0.1
        assert abs(p["beta"] - beta) < 0.1
        assert abs(p["gamma"] - gamma) < 0.5

    def test_rank1_hessian(self):
        from src.models.rank1_quadratic import Rank1Quadratic
        model = Rank1Quadratic()
        model.params_ = {"alpha": 1.0, "beta": 0.5, "gamma": 1.0, "delta": 0.0, "epsilon": 0.0}
        eigenvalue = model.get_hessian_eigenvalue()
        # H = 2vv^T, eigenvalue = 2||v||^2 = 2*(1^2 + 0.5^2 + (1/100)^2)
        expected = 2 * (1.0**2 + 0.5**2 + (1.0/100)**2)
        assert abs(eigenvalue - expected) < 1e-6

    def test_full_quadratic_params(self):
        from src.models.full_quadratic import FullQuadratic
        model = FullQuadratic()
        assert model.parameter_count == 10

    def test_cubic_params(self):
        from src.models.cubic import CubicPolynomial
        model = CubicPolynomial()
        assert model.parameter_count == 20


class TestFrozenSymbolic:
    def test_frozen_prediction(self):
        from src.models.frozen_symbolic import FrozenSymbolicModel
        model = FrozenSymbolicModel(
            "(sal * 2.0) + tmp - aou100",
            "2.0*sal + tmp - aou100",
            complexity=7,
        )
        S = np.array([35.0, 36.0])
        T = np.array([10.0, 15.0])
        A = np.array([100.0, 200.0])
        model.fit(S, T, A, np.array([0.0, 0.0]))
        prediction = model.predict(S, T, A)
        assert np.allclose(prediction, [79.0, 85.0])
        assert model.complexity == 7


class TestMetrics:
    def test_perfect_prediction(self):
        from src.analysis.metrics import compute_metrics
        y = np.array([1.0, 2.0, 3.0, 4.0])
        m = compute_metrics(y, y)
        assert m["RMSE"] == 0.0
        assert m["MAE"] == 0.0
        assert m["R2"] == 1.0
        assert m["Bias"] == 0.0

    def test_known_metrics(self):
        from src.analysis.metrics import compute_metrics
        y_true = np.array([1.0, 2.0, 3.0, 4.0])
        y_pred = np.array([1.1, 2.2, 2.8, 4.1])
        m = compute_metrics(y_true, y_pred)
        assert m["RMSE"] > 0
        assert m["R2"] > 0.9


class TestDerivatives:
    def test_derivative_signs(self):
        from src.analysis.derivatives import derivative_summary
        # With positive alpha, beta < 0, gamma > 0 and positive core
        # dS should be positive, dT negative, dA positive
        S = np.array([35.0, 36.0])
        T = np.array([10.0, 15.0])
        A = np.array([100.0, 200.0])
        ds = derivative_summary(S, T, A, alpha=1.5, beta=-0.3, gamma=2.0, delta=10.0)
        assert ds["frac_dS_positive_pct"] == 100.0
        assert ds["frac_dT_negative_pct"] == 100.0
        assert ds["frac_dA_positive_pct"] == 100.0


class TestEquivalence:
    def test_tost_known_equivalent(self):
        from src.validation.equivalence import tost_paired_ae
        rng = np.random.RandomState(42)
        n = 10000
        y_true = rng.normal(2000, 50, n)
        # Two nearly identical predictions (same noise level)
        y_pred_a = y_true + rng.normal(0, 2, n)
        y_pred_b = y_true + rng.normal(0, 2, n)
        result = tost_paired_ae(y_true, y_pred_a, y_pred_b, delta=5.0)
        assert result["equiv_by_test"]

    def test_tost_known_different(self):
        from src.validation.equivalence import tost_paired_ae
        rng = np.random.RandomState(42)
        n = 10000
        y_true = rng.normal(2000, 50, n)
        y_pred_a = y_true + rng.normal(0, 10, n)
        y_pred_b = y_true + 100 + rng.normal(0, 10, n)  # very different
        result = tost_paired_ae(y_true, y_pred_a, y_pred_b, delta=2.0)
        assert not result["equiv_by_test"]


class TestPareto:
    def test_pareto_frontier(self):
        from src.analysis.pareto import compute_pareto_frontier
        df = pd.DataFrame({
            "complexity": [1, 2, 3, 4, 5, 6, 7],
            "RMSE": [100, 80, 70, 75, 60, 65, 55],
        })
        pareto = compute_pareto_frontier(df)
        # Pareto-optimal: 1(100), 2(80), 3(70), 5(60), 7(55)
        assert len(pareto) == 5
        assert 4 not in pareto["complexity"].values  # 75 > 70 at complexity 3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
