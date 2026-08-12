from __future__ import annotations

import itertools
import math
import os

import numpy as np
import pytest

import topp


def _random_diagram(generator: np.random.Generator, size: int) -> np.ndarray:
    births = generator.uniform(-2.0, 2.0, size=size)
    return np.column_stack((births, births + generator.uniform(0.0, 2.0, size=size)))


@pytest.mark.oracle
@pytest.mark.skipif(os.environ.get("TOPP_RUN_ORACLE_TESTS") != "1", reason="oracle suite is opt-in")
def test_bottleneck_matches_gudhi_313() -> None:
    gudhi = pytest.importorskip("gudhi")
    generator = np.random.default_rng(0x70_50_50)
    for _ in range(100):
        first = _random_diagram(generator, int(generator.integers(0, 25)))
        second = _random_diagram(generator, int(generator.integers(0, 25)))
        expected = float(gudhi.bottleneck_distance(first, second, e=0.0))
        assert topp.bottleneck_distance(first, second) == expected


@pytest.mark.oracle
@pytest.mark.skipif(os.environ.get("TOPP_RUN_ORACLE_TESTS") != "1", reason="oracle suite is opt-in")
@pytest.mark.parametrize("order,internal_p", [(1, math.inf), (2, 2)])
def test_wasserstein_matches_gudhi_pot(
    order: int, internal_p: float
) -> None:
    pytest.importorskip("ot")
    gudhi_wasserstein = pytest.importorskip("gudhi.wasserstein")
    generator = np.random.default_rng(0x50_07)
    for _ in range(50):
        first = _random_diagram(generator, int(generator.integers(0, 12)))
        second = _random_diagram(generator, int(generator.integers(0, 12)))
        expected = float(
            gudhi_wasserstein.wasserstein_distance(
                first,
                second,
                order=order,
                internal_p=internal_p,
            )
        )
        actual = topp.wasserstein_distance(
            first, second, order=order, internal_p=internal_p
        )
        assert actual == pytest.approx(expected, rel=2e-12, abs=2e-12)


def _diagonal_power(point: np.ndarray, metric: int) -> float:
    persistence = float(point[1] - point[0])
    return persistence / 2 if metric == 1 else persistence * persistence / 2


def _cross_power(first: np.ndarray, second: np.ndarray, metric: int) -> float:
    delta = np.abs(first - second)
    return float(delta.max()) if metric == 1 else float(np.dot(delta, delta))


def _brute_force(first: np.ndarray, second: np.ndarray, metric: int) -> float:
    baseline = sum(_diagonal_power(point, metric) for point in first)
    baseline += sum(_diagonal_power(point, metric) for point in second)
    best_saving = 0.0
    for count in range(1, min(len(first), len(second)) + 1):
        for rows in itertools.combinations(range(len(first)), count):
            for columns in itertools.permutations(range(len(second)), count):
                saving = sum(
                    _diagonal_power(first[row], metric)
                    + _diagonal_power(second[column], metric)
                    - _cross_power(first[row], second[column], metric)
                    for row, column in zip(rows, columns, strict=True)
                )
                best_saving = max(best_saving, saving)
    powered = max(0.0, baseline - best_saving)
    return powered if metric == 1 else math.sqrt(powered)


@pytest.mark.parametrize("metric", [1, 2])
def test_wasserstein_matches_independent_bruteforce(metric: int) -> None:
    generator = np.random.default_rng(0x57_A55)
    kwargs = {} if metric == 1 else {"order": 2, "internal_p": 2}
    for _ in range(60):
        first = _random_diagram(generator, int(generator.integers(0, 5)))
        second = _random_diagram(generator, int(generator.integers(0, 5)))
        actual = topp.wasserstein_distance(first, second, **kwargs)
        expected = _brute_force(first, second, metric)
        assert actual == pytest.approx(expected, rel=2e-12, abs=2e-12)
