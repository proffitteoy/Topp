from __future__ import annotations

import concurrent.futures
import math

import numpy as np
import pytest

import topp


X = np.array([[0.0, 1.0], [0.25, 0.75]], dtype=np.float64)
Y = np.array([[0.0, 1.1], [0.3, 0.8]], dtype=np.float64)


def test_public_surface_and_version() -> None:
    assert topp.__version__ == "1.0.0"
    assert set(topp.__all__) == {
        "PreparedDiagram",
        "prepare_diagram",
        "bottleneck_distance",
        "wasserstein_distance",
        "bottleneck_distances",
        "wasserstein_distances",
        "bottleneck_within",
    }


def test_array_like_dtype_and_non_contiguous_inputs() -> None:
    base = np.array([[0, 99, 1], [1, 99, 3]], dtype=np.int32)
    view = base[:, ::2]
    assert not view.flags.c_contiguous
    expected = topp.bottleneck_distance([[0, 1], [1, 3]], X)
    assert topp.bottleneck_distance(view, X) == expected


def test_empty_and_prepared_metadata() -> None:
    empty = topp.prepare_diagram([])
    assert len(empty) == 0
    assert empty.n_points == 0
    assert empty.n_finite_points == 0
    assert topp.bottleneck_distance(empty, np.empty((0, 2))) == 0.0

    prepared = topp.prepare_diagram([[0, 0], [0, 1], [2, 2], [3, math.inf]])
    assert len(prepared) == 4
    assert prepared.n_points == 4
    assert prepared.n_finite_points == 1


def test_duplicates_and_essential_points() -> None:
    duplicates = [[0, 1], [0, 1], [0, 1]]
    assert topp.bottleneck_distance(duplicates, duplicates) == 0.0
    essential_a = [[1, math.inf], [-math.inf, 2], [-math.inf, math.inf]]
    essential_b = [[2, math.inf], [-math.inf, 4], [-math.inf, math.inf]]
    assert topp.bottleneck_distance(essential_a, essential_b) == 2.0
    assert topp.wasserstein_distance(essential_a, essential_b) == 3.0


@pytest.mark.parametrize(
    "diagram",
    [
        [[0, math.nan]],
        [[2, 1]],
        [[math.inf, math.inf]],
        [[math.inf, 2]],
        [[0, -math.inf]],
    ],
)
def test_invalid_points_raise_value_error(diagram: object) -> None:
    with pytest.raises(ValueError):
        topp.prepare_diagram(diagram)


def test_invalid_shapes_and_types() -> None:
    for diagram in ([0, 1], [[0, 1, 2]], np.zeros((2, 2, 1))):
        with pytest.raises(ValueError):
            topp.prepare_diagram(diagram)
    with pytest.raises(TypeError):
        topp.prepare_diagram("hello")


def test_wasserstein_metrics_and_unsupported_combinations() -> None:
    w1 = topp.wasserstein_distance(X, Y)
    w2 = topp.wasserstein_distance(X, Y, order=2, internal_p=2)
    assert w1 >= 0
    assert w2 >= 0
    with pytest.raises(NotImplementedError):
        topp.wasserstein_distance(X, Y, order=3, internal_p=2)


def test_prepared_and_batch_match_single_calls() -> None:
    query = topp.prepare_diagram(X)
    targets = [Y, topp.prepare_diagram([[0.0, 2.0]]), []]
    bottleneck = topp.bottleneck_distances(query, targets)
    expected_b = np.array([topp.bottleneck_distance(query, target) for target in targets])
    np.testing.assert_array_equal(bottleneck, expected_b)

    wasserstein = topp.wasserstein_distances(query, targets, order=2, internal_p=2)
    expected_w = np.array(
        [
            topp.wasserstein_distance(query, target, order=2, internal_p=2)
            for target in targets
        ]
    )
    np.testing.assert_allclose(wasserstein, expected_w, rtol=1e-13, atol=1e-13)


def test_batch_out_contract() -> None:
    out = np.empty(2, dtype=np.float64)
    returned = topp.bottleneck_distances(X, [X, Y], out=out)
    assert returned is out
    np.testing.assert_array_equal(out, [0.0, topp.bottleneck_distance(X, Y)])

    with pytest.raises(TypeError):
        topp.bottleneck_distances(X, [X], out=np.empty(1, dtype=np.float32))
    with pytest.raises(ValueError):
        topp.bottleneck_distances(X, [X], out=np.empty(2, dtype=np.float64))
    non_contiguous = np.empty(4, dtype=np.float64)[::2]
    with pytest.raises(ValueError):
        topp.bottleneck_distances(X, [X, Y], out=non_contiguous)
    readonly = np.empty(1, dtype=np.float64)
    readonly.flags.writeable = False
    with pytest.raises(ValueError):
        topp.bottleneck_distances(X, [X], out=readonly)


def test_empty_batch() -> None:
    result = topp.wasserstein_distances(X, [])
    assert result.dtype == np.float64
    assert result.shape == (0,)


def test_bottleneck_within_contract() -> None:
    distance = topp.bottleneck_distance(X, Y)
    assert topp.bottleneck_within(X, Y, distance)
    assert topp.bottleneck_within(X, Y, math.inf)
    if distance > 0:
        assert not topp.bottleneck_within(X, Y, np.nextafter(distance, -math.inf))
    for threshold in (-1, math.nan):
        with pytest.raises(ValueError):
            topp.bottleneck_within(X, Y, threshold)


def test_prepared_is_safe_for_concurrent_reads() -> None:
    prepared = topp.prepare_diagram(X)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = list(
            pool.map(lambda _: topp.bottleneck_distance(prepared, Y), range(16))
        )
    assert results == [results[0]] * len(results)
