from __future__ import annotations

import concurrent.futures
from decimal import Decimal
from fractions import Fraction
import math

import numpy as np
import pytest

import topp


X = np.array([[0.0, 1.0], [0.25, 0.75]], dtype=np.float64)
Y = np.array([[0.0, 1.1], [0.3, 0.8]], dtype=np.float64)


def test_public_surface_and_version() -> None:
    assert topp.__version__ == "1.0.1"
    assert set(topp.__all__) == {
        "DiagramLike",
        "PreparedDiagram",
        "prepare_diagram",
        "bottleneck_distance",
        "wasserstein_distance",
        "bottleneck_distances",
        "wasserstein_distances",
        "bottleneck_within",
    }
    assert topp.DiagramLike is not None


def test_array_like_dtype_and_non_contiguous_inputs() -> None:
    base = np.array([[0, 99, 1], [1, 99, 3]], dtype=np.int32)
    view = base[:, ::2]
    assert not view.flags.c_contiguous
    expected = topp.bottleneck_distance([[0, 1], [1, 3]], X)
    assert topp.bottleneck_distance(view, X) == expected


def test_empty_and_prepared_metadata() -> None:
    empty = topp.prepare_diagram([])
    assert topp.prepare_diagram(empty) is empty
    assert len(empty) == 0
    assert empty.n_points == 0
    assert empty.n_finite_points == 0
    assert topp.bottleneck_distance(empty, np.empty((0, 2))) == 0.0

    prepared = topp.prepare_diagram([[0, 0], [0, 1], [2, 2], [3, math.inf]])
    assert len(prepared) == 4
    assert prepared.n_points == 4
    assert prepared.n_finite_points == 1
    with pytest.raises(TypeError):
        topp.PreparedDiagram()


def test_decimal_and_fraction_inputs() -> None:
    diagram = np.array(
        [
            [Decimal("0.0"), Decimal("1.0")],
            [Fraction(1, 4), Fraction(3, 4)],
        ],
        dtype=object,
    )
    assert topp.bottleneck_distance(diagram, diagram) == 0.0


def test_decimal_essential_points() -> None:
    diagram = [
        [Decimal("1"), Decimal("Infinity")],
        [Decimal("-Infinity"), Decimal("2")],
        [Decimal("-Infinity"), Decimal("Infinity")],
    ]
    assert topp.bottleneck_distance(diagram, diagram) == 0.0


def test_prepared_diagram_owns_its_input_copy() -> None:
    source = X.copy()
    expected = source.copy()
    prepared = topp.prepare_diagram(source)
    source[:] = 0.0
    assert topp.bottleneck_distance(prepared, expected) == 0.0


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
    for diagram in (
        [0, 1],
        [[0, 1, 2]],
        [[0, 1], [2]],
        np.zeros((2, 2, 1)),
    ):
        with pytest.raises(ValueError):
            topp.prepare_diagram(diagram)
    with pytest.raises(TypeError):
        topp.prepare_diagram("hello")

    class BrokenArray:
        def __array__(self, *args: object, **kwargs: object) -> np.ndarray:
            raise TypeError("array conversion failed")

    with pytest.raises(TypeError, match="array-like numeric data"):
        topp.prepare_diagram(BrokenArray())  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "diagram",
    [
        np.array([[0.0 + 0.0j, 1.0 + 1.0j]]),
        [["0", "1"]],
        [[False, True]],
        np.array([[0, 1]], dtype="datetime64[D]"),
        np.ma.array([[0.0, 1.0]], mask=[[False, False]]),
        np.ma.array([[0.0, 1.0]], mask=[[False, True]]),
        [[0.0, object()]],
    ],
)
def test_non_real_diagram_values_raise_type_error(diagram: object) -> None:
    with pytest.raises(TypeError):
        topp.prepare_diagram(diagram)


@pytest.mark.parametrize(
    "diagram",
    [
        [[0, 10**400]],
        [[Decimal("0"), Decimal("1e10000")]],
        [[Fraction(0), Fraction(10**400, 1)]],
    ],
)
def test_finite_values_must_be_representable_as_float64(diagram: object) -> None:
    with pytest.raises(ValueError, match="representable as float64"):
        topp.prepare_diagram(diagram)


def test_wider_numpy_float_overflow_is_not_an_essential_point() -> None:
    longdouble_max = np.finfo(np.longdouble).max
    if longdouble_max <= np.finfo(np.float64).max:
        pytest.skip("platform longdouble has no wider finite range than float64")
    with pytest.raises(ValueError, match="representable as float64"):
        topp.prepare_diagram(np.array([[np.longdouble(0), longdouble_max]]))


def test_wasserstein_metrics_and_unsupported_combinations() -> None:
    w1 = topp.wasserstein_distance(X, Y)
    w2 = topp.wasserstein_distance(X, Y, order=2, internal_p=2)
    assert w1 >= 0
    assert w2 >= 0
    with pytest.raises(NotImplementedError):
        topp.wasserstein_distance(X, Y, order=3, internal_p=2)


def test_decimal_and_fraction_scalar_parameters() -> None:
    expected = topp.wasserstein_distance(X, Y, order=2, internal_p=2)
    actual = topp.wasserstein_distance(
        X, Y, order=Decimal("2"), internal_p=Fraction(2, 1)
    )
    assert actual == expected
    assert topp.bottleneck_within(X, Y, Decimal("Infinity"))
    assert topp.bottleneck_within(X, X, Fraction(0, 1))
    assert topp.wasserstein_distance(
        X, Y, order=np.int64(2), internal_p=np.float32(2)
    ) == expected


def test_wasserstein_parameters_are_validated_before_diagrams() -> None:
    with pytest.raises(NotImplementedError):
        topp.wasserstein_distance([0, 1], Y, order=3, internal_p=2)
    with pytest.raises(ValueError, match="positive and finite"):
        topp.wasserstein_distance([0, 1], Y, order=0)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"order": True},
        {"order": "1"},
        {"order": 1 + 0j},
        {"order": np.array(1.0)},
        {"internal_p": True},
        {"internal_p": "inf"},
        {"internal_p": 2 + 0j},
    ],
)
def test_wasserstein_parameter_types(kwargs: dict[str, object]) -> None:
    with pytest.raises(TypeError):
        topp.wasserstein_distance(X, Y, **kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"order": 0},
        {"order": math.nan},
        {"order": math.inf},
        {"order": 10**400},
        {"internal_p": 0},
        {"internal_p": math.nan},
        {"internal_p": -math.inf},
        {"internal_p": Decimal("1e10000")},
    ],
)
def test_wasserstein_parameter_values(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        topp.wasserstein_distance(X, Y, **kwargs)


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

    with pytest.raises(ValueError, match="one-dimensional"):
        topp.bottleneck_distances(X, [X], out=np.empty((1, 1), dtype=np.float64))

    buffer = bytearray(17)
    unaligned = np.ndarray((2,), dtype=np.float64, buffer=buffer, offset=1)
    assert unaligned.flags.c_contiguous and not unaligned.flags.aligned
    with pytest.raises(ValueError, match="aligned"):
        topp.bottleneck_distances(X, [X, Y], out=unaligned)


def test_native_batch_rejects_unaligned_output() -> None:
    from topp import _core

    query = topp.prepare_diagram(X)
    buffer = bytearray(17)
    unaligned = np.ndarray((2,), dtype=np.float64, buffer=buffer, offset=1)
    with pytest.raises(ValueError, match="aligned"):
        _core.bottleneck_distances(query, [query, query], unaligned)


def test_batch_validation_does_not_consume_targets_early() -> None:
    seen: list[int] = []

    def targets():
        for index in range(3):
            seen.append(index)
            yield X

    with pytest.raises(NotImplementedError):
        topp.wasserstein_distances(X, targets(), order=3, internal_p=2)
    assert seen == []

    with pytest.raises(ValueError, match="shape"):
        topp.bottleneck_distances([0, 1], targets())
    assert seen == []

    with pytest.raises(TypeError, match="dtype float64"):
        topp.bottleneck_distances(X, targets(), out=np.empty(3, dtype=np.float32))
    assert seen == []


def test_each_invalid_output_layout_is_rejected_before_targets() -> None:
    readonly = np.empty(2, dtype=np.float64)
    readonly.flags.writeable = False
    unaligned_buffer = bytearray(17)
    invalid_outputs = [
        [0.0, 0.0],
        np.empty(2, dtype=np.float32),
        np.empty((1, 2), dtype=np.float64),
        np.empty(4, dtype=np.float64)[::2],
        readonly,
        np.ndarray((2,), dtype=np.float64, buffer=unaligned_buffer, offset=1),
    ]

    for out in invalid_outputs:
        seen: list[int] = []

        def targets():
            seen.append(0)
            yield X

        with pytest.raises((TypeError, ValueError)):
            topp.bottleneck_distances(X, targets(), out=out)  # type: ignore[arg-type]
        assert seen == []


def test_batch_target_errors_include_index() -> None:
    with pytest.raises(TypeError, match=r"diagrams\[1\]") as caught:
        topp.bottleneck_distances(X, [X, [[0.0, object()]]])
    assert isinstance(caught.value.__cause__, TypeError)

    with pytest.raises(ValueError, match=r"diagrams\[1\].*shape") as caught_value:
        topp.bottleneck_distances(X, [X, [0.0, 1.0]])
    assert isinstance(caught_value.value.__cause__, ValueError)


def test_batch_requires_an_iterable_of_targets() -> None:
    with pytest.raises(TypeError, match="iterable of persistence diagrams"):
        topp.bottleneck_distances(X, 42)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "batch",
    [topp.bottleneck_distances, topp.wasserstein_distances],
)
def test_batch_consumes_each_generator_target_once(batch: object) -> None:
    seen: list[int] = []

    def targets():
        for index, target in enumerate((X, Y)):
            seen.append(index)
            yield target

    result = batch(X, targets())  # type: ignore[operator]
    assert result.shape == (2,)
    assert seen == [0, 1]


def test_generator_errors_are_not_relabelled_as_target_errors() -> None:
    def targets():
        yield X
        raise TypeError("generator exploded")

    with pytest.raises(TypeError, match="generator exploded"):
        topp.bottleneck_distances(X, targets())


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


@pytest.mark.parametrize("threshold", [True, "0.1", np.array(0.1)])
def test_bottleneck_within_rejects_non_scalar_real_types(threshold: object) -> None:
    with pytest.raises(TypeError, match="real number"):
        topp.bottleneck_within(X, Y, threshold)


@pytest.mark.parametrize(
    "threshold", [10**400, Decimal("1e10000"), Fraction(10**400, 1)]
)
def test_bottleneck_within_rejects_float64_overflow(threshold: object) -> None:
    with pytest.raises(ValueError, match="representable as float64"):
        topp.bottleneck_within(X, Y, threshold)


def test_prepared_is_safe_for_concurrent_reads() -> None:
    prepared = topp.prepare_diagram(X)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results = list(
            pool.map(lambda _: topp.bottleneck_distance(prepared, Y), range(16))
        )
    assert results == [results[0]] * len(results)
