from __future__ import annotations

from collections.abc import Iterable, Sequence
from decimal import Decimal
from fractions import Fraction
import math
from numbers import Real
from typing import TypeAlias
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _core
from ._core import PreparedDiagram

_RealScalar: TypeAlias = (
    int | float | Decimal | Fraction | np.integer | np.floating
)
_DiagramArrayLike: TypeAlias = ArrayLike | Sequence[Sequence[_RealScalar]]
DiagramLike: TypeAlias = _DiagramArrayLike | PreparedDiagram


def _is_real_scalar(value: object) -> bool:
    return not isinstance(value, (bool, np.bool_)) and (
        isinstance(value, Real) or isinstance(value, Decimal)
    )


def _source_is_finite(value: object) -> bool:
    if isinstance(value, Decimal):
        return value.is_finite()
    try:
        return math.isfinite(value)  # type: ignore[arg-type]
    except OverflowError:
        # Arbitrary-precision integers and fractions are finite even when they
        # are too large for Python's float conversion.
        return True


def _validate_real_array(values: NDArray[np.generic]) -> None:
    if values.dtype.kind in "iuf":
        return
    if values.dtype.kind == "O":
        if all(_is_real_scalar(value) for value in values.flat):
            return
    raise TypeError("diagram must contain only real numeric values")


def _reject_finite_overflow(
    source: NDArray[np.generic], converted: NDArray[np.float64]
) -> None:
    infinite = np.isinf(converted)
    if not infinite.any():
        return
    if source.dtype.kind == "f":
        if np.any(np.isfinite(source) & infinite):
            raise ValueError("diagram values must be representable as float64")
        return
    if source.dtype.kind == "O":
        for original, became_infinite in zip(
            source.flat, infinite.flat, strict=True
        ):
            if became_infinite and _source_is_finite(original):
                raise ValueError("diagram values must be representable as float64")


def _as_diagram(diagram: _DiagramArrayLike) -> NDArray[np.float64]:
    if np.ma.isMaskedArray(diagram):
        raise TypeError("diagram must not be a masked array")
    try:
        source = np.asarray(diagram)
    except TypeError as error:
        raise TypeError("diagram must be array-like numeric data") from error
    except ValueError as error:
        raise ValueError("diagram must have shape (n, 2)") from error
    _validate_real_array(source)
    if source.size == 0 and source.ndim == 1:
        source = source.reshape(0, 2)
    if source.ndim != 2 or source.shape[1] != 2:
        raise ValueError("diagram must have shape (n, 2)")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            values = np.asarray(source, dtype=np.float64)
    except (OverflowError, RuntimeWarning, TypeError, ValueError) as error:
        raise ValueError("diagram values must be representable as float64") from error
    _reject_finite_overflow(source, values)
    if np.isnan(values).any():
        raise ValueError("diagram must not contain NaN")

    births = values[:, 0]
    deaths = values[:, 1]
    if np.isposinf(births).any():
        raise ValueError("birth=+inf is not a supported essential-point form")
    if np.isneginf(deaths).any():
        raise ValueError("death=-inf is not a supported essential-point form")
    if np.any(births > deaths):
        raise ValueError("diagram points must satisfy birth <= death")
    return np.ascontiguousarray(values, dtype=np.float64)


def _prepare(diagram: DiagramLike) -> PreparedDiagram:
    if isinstance(diagram, PreparedDiagram):
        return diagram
    return _core.prepare_diagram(_as_diagram(diagram))


def prepare_diagram(diagram: DiagramLike) -> PreparedDiagram:
    """Validate and precompute an immutable persistence diagram."""
    return _prepare(diagram)


def bottleneck_distance(diagram_a: DiagramLike, diagram_b: DiagramLike) -> float:
    """Return the exact L-infinity bottleneck distance."""
    return float(_core.bottleneck_distance(_prepare(diagram_a), _prepare(diagram_b)))


def _metric(order: _RealScalar, internal_p: _RealScalar) -> int:
    order_value = _real_parameter(order, "order")
    internal_p_value = _real_parameter(internal_p, "internal_p")
    if not np.isfinite(order_value) or order_value <= 0:
        raise ValueError("order must be positive and finite")
    if np.isnan(internal_p_value) or internal_p_value <= 0:
        raise ValueError("internal_p must be positive and not NaN")
    if order_value == 1.0 and np.isposinf(internal_p_value):
        return 0
    if order_value == 2.0 and internal_p_value == 2.0:
        return 1
    raise NotImplementedError(
        "topp supports only (order=1, internal_p=inf) and "
        "(order=2, internal_p=2)"
    )


def wasserstein_distance(
    diagram_a: DiagramLike,
    diagram_b: DiagramLike,
    *,
    order: _RealScalar = 1,
    internal_p: _RealScalar = np.inf,
) -> float:
    """Return an exact supported Wasserstein distance."""
    metric = _metric(order, internal_p)
    return float(
        _core.wasserstein_distance(
            _prepare(diagram_a), _prepare(diagram_b), metric
        )
    )


def _validate_output_layout(out: NDArray[np.float64] | None) -> None:
    if out is None:
        return
    if not isinstance(out, np.ndarray):
        raise TypeError("out must be a NumPy ndarray")
    if out.dtype != np.dtype(np.float64):
        raise TypeError("out must have dtype float64")
    if out.ndim != 1:
        raise ValueError("out must be one-dimensional")
    if not out.flags.c_contiguous:
        raise ValueError("out must be C-contiguous")
    if not out.flags.aligned:
        raise ValueError("out must be aligned")
    if not out.flags.writeable:
        raise ValueError("out must be writeable")


def _output(out: NDArray[np.float64] | None, size: int) -> NDArray[np.float64]:
    if out is None:
        return np.empty(size, dtype=np.float64)
    _validate_output_layout(out)
    if out.shape != (size,):
        raise ValueError("out must have shape (len(diagrams),)")
    return out


def _targets(diagrams: Iterable[DiagramLike]) -> list[PreparedDiagram]:
    try:
        iterator = iter(diagrams)
    except TypeError as error:
        raise TypeError("diagrams must be an iterable of persistence diagrams") from error
    prepared: list[PreparedDiagram] = []
    for index, diagram in enumerate(iterator):
        try:
            prepared.append(_prepare(diagram))
        except (TypeError, ValueError) as error:
            raise type(error)(f"diagrams[{index}]: {error}") from error
    return prepared


def bottleneck_distances(
    query: DiagramLike,
    diagrams: Iterable[DiagramLike],
    *,
    out: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute exact bottleneck distances from one query to many targets."""
    prepared_query = _prepare(query)
    _validate_output_layout(out)
    prepared = _targets(diagrams)
    result = _output(out, len(prepared))
    _core.bottleneck_distances(prepared_query, prepared, result)
    return result


def wasserstein_distances(
    query: DiagramLike,
    diagrams: Iterable[DiagramLike],
    *,
    order: _RealScalar = 1,
    internal_p: _RealScalar = np.inf,
    out: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute exact supported Wasserstein distances from one query to many targets."""
    metric = _metric(order, internal_p)
    prepared_query = _prepare(query)
    _validate_output_layout(out)
    prepared = _targets(diagrams)
    result = _output(out, len(prepared))
    _core.wasserstein_distances(prepared_query, prepared, metric, result)
    return result


def _real_parameter(value: _RealScalar, name: str) -> float:
    if not _is_real_scalar(value):
        raise TypeError(f"{name} must be a real number")
    try:
        converted = float(value)
    except (OverflowError, TypeError, ValueError) as error:
        raise ValueError(f"{name} must be representable as float64") from error
    if np.isinf(converted) and _source_is_finite(value):
        raise ValueError(f"{name} must be representable as float64")
    return converted


def bottleneck_within(
    diagram_a: DiagramLike, diagram_b: DiagramLike, threshold: _RealScalar
) -> bool:
    """Return whether the exact bottleneck distance is at most ``threshold``."""
    value = _real_parameter(threshold, "threshold")
    if np.isnan(value) or value < 0:
        raise ValueError("threshold must be non-negative and not NaN")
    return bool(_core.bottleneck_within(_prepare(diagram_a), _prepare(diagram_b), value))
