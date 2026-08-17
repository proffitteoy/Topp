from __future__ import annotations

from collections.abc import Iterable
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _core

PreparedDiagram = _core.PreparedDiagram
DiagramLike: TypeAlias = ArrayLike | PreparedDiagram


def _as_diagram(diagram: ArrayLike) -> NDArray[np.float64]:
    try:
        values = np.asarray(diagram, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError("diagram must be array-like numeric data") from error
    if values.size == 0 and values.ndim == 1:
        values = values.reshape(0, 2)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("diagram must have shape (n, 2)")
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


def prepare_diagram(diagram: ArrayLike) -> PreparedDiagram:
    """Validate and precompute an immutable persistence diagram."""
    return _prepare(diagram)


def bottleneck_distance(diagram_a: DiagramLike, diagram_b: DiagramLike) -> float:
    """Return the exact L-infinity bottleneck distance."""
    return float(_core.bottleneck_distance(_prepare(diagram_a), _prepare(diagram_b)))


def _metric(order: float, internal_p: float) -> int:
    if float(order) == 1.0 and np.isposinf(internal_p):
        return 0
    if float(order) == 2.0 and float(internal_p) == 2.0:
        return 1
    raise NotImplementedError(
        "topp supports only (order=1, internal_p=inf) and "
        "(order=2, internal_p=2)"
    )


def wasserstein_distance(
    diagram_a: DiagramLike,
    diagram_b: DiagramLike,
    *,
    order: float = 1,
    internal_p: float = np.inf,
) -> float:
    """Return an exact supported Wasserstein distance."""
    return float(
        _core.wasserstein_distance(
            _prepare(diagram_a), _prepare(diagram_b), _metric(order, internal_p)
        )
    )


def _output(out: NDArray[np.float64] | None, size: int) -> NDArray[np.float64]:
    if out is None:
        return np.empty(size, dtype=np.float64)
    if not isinstance(out, np.ndarray):
        raise TypeError("out must be a NumPy ndarray")
    if out.dtype != np.dtype(np.float64):
        raise TypeError("out must have dtype float64")
    if out.shape != (size,):
        raise ValueError("out must have shape (len(diagrams),)")
    if not out.flags.c_contiguous:
        raise ValueError("out must be C-contiguous")
    if not out.flags.writeable:
        raise ValueError("out must be writeable")
    return out


def _targets(diagrams: Iterable[DiagramLike]) -> list[PreparedDiagram]:
    try:
        return [_prepare(diagram) for diagram in diagrams]
    except TypeError as error:
        raise TypeError("diagrams must be an iterable of persistence diagrams") from error


def bottleneck_distances(
    query: DiagramLike,
    diagrams: Iterable[DiagramLike],
    *,
    out: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute exact bottleneck distances from one query to many targets."""
    prepared = _targets(diagrams)
    result = _output(out, len(prepared))
    _core.bottleneck_distances(_prepare(query), prepared, result)
    return result


def wasserstein_distances(
    query: DiagramLike,
    diagrams: Iterable[DiagramLike],
    *,
    order: float = 1,
    internal_p: float = np.inf,
    out: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Compute exact supported Wasserstein distances from one query to many targets."""
    prepared = _targets(diagrams)
    result = _output(out, len(prepared))
    _core.wasserstein_distances(
        _prepare(query), prepared, _metric(order, internal_p), result
    )
    return result


def bottleneck_within(
    diagram_a: DiagramLike, diagram_b: DiagramLike, threshold: float
) -> bool:
    """Return whether the exact bottleneck distance is at most ``threshold``."""
    try:
        value = float(threshold)
    except (TypeError, ValueError) as error:
        raise TypeError("threshold must be a real number") from error
    if np.isnan(value) or value < 0:
        raise ValueError("threshold must be non-negative and not NaN")
    return bool(_core.bottleneck_within(_prepare(diagram_a), _prepare(diagram_b), value))
