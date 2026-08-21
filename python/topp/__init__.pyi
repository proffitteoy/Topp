from collections.abc import Iterable, Sequence
from decimal import Decimal
from fractions import Fraction
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._core import PreparedDiagram as PreparedDiagram

_RealScalar: TypeAlias = (
    int | float | Decimal | Fraction | np.integer[Any] | np.floating[Any]
)
DiagramLike: TypeAlias = (
    ArrayLike | Sequence[Sequence[_RealScalar]] | PreparedDiagram
)
__version__: str
__all__ = [
    "DiagramLike",
    "PreparedDiagram",
    "bottleneck_distance",
    "bottleneck_distances",
    "bottleneck_within",
    "prepare_diagram",
    "wasserstein_distance",
    "wasserstein_distances",
]

def prepare_diagram(diagram: DiagramLike) -> PreparedDiagram: ...
def bottleneck_distance(diagram_a: DiagramLike, diagram_b: DiagramLike) -> float: ...
def wasserstein_distance(
    diagram_a: DiagramLike,
    diagram_b: DiagramLike,
    *,
    order: _RealScalar = ...,
    internal_p: _RealScalar = ...,
) -> float: ...
def bottleneck_distances(
    query: DiagramLike,
    diagrams: Iterable[DiagramLike],
    *,
    out: NDArray[np.float64] | None = ...,
) -> NDArray[np.float64]: ...
def wasserstein_distances(
    query: DiagramLike,
    diagrams: Iterable[DiagramLike],
    *,
    order: _RealScalar = ...,
    internal_p: _RealScalar = ...,
    out: NDArray[np.float64] | None = ...,
) -> NDArray[np.float64]: ...
def bottleneck_within(
    diagram_a: DiagramLike, diagram_b: DiagramLike, threshold: _RealScalar
) -> bool: ...
