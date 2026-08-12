from collections.abc import Iterable
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

class PreparedDiagram:
    @property
    def n_points(self) -> int: ...
    @property
    def n_finite_points(self) -> int: ...
    def __len__(self) -> int: ...

DiagramLike: TypeAlias = ArrayLike | PreparedDiagram
__version__: str

def prepare_diagram(diagram: ArrayLike) -> PreparedDiagram: ...
def bottleneck_distance(diagram_a: DiagramLike, diagram_b: DiagramLike) -> float: ...
def wasserstein_distance(
    diagram_a: DiagramLike,
    diagram_b: DiagramLike,
    *,
    order: float = ...,
    internal_p: float = ...,
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
    order: float = ...,
    internal_p: float = ...,
    out: NDArray[np.float64] | None = ...,
) -> NDArray[np.float64]: ...
def bottleneck_within(
    diagram_a: DiagramLike, diagram_b: DiagramLike, threshold: float
) -> bool: ...
