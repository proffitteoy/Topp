"""Exact persistence-diagram distances."""

from ._api import (
    DiagramLike,
    PreparedDiagram,
    bottleneck_distance,
    bottleneck_distances,
    bottleneck_within,
    prepare_diagram,
    wasserstein_distance,
    wasserstein_distances,
)

__version__ = "1.0.1"

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
