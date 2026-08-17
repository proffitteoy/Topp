"""Exact persistence-diagram distances."""

from ._api import (
    PreparedDiagram,
    bottleneck_distance,
    bottleneck_distances,
    bottleneck_within,
    prepare_diagram,
    wasserstein_distance,
    wasserstein_distances,
)

__version__ = "1.0.0"

__all__ = [
    "PreparedDiagram",
    "bottleneck_distance",
    "bottleneck_distances",
    "bottleneck_within",
    "prepare_diagram",
    "wasserstein_distance",
    "wasserstein_distances",
]
