from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
import topp


def prepare(value: topp.DiagramLike) -> topp.PreparedDiagram:
    return topp.prepare_diagram(value)


diagram: list[list[Decimal | Fraction]] = [
    [Decimal("0"), Fraction(1, 1)]
]
prepared = prepare(diagram)
same = topp.prepare_diagram(prepared)
distance: float = topp.bottleneck_distance(prepared, same)
batch: NDArray[np.float64] = topp.bottleneck_distances(prepared, [same])
within: bool = topp.bottleneck_within(prepared, same, 0.0)

assert distance == 0.0
assert batch.shape == (1,)
assert within

if TYPE_CHECKING:
    topp.PreparedDiagram()  # type: ignore[call-arg]
