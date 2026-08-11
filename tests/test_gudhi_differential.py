from __future__ import annotations

import ctypes
import itertools
import math
import sys
from pathlib import Path

import gudhi
import numpy as np


def load_core(path: Path) -> ctypes._FuncPointer:
    library = ctypes.CDLL(str(path))
    function = library.bottleneck_core_distance
    function.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
    ]
    function.restype = ctypes.c_double
    # Keep the library alive for the lifetime of the function object.
    function._library = library
    return function


def as_diagram(values: list[list[float]] | np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).reshape(-1, 2)


def core_distance(
    function: ctypes._FuncPointer,
    first: np.ndarray,
    second: np.ndarray,
    config: tuple[int, int, int, int, int, int],
) -> float:
    pointer = ctypes.POINTER(ctypes.c_double)
    return float(
        function(
            first.ctypes.data_as(pointer),
            first.shape[0],
            second.ctypes.data_as(pointer),
            second.shape[0],
            *config,
        )
    )


def assert_same(actual: float, expected: float, label: str) -> None:
    if actual == expected or (math.isnan(actual) and math.isnan(expected)):
        return
    raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    library_path = (
        Path(sys.argv[1]) if len(sys.argv) > 1 else root / "build/manual/bottleneck_core_c.dll"
    )
    function = load_core(library_path)
    configs = list(itertools.product(range(4), range(6), range(5), range(10), range(10), range(2)))
    infinity = float("inf")
    cases = [
        (as_diagram([]), as_diagram([])),
        (as_diagram([]), as_diagram([[0.0, 1.0]])),
        (as_diagram([[0.0, 1.0]]), as_diagram([[0.25, 1.25]])),
        (as_diagram([[0.0, 0.0], [0.0, 1.0]]), as_diagram([[0.0, 1.0]])),
        (as_diagram([[1.0, infinity]]), as_diagram([[1.75, infinity]])),
        (as_diagram([[-infinity, 1.0]]), as_diagram([[-infinity, 2.0]])),
        (as_diagram([[0.0, 1.0], [0.0, 1.0]]), as_diagram([[0.0, 1.0]])),
    ]

    generator = np.random.default_rng(0xB0771E)
    for _ in range(500):
        first_size = int(generator.integers(0, 25))
        second_size = int(generator.integers(0, 25))
        first_birth = generator.uniform(-2.0, 2.0, size=first_size)
        second_birth = generator.uniform(-2.0, 2.0, size=second_size)
        first = np.column_stack(
            (first_birth, first_birth + generator.uniform(0.0, 2.0, size=first_size))
        )
        second = np.column_stack(
            (second_birth, second_birth + generator.uniform(0.0, 2.0, size=second_size))
        )
        cases.append((as_diagram(first), as_diagram(second)))

    comparisons = 0
    for case_index, (first, second) in enumerate(cases):
        expected = float(gudhi.bottleneck_distance(first, second, e=0.0))
        for config in configs:
            actual = core_distance(function, first, second, config)
            assert_same(actual, expected, f"case={case_index}, config={config}")
            comparisons += 1

    print(
        f"gudhi differential: {len(cases)} cases, {len(configs)} configs, "
        f"{comparisons} exact comparisons passed"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
