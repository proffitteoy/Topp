# Topp Python API

[中文](API.md) · [Usage](USAGE.en.md)

All public objects are imported from `topp`. `topp._core` is private and has no compatibility guarantee.

## Input type

`DiagramLike = ArrayLike | PreparedDiagram`. Array-like inputs must convert to a `float64` array of shape `(n, 2)`. See [Usage](USAGE.en.md#input-format) for strict validation rules.

## `PreparedDiagram`

An immutable preprocessed object returned by `prepare_diagram`.

- `n_points`: rows in the valid original input, including diagonal and essential points;
- `n_finite_points`: finite points strictly above the diagonal;
- `len(prepared)`: equal to `n_points`.

## `prepare_diagram(diagram) -> PreparedDiagram`

Validates and copies the input, then creates caches used by batch operations. Later input mutations cannot affect the prepared object.

## `bottleneck_distance(diagram_a, diagram_b) -> float`

Returns the exact Bottleneck distance with `L∞` between points and infinite diagonal multiplicity. The result may be `np.inf`.

## `wasserstein_distance(diagram_a, diagram_b, *, order=1, internal_p=np.inf) -> float`

Supported combinations:

- `order=1, internal_p=np.inf`: exact `W1-L∞`;
- `order=2, internal_p=2`: exact `W2-L2`.

Other combinations raise `NotImplementedError`.

## `bottleneck_distances(query, diagrams, *, out=None) -> np.ndarray`

Computes one-to-many exact Bottleneck distances in one native call. Returns a `float64` array of shape `(m,)`. A valid `out` is filled and returned unchanged.

## `wasserstein_distances(query, diagrams, *, order=1, internal_p=np.inf, out=None) -> np.ndarray`

Computes one-to-many Wasserstein distances in one native call and reuses a workspace across the batch. Metric parameters match the pairwise function.

## `bottleneck_within(diagram_a, diagram_b, threshold) -> bool`

Exact decision API. Returns `True` exactly when the Bottleneck distance is `<= threshold`. The threshold must be non-negative and not NaN; `+inf` is valid.

## Exceptions

- `TypeError`: non-convertible inputs, an invalid `out` type, or a non-numeric threshold;
- `ValueError`: invalid shape, point semantics, threshold, or `out` layout;
- `NotImplementedError`: unsupported Wasserstein parameters;
- `MemoryError`: native allocation failure.

## Threads and the GIL

Native distance calculations and batch loops release the Python GIL. `PreparedDiagram` is immutable and safe for concurrent reads; callers must synchronize concurrent writes to `out`.
