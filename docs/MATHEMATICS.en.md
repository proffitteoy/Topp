# Mathematical conventions

[中文](MATHEMATICS.md) · [README](../README.en.md) · [API](API.en.md)

This page defines the distance semantics of the Topp `v0.1.0` public Python API. **Exact** means that the algorithm introduces no approximation parameter or approximation tolerance; inputs and results still use IEEE 754 `float64`.

## Persistence diagrams

A diagram is a multiset of points `x = (birth, death)` with `birth <= death`. Finite points strictly above the diagonal form `X_f`. Repeated rows are distinct point instances, and their multiplicity participates in matching.

A finite diagonal point `(a, a)` has zero persistence and is ignored by distance calculations. The diagonal

```text
Δ = {(t, t) | t ∈ R}
```

has infinite multiplicity in a matching, so the finite parts of two diagrams may have different sizes.

## Bottleneck distance

For the finite parts, Topp uses `L∞` between points and an infinitely repeated diagonal:

```text
d_B(X, Y) = inf_γ sup_x ||x - γ(x)||∞,
```

where `γ` is a multiplicity-preserving bijection between `X_f ∪ Δ` and `Y_f ∪ Δ`. The cost of matching a finite point `x = (b, d)` to the diagonal is

```text
dist∞(x, Δ) = (d - b) / 2.
```

`bottleneck_within(X, Y, t)` decides exactly whether `d_B(X, Y) <= t`; it does not compute an approximate distance first.

## Wasserstein distance

Topp uses the following `W_{q,p}` convention:

```text
W_{q,p}(X, Y) = (inf_γ Σ_x ||x - γ(x)||p^q)^(1/q),
```

where `order = q`, `internal_p = p`, and the matching again includes an infinitely repeated diagonal. Version `0.1.0` implements only:

- `order=1, internal_p=np.inf`, or `W_{1,∞}`;
- `order=2, internal_p=2`, or `W_{2,2}`.

For a finite point `x = (b, d)`:

```text
dist∞(x, Δ) = (d - b) / 2
dist2(x, Δ) = (d - b) / √2
```

The contribution of this diagonal assignment to the squared sum in `W_{2,2}` is therefore `(d - b)^2 / 2`. Other `(order, internal_p)` pairs raise `NotImplementedError`; they are never mapped to the nearest supported metric.

## Essential points

Three essential-point forms are supported:

| Type | Form | Cost between points of the same type |
|---|---|---|
| positive essential | `(b, +inf)` | `|b₁ - b₂|` |
| negative essential | `(-inf, d)` | `|d₁ - d₂|` |
| fully essential | `(-inf, +inf)` | `0` |

An essential point can match only a point of the same type, never the diagonal. If the multiplicity of any type differs between the two diagrams, both Bottleneck and Wasserstein return `inf`.

When multiplicities agree, the finite coordinates within each type are paired in sorted order. Bottleneck takes the maximum of the resulting essential costs and the finite-part distance. `W_{1,∞}` sums the essential costs; `W_{2,2}` adds their squares to the total before taking the square root. Fully essential points only require equal counts.

## Boundary behavior

- the distance between two empty diagrams is `0`;
- repeated off-diagonal points retain their multiplicity and are not treated as a set;
- diagonal points may appear in the input and in `PreparedDiagram.n_points`, but not in `n_finite_points`;
- NaN, `birth > death`, `birth=+inf`, `death=-inf`, and invalid infinity forms raise `ValueError` at the Python layer;
- a distance may be `inf`, but valid input is never silently swapped, removed, or repaired.

See the [API reference](API.en.md) for call signatures and exception types.
