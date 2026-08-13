# Mathematical conventions

[中文](https://proffitteoy.github.io/Topp/MATHEMATICS.html) · [Repository README](https://github.com/proffitteoy/Topp/blob/main/README.en.md) · [API](API.md)

This page defines the distance semantics of the Topp `v0.1.0` public Python API. **Exact** means that the algorithm introduces no approximation parameter or approximation tolerance; inputs and results still use IEEE 754 `float64`.

## Persistence diagrams

A diagram is a multiset of points $x=(b,d)$ with $b\leq d$. Finite points strictly above the diagonal form $X_f$. Repeated rows are distinct point instances, and their multiplicity participates in matching.

A finite diagonal point $(a,a)$ has zero persistence and is ignored by distance calculations. The diagonal is

$$
\Delta=\{(t,t):t\in\mathbb{R}\}.
$$

has infinite multiplicity in a matching, so the finite parts of two diagrams may have different sizes.

## Bottleneck distance

For the finite parts, Topp uses $L^\infty$ between points and an infinitely repeated diagonal:

$$
d_B(X,Y)
=
\inf_{\gamma:X_f\cup\Delta\to Y_f\cup\Delta}
\sup_{x\in X_f\cup\Delta}
\lVert x-\gamma(x)\rVert_\infty,
$$

where $\gamma$ is a multiplicity-preserving bijection between $X_f\cup\Delta$ and $Y_f\cup\Delta$. The cost of matching a finite point $x=(b,d)$ to the diagonal is

$$
d_\infty(x,\Delta)
=
\inf_{z\in\Delta}\lVert x-z\rVert_\infty
=
\frac{d-b}{2}.
$$

`bottleneck_within(X, Y, t)` decides exactly whether `d_B(X, Y) <= t`; it does not compute an approximate distance first.

## Wasserstein distance

For $1\leq p,q<\infty$, Topp uses the convention

$$
W_{q,p}(X,Y)
=
\left(
\inf_{\gamma:X_f\cup\Delta\to Y_f\cup\Delta}
\sum_{x\in X_f\cup\Delta}
\lVert x-\gamma(x)\rVert_p^q
\right)^{1/q}.
$$

In the Python API, `order` is $q$ and `internal_p` is $p$. The matching again includes an infinitely repeated diagonal. Version `0.1.0` implements only:

- `order=1, internal_p=np.inf`, or $W_{1,\infty}$;
- `order=2, internal_p=2`, or $W_{2,2}$.

For a finite point $x=(b,d)$, the two supported diagonal costs are

$$
d_\infty(x,\Delta)=\frac{d-b}{2},
\qquad
d_2(x,\Delta)=\frac{d-b}{\sqrt{2}}.
$$

The contribution of this diagonal assignment to the squared sum in $W_{2,2}$ is therefore $(d-b)^2/2$. Other `(order, internal_p)` pairs raise `NotImplementedError`; they are never mapped to the nearest supported metric.

## Essential points

Three essential-point forms are supported:

| Type | Form | Cost between points of the same type |
|---|---|---|
| positive essential | $(b,+\infty)$ | $|b_1-b_2|$ |
| negative essential | $(-\infty,d)$ | $|d_1-d_2|$ |
| fully essential | $(-\infty,+\infty)$ | $0$ |

An essential point can match only a point of the same type, never the diagonal. If the multiplicity of any type differs between the two diagrams, both Bottleneck and Wasserstein return `inf`.

When multiplicities agree, the finite coordinates within each type are paired in sorted order. If the resulting essential costs are $e_1,\ldots,e_k$, then Bottleneck contributes $\max_i e_i$, $W_{1,\infty}$ contributes $\sum_i e_i$, and $W_{2,2}$ contributes $\sum_i e_i^2$ under the final square root. Fully essential points only require equal counts.

## Boundary behavior

- the distance between two empty diagrams is `0`;
- repeated off-diagonal points retain their multiplicity and are not treated as a set;
- diagonal points may appear in the input and in `PreparedDiagram.n_points`, but not in `n_finite_points`;
- NaN, `birth > death`, `birth=+inf`, `death=-inf`, and invalid infinity forms raise `ValueError` at the Python layer;
- a distance may be `inf`, but valid input is never silently swapped, removed, or repaired.

See the [API reference](API.md) for call signatures and exception types.
