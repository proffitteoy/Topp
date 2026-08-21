# Topp v1.0.1

Topp 1.0.1 tightens the stable Python input and error contracts without changing
the exact C++ distance kernels.

## Highlights

- Strict real-valued diagram inputs with explicit rejection of complex,
  masked, boolean, string, and other non-real data.
- Consistent `TypeError`, `ValueError`, and `NotImplementedError` categories for
  diagram, scalar-parameter, and Wasserstein-mode failures.
- Indexed batch-target errors and validation ordering that does not consume a
  generator when the metric, query, or output layout is already invalid.
- Runtime and typing support for the public `DiagramLike` alias.
- Aligned caller-output enforcement in the Python and native binding layers.
- Linux x86_64 manylinux wheels plus Linux C++/Python, sdist, and TestPyPI
  installation coverage for CPython 3.10–3.14.

## Compatibility

Integer and floating-point arrays, Python numeric sequences, `Decimal`,
`Fraction`, empty diagrams, documented essential points, prepared diagrams, and
all supported exact distance modes remain supported. Callers that previously
relied on coercing boolean or numeric-string diagrams must convert them to real
numeric values before calling Topp.

## Validation target

- Windows x64 and Linux x86_64 wheels for CPython 3.10–3.14 and an sdist.
- Python API regressions, strict mypy checks, and installed-package stub checks.
- Existing C++ exact Bottleneck and Wasserstein regressions.
- Strict Chinese and English Sphinx builds and package smoke tests.

## Known limits

- macOS remains outside the release-CI platform matrix.
- Linux wheels use the portable scalar kernel; Windows wheels may use the
  runtime-dispatched AVX2 implementation.
- Wasserstein supports only exact `W1-L∞` and `W2-L2`.
