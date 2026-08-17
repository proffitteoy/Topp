# Topp v1.0.0

Topp 1.0 freezes the public Python API and the validated adaptive exact kernels as the stable release baseline.

## Highlights

- Exact Bottleneck distance with prepared, one-to-many, caller-output, and threshold-decision APIs.
- Exact `W1-L∞` and `W2-L2` Wasserstein distance with prepared diagrams, native batch calls, and reusable workspace.
- Adaptive Bottleneck routing for duplicate-heavy, separated, near-diagonal, sparse, and general inputs.
- Adaptive Wasserstein candidate, graph, matching, and component paths selected from the completed kernel experiments.
- Chinese-default bilingual documentation with explicit input, mathematical, platform, and benchmark boundaries.

## Compatibility

No migration is required from `0.1.0`. The documented Python API now follows the `1.x` compatibility contract. The C++ headers and ABI remain maintenance interfaces and are not covered by that promise.

## Validation target

- Windows x64 wheels for CPython 3.10–3.14 and an sdist.
- C++ exact Bottleneck and Wasserstein regressions with AVX2 and scalar builds.
- Python API, brute-force, GUDHI, and POT oracle checks.
- Strict Chinese and English Sphinx builds and fresh package smoke tests.

## Known limits

- Linux and macOS are not release-CI platforms.
- Wasserstein supports only exact `W1-L∞` and `W2-L2`.
- Public benchmark summaries are bounded to the documented machine, inputs, and call patterns.
