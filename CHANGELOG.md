# Changelog

## 1.0.1 - 2026-08-21

Python API validation and error-contract patch release.

### Changed

- Diagram inputs now require real numeric values. Complex arrays, masked arrays,
  booleans, numeric strings, and other non-real values are rejected instead of
  being silently coerced.
- `DiagramLike` is now a documented top-level runtime and typing export, and
  `prepare_diagram` explicitly accepts an existing `PreparedDiagram`.
- Linux x86_64 joins Windows x64 as a wheel, C++/Python CI, sdist, and
  TestPyPI-verified release platform for CPython 3.10–3.14.

### Fixed

- Finite values that cannot be represented as `float64` now raise `ValueError`
  instead of overflowing to essential points or leaking `OverflowError`.
- Batch validation checks metric, query, and output layout before consuming
  targets, and target failures report their `diagrams[index]` location without
  hiding the original exception category.
- Caller-provided output arrays must be aligned in both the Python wrapper and
  native binding.

## 1.0.0 - 2026-08-17

Stable Python API release with the optimized exact kernels frozen as the default implementation.

### Added

- A Chinese-default, independently accessible bilingual documentation site with installation, task guides, API, mathematics, platform support, and bounded benchmark results.
- Exact duplicate-capacity, no-cross, mandatory-flow, geometric-refinement, and quickselect routes behind the Bottleneck adaptive dispatcher.
- Adaptive Wasserstein candidate generation, graph representation, matching, component decomposition, and native batch/workspace reuse.

### Changed

- The documented Python API now follows the `1.x` compatibility contract; no migration is required from `0.1.0`.
- Package metadata now declares production/stable status, while the C++ headers and ABI remain maintenance interfaces without a stability guarantee.
- Release automation validates that the Git tag matches package metadata and creates a stable GitHub release.
- Research material is frozen around the selected default paths; superseded proposals and one-off negative-experiment harnesses are removed from the active tree.

### Known limitations

- Published wheels cover Windows x64 and CPython 3.10–3.14; Linux and macOS source builds are not release-CI platforms.
- Wasserstein supports exact `W1-L∞` and `W2-L2`; other `(order, internal_p)` pairs raise `NotImplementedError`.
- The C++ interface and ABI are not part of the stable `1.x` Python compatibility contract.

## 0.1.0 - 2026-08-12

Public preview release.

### Added

- Exact Bottleneck distance with prepared, one-to-many, caller-output, and threshold-decision APIs.
- Exact `W1-L∞` and `W2-L2` Wasserstein distances with prepared and native one-to-many APIs.
- Strict Python input validation, immutable `PreparedDiagram`, GIL-releasing native calls, type information, and Windows x64 wheels for CPython 3.10–3.14.
- Bilingual user, API, and development documentation.

### Experimental

- C++ candidate, graph, matching, component, duplicate-compression, and incremental-pricing strategies remain available for kernel research but are not exposed in the public Python API.

### Known limitations

- Only Windows x64 wheels are published initially; Linux and macOS source builds are not yet covered by CI.
- Generic Wasserstein `(order, internal_p)` combinations are not implemented.
- The C++ interface and ABI are intended for maintenance and experimentation and are not stability commitments.
- Wasserstein optimization is ongoing.
