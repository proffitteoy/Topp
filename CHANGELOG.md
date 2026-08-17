# Changelog

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
