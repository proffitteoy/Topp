# Changelog

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

- Only Windows x64 wheels are published initially; other platforms require an sdist build.
- Generic Wasserstein `(order, internal_p)` combinations are not implemented.
- The C++ interface and ABI are intended for maintenance and experimentation and are not stability commitments.
- Wasserstein optimization is ongoing.
