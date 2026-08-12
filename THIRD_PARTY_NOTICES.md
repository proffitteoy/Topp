# Third-party notices

Topp's runtime kernel is implemented in this repository and has no GUDHI or Hera runtime dependency.

## GUDHI

GUDHI 3.13.0 is used as an optional exact Bottleneck test oracle and as a semantic reference. The historical file `docs/research/patches/0001-gudhi-3.13-small-n-neighbors.patch` describes modifications against GUDHI source and is retained only for research provenance.

GUDHI is licensed under the MIT License:

> Copyright (c) 2014 The GUDHI developers.

The GUDHI copyright and permission notice apply to substantial portions of GUDHI represented by that patch. See <https://github.com/GUDHI/gudhi-devel/blob/master/LICENSE>.

## NumPy, pybind11, and scikit-build-core

NumPy is Topp's runtime dependency. pybind11 and scikit-build-core are build dependencies. They are distributed separately under their respective licenses and are not vendored in this repository.
