# Installation

[中文](https://proffitteoy.github.io/Topp/getting-started/installation.html)

## Windows and Linux wheels

Topp publishes wheels for Windows x64, Linux x86_64, and CPython 3.10–3.14:

```console
py -m pip install topp
```

Confirm the installed version:

```pycon
>>> import topp
>>> topp.__version__
'1.0.1'
```

NumPy 1.23 or newer is installed automatically and is Topp's only runtime dependency.

## Virtual environment

Using an isolated environment is recommended:

```console
py -m venv .venv
.venv\Scripts\activate
py -m pip install topp
```

## Source build

A source build requires CMake 3.24 or newer and a C++20 compiler:

```console
git clone https://github.com/proffitteoy/Topp.git
cd Topp
py -m pip install -v .
```

Linux is covered by release CI; macOS is not. A successful local macOS compilation is not currently a published support guarantee.

## Smoke test

```python
import topp

assert topp.bottleneck_distance([[0.0, 1.0]], []) == 0.5
print("Topp is ready")
```

Expected output:

```text
Topp is ready
```
