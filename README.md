# MeijerG.jl

[![CI](https://github.com/TSGut/MeijerG.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/TSGut/MeijerG.jl/actions/workflows/tests.yml)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://tsgut.github.io/MeijerG.jl/)

A Julia package for computing the **Meijer G-function** defined by the Mellin-Barnes integral
$$
G_{p,q}^{m,n}\left(z \hspace{1mm} \left| {a_1,\dots,a_p}\atop{b_1,\dots,b_q}\right) \right.
= \frac{1}{2\pi i}\int_{\mathcal{L}}
\frac{\prod_{j=1}^{m}\Gamma(b_j-s)\prod_{j=1}^{n}\Gamma(1-a_j+s)}
  {\prod_{j=m+1}^{q}\Gamma(1-b_j+s)\prod_{j=n+1}^{p}\Gamma(a_j-s)}
 z^s ds.
$$
The Meijer-G function is a single unifying
function that encompasses almost many classical special functions as particular
cases and appears, even in its more general forms, in many applications.

The package is built on top of
[HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl)
and [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
and supports **arbitrary Julia numeric types**: pass `Float64`, `Complex{Float64}` for standard
double precision or `BigFloat` (with any precision) for arbitrary-precision
results.

For an explanation of how this package computes the Meijer-G function see the [Math notes](docs/src/math.md). in the docs.

---

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [API reference](#api-reference)
- [Limitations](#limitations)

---

## Installation

The package is not yet registered; install directly from this repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/TSGut/MeijerG.jl")
```

---

## Examples

```julia
using MeijerG

# exp(x) = G_{0,1}^{1,0}(-x | _ ; 0)
meijerg((), (), (0,), (), -0.3)        # ≈ exp(0.3) ≈ 1.34985...

# sin(y) = √π · G_{0,2}^{1,0}(y²/4 | _ ; 1/2, 0)
sqrt(pi) * meijerg((), (), (0.5,), (0,), 0.49)  # ≈ sin(0.7) ≈ 0.64421...

# Arbitrary precision via BigFloat
setprecision(512) do
    x = big"0.3"
    meijerg((), (), (big"0.0",), (), -x)   # exp(x) to ~150 decimal places
end

# Split-parameter form: G_{2,2}^{1,1}(z | a_L, a_R ; b_L, b_R)
meijerg((0.25,), (1.75,), (0.5,), (1.25,), 2.0)

# Power-user API: pure Slater residue evaluation, no reductions or perturbation
meijerg_slater((0.25,), (1.75,), (0.5,), (1.25,), 2.0)
```

---

## Supported methods

### `meijerg(a, b, m, n, z)`

```
meijerg(a, b, m, n, z) -> Number
```

Compute $G_{p,q}^{m,n}(z)$ where `a` and `b` are the complete upper and lower
parameter vectors (or tuples), `m` and `n` are the Meijer G indices
satisfying `0 ≤ m ≤ length(b)` and `0 ≤ n ≤ length(a)`, and `z` is the
evaluation point.

This is the recommended user-facing API. It applies explicit reductions to simpler functions when
available, uses perturbation for confluent-pole cases, and otherwise delegates to `meijerg_slater`.

### `meijerg(a_left, a_right, b_left, b_right, z)`

```
meijerg(a_left, a_right, b_left, b_right, z) -> Number
```

Split-parameter form following the Wolfram Mathematica and mpmath convention.  Here
`n = length(a_left)` and `m = length(b_left)`.  The full parameter tuples
are assembled as `a = [a_left; a_right]` and `b = [b_left; b_right]`.

Both forms accept `Tuple` or `AbstractVector` for the parameter arguments.

### `meijerg_slater(a, b, m, n, z)`

```
meijerg_slater(a, b, m, n, z) -> Number
```

Pure Slater residue evaluation with no explicit reductions and no perturbation handling.
This is intended for advanced users who want direct access to the residue-expansion algorithm.