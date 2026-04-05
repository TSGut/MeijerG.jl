# MeijerG.jl

[![CI](https://github.com/TSGut/MeijerG.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/TSGut/MeijerG.jl/actions/workflows/tests.yml)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://tsgut.github.io/MeijerG.jl/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19430427.svg)](https://doi.org/10.5281/zenodo.19430427)

A Julia package for computing the **Meijer G-function** defined by the Mellin-Barnes integral

$$
G_{p,q}^{m,n}\left(z \hspace{1mm} \left| {a_1,\dots,a_p}\atop{b_1,\dots,b_q}\right) \right.
= \frac{1}{2\pi i}\int_{\mathcal{L}}
\frac{\prod_{j=1}^{m}\Gamma(b_j-s)\prod_{j=1}^{n}\Gamma(1-a_j+s)}
  {\prod_{j=m+1}^{q}\Gamma(1-b_j+s)\prod_{j=n+1}^{p}\Gamma(a_j-s)}
 z^s ds.
$$

The Meijer-G function is a single unifying
function that encompasses many classical special functions as particular
cases and appears, even in its more general forms, in many applications.

The package is built on top of
[HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl)
and [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
and supports **arbitrary Julia numeric types**: pass `Float64`, `Complex{Float64}` for standard
double precision or `BigFloat`, `Complex{BigFloat}` for arbitrary-precision results.

The default algorithm uses **special-case reductions** when available, **Slater expansions** when safe, **adaptive perturbation** for confluent-pole robustness and finally deploys numerical contour integration as a fallback. See the [math notes](https://timon.gutleb.com/MeijerG.jl/dev/math/) in the docs for details.

<table>
<tr>
<td colspan="2" align="center">
Complex phase portraits of $G_{3,3}^{2,1}\!\left(z\,\middle|\,\frac{1}{7}+0.2\cos(2\pi t),\frac{2}{7},\frac{4}{7};\frac{3}{7},\frac{5}{7},\frac{6}{7}\right)$ for $t \in [0,1)$
</td>
</tr>
<tr>
<td align="center">
<img src="https://raw.githubusercontent.com/TSGut/tsgut.github.io/main/images/meijerg_animation_nist.gif" alt="NIST Standard Animation" width="400" />
</td>
<td align="center">
<img src="https://raw.githubusercontent.com/TSGut/tsgut.github.io/main/images/meijerg_animation_viridis.gif" alt="Viridis (Periodic) Animation" width="400" />
</td>
</tr>
<tr>
<td align="center"><strong>NIST Standard</strong></td>
<td align="center"><strong>Viridis (Periodic)</strong></td>
</tr>
</table>

## Contents

- [Installation](#installation)
- [Examples](#examples)
- [Formatter](#formatter)
- [Interface](#interface)
- [Power-user API](#power-user-api)
- [References](#references)

## Installation

The package is not yet registered; install directly from this repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/TSGut/MeijerG.jl")
```

## Examples

```julia
using MeijerG

# exp(x) = G_{0,1}^{1,0}(-x | _ ; 0)
meijerg((), (), (0,), (), -0.3)

# sin(y) = √π · G_{0,2}^{1,0}(y²/4 | _ ; 1/2, 0)
sqrt(pi) * meijerg((), (), (0.5,), (0,), 0.49)

# Arbitrary precision via BigFloat
setprecision(512) do
    x = big"0.3"
    meijerg((), (), (big"0.0",), (), -x)
end

# Split-parameter form: G_{2,2}^{1,1}(z | a_L, a_R ; b_L, b_R)
meijerg((0.25,), (1.75,), (0.5,), (1.25,), 2.0)

# Power-user API: pure Slater residue evaluation, no reductions or perturbation
meijerg_slater((0.25,), (1.75,), (0.5,), (1.25,), 2.0)

# Power-user API: direct Mellin-Barnes contour evaluation
meijerg_contour((0.25,), (1.75,), (0.5,), (1.25,), 2.0)
```

## Formatter

The `@formatter` macro exported by `MeijerG.jl` provides a convenient way to convert Meijer G-function syntax into widely used formats:

```julia
julia> plain, latex, math = @formatter meijerg((1, 2), (3, 4), 1, 1, 5);

==============================================
Plaintext: G_{2,2}^{1,1}(5 | 1, 2 ; 3, 4)

LaTeX: G_{2,2}^{1,1}\left(5\;\middle|\begin{matrix}1, 2\\3, 4\end{matrix}\right)

Mathematica: MeijerG[{{1}, {2}}, {{3}, {4}}, 5]
==============================================
```

## Interface

```julia
meijerg(a, b, m, n, z) -> Number
```

Compute $G_{p,q}^{m,n}(z | {{a}\atop{b}})$  where `a` and `b` are the complete upper and lower
parameter vectors (or tuples), `m` and `n` are the Meijer G indices
satisfying `0 ≤ m ≤ length(b)` and `0 ≤ n ≤ length(a)`, and `z` is the
evaluation point.

```julia
meijerg(a_left, a_right, b_left, b_right, z) -> Number
```

Split-parameter form following the Wolfram Mathematica convention.  Here `n = length(a_left)` and `m = length(b_left)`.  The full parameter tuples are assembled as `a = [a_left; a_right]` and `b = [b_left; b_right]`.

Both forms accept `Tuple` or `AbstractVector` for the parameter arguments.

## Power-user API

```julia
meijerg_slater(a, b, m, n, z) -> Number
meijerg_slater(a_left, a_right, b_left, b_right, z) -> Number
```

Pure Slater residue evaluation without reductions or perturbation. This will fail on confluent-pole cases (e.g., integer-separated parameters). Use only when you 
specifically want raw Slater residue sums.

```julia
meijerg_contour(a, b, m, n, z; kwargs...) -> Number
meijerg_contour(a_left, a_right, b_left, b_right, z; kwargs...) -> Number
```

Direct numerical Mellin-Barnes contour integration.
In production use, `meijerg()` typically handles these cases via adaptive perturbation instead. Use `meijerg_contour` for specific workflows if you know what you are doing.

## References

- [DLMF (Meijer G section)](https://dlmf.nist.gov/16)
- [Wolfram Functions Site (Meijer G)](https://functions.wolfram.com/HypergeometricFunctions/MeijerG/)
- [HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl)
- [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
- [ComplexPhasePortrait.jl](https://github.com/JuliaHolomorphic/ComplexPhasePortrait.jl) (for the generation of complex phase portraits)
