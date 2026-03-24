# MeijerG.jl

[![CI](https://github.com/TSGut/MeijerG.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/TSGut/MeijerG.jl/actions/workflows/tests.yml)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://tsgut.github.io/MeijerG.jl/)

A Julia package for computing the **Meijer G-function** — a single unifying
function that encompasses almost all classical special functions as particular
cases.  

The package is built on top of
[HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl)
and [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
and supports **arbitrary Julia numeric types**: pass `Float64` for standard
double precision or `BigFloat` (with any precision) for arbitrary-precision
results.

---

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [API reference](#api-reference)
- [Mathematical background](#mathematical-background)
  - [Definition](#definition)
  - [Slater's lower expansion](#slaters-lower-expansion)
  - [Upper expansion](#upper-expansion)
  - [Choosing the expansion](#choosing-the-expansion)
  - [Order reduction](#order-reduction)
  - [Supported special cases](#supported-special-cases)
- [Limitations](#limitations)

---

## Installation

The package is not yet registered; install directly from this repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/TSGut/MeijerG.jl")
```

---

## Quick start

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
```

---

## API reference

### `meijerg(a, b, m, n, z)`

```
meijerg(a, b, m, n, z) -> Number
```

Compute $G_{p,q}^{m,n}(z)$ where `a` and `b` are the complete upper and lower
parameter vectors (or tuples), `m` and `n` are the Meijer G indices
satisfying `0 ≤ m ≤ length(b)` and `0 ≤ n ≤ length(a)`, and `z` is the
evaluation point.

### `meijerg(a_left, a_right, b_left, b_right, z)`

```
meijerg(a_left, a_right, b_left, b_right, z) -> Number
```

Split-parameter form following the Mathematica/mpmath convention.  Here
`n = length(a_left)` and `m = length(b_left)`.  The full parameter tuples
are assembled as `a = [a_left; a_right]` and `b = [b_left; b_right]`.

Both forms accept `Tuple` or `AbstractVector` for the parameter arguments.

**Type promotion.** All parameters and `z` are promoted to a common type via
Julia's standard `promote` mechanism before any computation, so mixing
`Float64` and `BigFloat` inputs works as expected.

---

## Mathematical background

### Definition

The Meijer G-function is defined by the Mellin–Barnes contour integral (DLMF §16.17):

$$G_{p,q}^{m,n}\!\left(z\;\middle|\;\begin{matrix}a_1,\dots,a_p\\b_1,\dots,b_q\end{matrix}\right)
= \frac{1}{2\pi i}\int_{\mathcal{L}}
\frac{\displaystyle\prod_{j=1}^{m}\Gamma(b_j-s)\;\prod_{j=1}^{n}\Gamma(1-a_j+s)}
     {\displaystyle\prod_{j=m+1}^{q}\Gamma(1-b_j+s)\;\prod_{j=n+1}^{p}\Gamma(a_j-s)}
\,z^s\,ds$$

with $0\le m\le q$, $0\le n\le p$, and $z\ne 0$.  The contour $\mathcal{L}$
separates the poles of the $\Gamma(b_j - s)$ factors (lying to the right of
$\mathcal{L}$) from the poles of the $\Gamma(1 - a_j + s)$ factors (lying to
the left), and the integral converges under mild conditions on $p$, $q$ and
$|z|$.

Rather than evaluating this integral numerically, this package uses the
**residue theorem** to replace the contour integral by a finite sum of
generalised hypergeometric functions.  These series representations are exact
and can be evaluated to any precision.

---

### Slater's lower expansion

When $p < q$, or when $p = q$ and $|z| \le 1$, the integral can be closed by
a contour encircling the poles of $\Gamma(b_j - s)$ for $j = 1,\dots,m$.
Summing residues at those poles yields **Slater's theorem** (DLMF 16.17.2):

$$\boxed{G_{p,q}^{m,n}\!\left(z\,\middle|\,\mathbf{a};\,\mathbf{b}\right)
= \sum_{k=1}^{m}
  \underbrace{
    \frac{\displaystyle\prod_{\substack{j=1\\j\ne k}}^{m}\!\!\Gamma(b_j-b_k)
         \;\prod_{j=1}^{n}\Gamma(1+b_k-a_j)}
         {\displaystyle\prod_{j=m+1}^{q}\!\!\Gamma(1+b_k-b_j)
         \;\prod_{j=n+1}^{p}\!\!\Gamma(a_j-b_k)}
  }_{=:\;A_k}
  \;z^{b_k}\;
  {}_{p}F_{q-1}\!\left(\begin{matrix}1+b_k-a_1,\dots,1+b_k-a_p\\
    (1+b_k-b_j)_{j\ne k}\end{matrix}
    \;\middle|\;(-1)^{p-m-n}z\right)}$$

Each term is a power $z^{b_k}$ multiplied by a product of $\Gamma$ values
(the coefficient $A_k$) and a single generalised hypergeometric function
${}_{p}F_{q-1}$ evaluated at the **sign-corrected argument**
$(-1)^{p-m-n}\,z$.

---

### Upper expansion

When $p > q$, or when $p = q$ and $|z| > 1$, the contour is instead closed
around the poles of $\Gamma(1 - a_j + s)$ for $j = 1,\dots,n$, giving an
analogous series in powers of $z^{-1}$:

$$\boxed{G_{p,q}^{m,n}\!\left(z\,\middle|\,\mathbf{a};\,\mathbf{b}\right)
= \sum_{h=1}^{n}
  \underbrace{
    \frac{\displaystyle\prod_{\substack{j=1\\j\ne h}}^{n}\!\!\Gamma(a_h-a_j)
         \;\prod_{j=1}^{m}\Gamma(1-a_h+b_j)}
         {\displaystyle\prod_{j=n+1}^{p}\!\!\Gamma(1-a_h+a_j)
         \;\prod_{j=m+1}^{q}\!\!\Gamma(a_h-b_j)}
  }_{=:\;B_h}
  \;z^{a_h-1}\;
  {}_{q}F_{p-1}\!\left(\begin{matrix}1-a_h+b_1,\dots,1-a_h+b_q\\
    (1-a_h+a_j)_{j\ne h}\end{matrix}
    \;\middle|\;(-1)^{q-m-n}z^{-1}\right)}$$

---

### Choosing the expansion

The package selects the appropriate expansion automatically:

| Condition | Expansion used |
|---|---|
| $p < q$ | lower (always) |
| $p > q$ | upper (always) |
| $p = q$, $\|z\| \le 1$ | lower |
| $p = q$, $\|z\| > 1$ | upper |

---

### Order reduction

Before computing, the package removes cancelling $\Gamma$ factors from
numerator and denominator.  Specifically:

- if $a_k = b_j$ for some $k \le n$ and $j > m$, both parameters are removed
  and $(p,q,n)$ each decrease by one;
- if $a_k = b_j$ for some $k > n$ and $j \le m$, both are removed and
  $(p,q,m)$ each decrease by one.

This mirrors the standard algebraic simplification rule for the G-function
integrand.

---

### Supported special cases

Many well-known functions are particular cases of the Meijer G-function.
Some simple examples, all verified in the test suite:

| Function | G-function identity |
|---|---|
| $e^x$ | $G_{0,1}^{1,0}\!\left(-x\;\middle|\;\begin{smallmatrix}-\\0\end{smallmatrix}\right)$ |
| $\sin x$ | $\sqrt{\pi}\;G_{0,2}^{1,0}\!\left(\tfrac{x^2}{4}\;\middle|\;\begin{smallmatrix}-\\\tfrac{1}{2},0\end{smallmatrix}\right)$ |
| $\cos x$ | $\sqrt{\pi}\;G_{0,2}^{1,0}\!\left(\tfrac{x^2}{4}\;\middle|\;\begin{smallmatrix}-\\0,\tfrac{1}{2}\end{smallmatrix}\right)$ |
| $J_\nu(x)$ | $G_{0,2}^{1,0}\!\left(\tfrac{x^2}{4}\;\middle|\;\begin{smallmatrix}-\\\tfrac{\nu}{2},-\tfrac{\nu}{2}\end{smallmatrix}\right)$ |
| $K_\nu(x)$ | $\tfrac{1}{2}G_{0,2}^{2,0}\!\left(\tfrac{x^2}{4}\;\middle|\;\begin{smallmatrix}-\\\tfrac{\nu}{2},-\tfrac{\nu}{2}\end{smallmatrix}\right)$ |
| $\gamma(\alpha,x)$ | $G_{1,2}^{1,1}\!\left(x\;\middle|\;\begin{smallmatrix}1\\\alpha,0\end{smallmatrix}\right)$ |
| $\Gamma(\alpha,x)$ | $G_{1,2}^{2,0}\!\left(x\;\middle|\;\begin{smallmatrix}1\\\alpha,0\end{smallmatrix}\right)$ |

---

## Limitations

- **Simple poles only.** The current implementation requires that no two of
  $b_1,\dots,b_m$ differ by an integer (lower expansion) and no two of
  $a_1,\dots,a_n$ differ by an integer (upper expansion).  Confluent poles
  introduce logarithmic terms; passing parameters that violate this condition
  raises a `DomainError`.
- **Non-zero argument.** $z = 0$ is not supported and raises a `DomainError`.

