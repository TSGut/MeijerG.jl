# Reduction System

`MeijerG.jl` includes a reduction system through the main public function `meijerg(...)`.
This function attempts to recognize special cases and reduce them to elementary or well-known
special functions, then falls back to the pure residue-expansion evaluator `meijerg_slater(...)`
if no match is found.

## Overview

The Meijer G-function is extremely general and encodes many classical special functions
as special cases. When these cases are recognized, direct evaluation via existing functions
(exponential, trigonometric, Bessel, etc.) is typically faster and more numerically accurate
than generic residue summation.

The reduction system currently handles the following categories:

1. **Order cancellation**: Remove redundant parameters
2. **Elementary functions**: Exponential, sine, cosine
3. **Bessel-family reductions**: Bessel J, Bessel K, and the order-0 I continuation
4. **Error/incomplete-gamma reductions**: erf, erfc, lower gamma, upper gamma
5. **Orthogonal polynomials**: Laguerre, Hermite, Jacobi
6. **Logarithmic confluent-pole reduction**: `log1p(z)/z`

## API: `meijerg`

```julia
result = meijerg(a, b, m, n, z)
result = meijerg(a_left, a_right, b_left, b_right, z)
```

If a reduction rule matches, the result is returned directly.
If a confluent-pole case is detected, `meijerg` delegates to the private perturbation routine.
Otherwise, `meijerg` delegates to `meijerg_slater` for pure residue evaluation.

## Category 1: Order Cancellation

When parameters appear in both "left side" and "right side" of the Meijer G-function,
they cancel identically in the mathematical definition.

Specifically:
- If $a_i \in \{a_{\text{right}}\}$ and $a_i \in \{b_{\text{left}}\}$, they cancel
- If $a_i \in \{a_{\text{left}}\}$ and $a_i \in \{b_{\text{right}}\}$, they cancel

The reduced order is computed recursively, so multiple cascades are handled correctly.

### Example

```julia
# G_{2,2}^{1,1}(z | 1.0, 0.25 ; 0.25, 1.0)
# a_left = (1.0), a_right = (0.25)
# b_left = (0.25), b_right = (1.0)
#
# Match: a_left[1] = b_right[1] = 1.0  → cancel
# Match: a_right[1] = b_left[1] = 0.25 → cancel (after first reduction)
# Result: empty, return 1.0

result = meijerg((1.0, 0.25), (0.25, 1.0), 1, 1, 0.7)
# is equivalent to:
result = meijerg((), (), (), (), 0.7)  # both a and b empty → returns 1
```

## Category 2: Elementary Functions

### Exponential: $e^x$

The Meijer G-function

```math
G_{0,1}^{1,0}\!\left(z\;\middle|\;\begin{matrix} - \\ 0 \end{matrix}\right) = e^{-z}
```

**Mapping**: `meijerg((), (), (0,), (), -x)` → `exp(x)`

**Usage**:
```julia
meijerg((), (), (0,), (), -0.5)  # ≈ exp(0.5)
```

### Sine: $\sin(x)$

The Meijer G-function

```math
G_{0,2}^{1,0}\!\left(z\;\middle|\;\begin{matrix} - \\ 1/2, 0 \end{matrix}\right) = \frac{\sin(2\sqrt{z})}{\sqrt{\pi}}
```

**Mapping**: `meijerg((), (), (0.5,), (0,), z^2/4)` with scaling → `sin(z)`

Concretely:
$$\sin(y) = \sqrt{\pi} \cdot G_{0,2}^{1,0}\!\left(\frac{y^2}{4}\;\middle|\;\begin{matrix} - \\ 1/2, 0 \end{matrix}\right)$$

**Usage**:
```julia
y = 0.7
sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2/4)
```

### Cosine: $\cos(x)$

Similarly, the Meijer G-function

```math
G_{0,2}^{1,0}\!\left(z\;\middle|\;\begin{matrix} - \\ 0, 1/2 \end{matrix}\right) = \frac{\cos(2\sqrt{z})}{\sqrt{\pi}}
```

**Mapping**: `meijerg((), (), (0,), (0.5,), z^2/4)` with scaling → `cos(z)`

$$\cos(y) = \sqrt{\pi} \cdot G_{0,2}^{1,0}\!\left(\frac{y^2}{4}\;\middle|\;\begin{matrix} - \\ 0, 1/2 \end{matrix}\right)$$

**Usage**:
```julia
y = 0.7
sqrt(pi) * meijerg((), (), (0,), (0.5,), y^2/4)
```

## Category 3: Bessel Functions

### Bessel J: $J_\nu(x)$

```math
G_{0,2}^{1,0}\!\left(z\;\middle|\;\begin{matrix} - \\ \nu/2, -\nu/2 \end{matrix}\right) = J_\nu(2\sqrt{z})
```

**Mapping**: `meijerg((), (), (ν/2,), (-ν/2,), z)` → `besselj(ν, 2*sqrt(z))`

This rule is applied only for strictly positive real `z` on the real path. Use `z + 0im` to force complex-branch evaluation.

### Bessel I (order-0 continuation on negative reals)

For the special order-zero case, the implementation uses

```math
G_{0,2}^{1,0}\!\left(z\;\middle|\;\begin{matrix} - \\ 0, 0 \end{matrix}\right) = I_0(2\sqrt{-z}), \quad z < 0 \text{ real}
```

to preserve expected real continuation behavior on negative real inputs.

### Modified Bessel K: $K_\nu(x)$

The Meijer G-function

```math
G_{0,2}^{2,0}\!\left(z\;\middle|\;\begin{matrix} - \\ \nu/2, -\nu/2 \end{matrix}\right) = 2K_\nu(2\sqrt{z})
```

The implementation supports the symmetric identity and the shifted-parameter form

```math
G_{0,2}^{2,0}\!\left(z\;\middle|\;\begin{matrix} - \\ b_1, b_2 \end{matrix}\right)
= z^{-c} 2K_\nu(2\sqrt{z}),
\quad c = \frac{b_1+b_2}{2},\; \nu = |b_1-b_2|
```

**Mapping**: `meijerg((), (), (b1, b2), (), z)` → `z^(-c) * 2*besselk(ν, 2*sqrt(z))`

**Usage**:
```julia
ν = 0.8
z = 1.2
lhs = meijerg((), (), (ν/2, -ν/2), (), z)
rhs = 2 * besselk(ν, 2*sqrt(z))
lhs ≈ rhs  # true
```

## Category 4: Error Functions and Incomplete Gamma

### Error Function: $\mathrm{erf}(x)$

```math
G_{1,2}^{1,1}\!\left(z\;\middle|\;\begin{matrix}1\\1/2,0\end{matrix}\right)=\sqrt{\pi}\,\mathrm{erf}(\sqrt{z})
```

**Mapping**: `meijerg((1,), (), (0.5,), (0,), z)` → `sqrt(pi)*erf(sqrt(z))`

### Complementary Error Function: $\mathrm{erfc}(x)$

```math
G_{1,2}^{2,1}\!\left(z\;\middle|\;\begin{matrix}1\\1/2,0\end{matrix}\right)=\sqrt{\pi}\,\mathrm{erfc}(\sqrt{z})
```

**Mapping**: `meijerg((1,), (), (0.5, 0), (), z)` → `sqrt(pi)*erfc(sqrt(z))`

### Lower Incomplete Gamma: $\gamma(a,z)$

```math
G_{1,2}^{1,1}\!\left(z\;\middle|\;\begin{matrix}1\\a,0\end{matrix}\right)=\gamma(a,z)
```

**Mapping**: `meijerg((1,), (), (a,), (0,), z)` → `gamma(a) - gamma(a, z)`

### Upper Incomplete Gamma: $\Gamma(a,z)$

```math
G_{1,2}^{2,1}\!\left(z\;\middle|\;\begin{matrix}1\\a,0\end{matrix}\right)=\Gamma(a,z)
```

**Mapping**: `meijerg((1,), (), (a, 0), (), z)` → `gamma(a, z)`

## Category 5: Orthogonal Polynomials

### Generalized Laguerre: $L_n^{(\alpha)}(z)$

```math
G_{1,2}^{1,0}\!\left(z\;\middle|\;\begin{matrix}n+1\\0,-\alpha\end{matrix}\right)
= \frac{L_n^{(\alpha)}(z)}{\Gamma(n+\alpha+1)},\quad n\in\mathbb{Z}_{\ge 0}
```

**Mapping**: `meijerg((n+1,), (), (0,), (-α,), z)` → `laguerrel(n, α, z) / gamma(n+α+1)`

### Hermite (via Laguerre specializations)

The implementation recognizes the two half-order channels:
- even channel: `b = (0, 1/2)`
- odd channel: `b = (0, -1/2)`

and maps through Laguerre forms consistent with standard Hermite-Laguerre identities.

### Jacobi: $P_n^{(\alpha,\beta)}(x)$

```math
P_n^{(\alpha,\beta)}(x)
=\Gamma(n+\alpha+1)\Gamma(-n-\alpha-\beta)
G_{2,2}^{1,0}\!\left(\frac{x-1}{2}\;\middle|\;\begin{matrix}n+1,-n-\alpha-\beta\\0,-\alpha\end{matrix}\right)
```

With `x = 2z + 1` in the implementation,
`meijerg((n+1, -n-α-β), (0, -α), 1, 0, z)` maps to the corresponding `jacobip` form.

## Category 6: Logarithmic Functions (Confluent Poles)

Confluent-pole Meijer G-functions arise when parameters have integer differences, leading to
logarithmic singularities in the residue expansion. While the full evaluator handles these via
parameter perturbation (4× evaluation + Lagrange extrapolation), direct computation using
standard functions is faster and more accurate.

### Logarithm of (1+z): $\log(1 + z) / z$

The Meijer G-function with confluent poles

```math
G_{2,2}^{1,2}\!\left(z\;\middle|\;\begin{matrix} 1, 1 \\ 1, 0 \end{matrix}\right) = \frac{\log(1+z)}{z}
```

This represents the limiting case where the order-1 parameter difference between numerator
and denominator creates a simple logarithmic pole (DLMF 15.8.2, Gradshteyn & Ryzhik 9.353).

**Mapping**: `meijerg((1, 1), (1, 0), 1, 2, z)` → `log1p(z) / z`

**Usage**:
```julia
z = 0.3
result = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z)
# ≈ log(1.3) / 0.3

# High precision
z_big = BigFloat("0.3")
result_big = meijerg((BigFloat(1), BigFloat(1)), (BigFloat(1), BigFloat(0)), 1, 2, z_big)
# Preserves BigFloat precision throughout
```

## What about the many G → hypergeometric reductions?

The Meijer G-function has many reduction formulas to generalized hypergeometric
functions ${}_pF_q$ or linear combinations thereof. However, the residue expansion used in this package already is this hypergeometric reduction.

To see why, consider the lower expansion for a simple case like
$G_{0,2}^{1,0}(z \mid -;\, \nu/2, -\nu/2)$ (related to Bessel J). With $m=1$, the
residue loop executes once ($k=1$) and produces:

```math
\text{result} = A_1 \cdot z^{b_1} \cdot {}_{0}F_1\!\left(\begin{array}{c}-\\ \beta_1\end{array}\;\middle|\; (-1)^{p-m-n}z\right)
```

That single `pFq(α, β, argument)` call **is** the hypergeometric representation found in many formula collections, derived
analytically from the pole structure.

## References

- **DLMF**: https://dlmf.nist.gov/16.17 (Meijer G-function, special cases)
- **Slater, L.J.** (1966). *Generalized Hypergeometric Functions*. Cambridge University Press.