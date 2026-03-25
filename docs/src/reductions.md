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

The reduction system currently handles three categories:

1. **Order cancellation**: Remove redundant parameters
2. **Elementary functions**: Exponential, sine, cosine
3. **Bessel functions**: Modified Bessel function of the second kind

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
# G_{2,2}^{1,1}(z | 1.0, 0.8 ; 0.25, 1.0)
# a_left = (1.0), a_right = (0.8)
# b_left = (0.25), b_right = (1.0)
# 
# Match: a_left[1] = b_right[1] = 1.0  → cancel
# Match: a_right[1] = b_left[1] = 0.25 → cancel (after first reduction)
# Result: empty, return 1.0

result = meijerg((1.0, 0.8), (0.25, 1.0), 1, 1, z)
# is equivalent to:
result = meijerg((), (), (), (), z)  # both a and b empty → returns 1
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
sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2/4)  # ≈ sin(0.7)
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
sqrt(pi) * meijerg((), (), (0,), (0.5,), y^2/4)  # ≈ cos(0.7)
```

## Category 3: Bessel Functions

### Modified Bessel K: $K_\nu(x)$

The Meijer G-function

```math
G_{0,2}^{2,0}\!\left(z\;\middle|\;\begin{matrix} - \\ \nu/2, -\nu/2 \end{matrix}\right) = 2K_\nu(2\sqrt{z})
```

This is one of the most common Meijer G-function identities, arising in many applied contexts
(probability theory, physics, signal processing).

**Mapping**: `meijerg((), (), (nu/2, -nu/2), (), z)` → `2*besselk(nu, 2*sqrt(z))`

**Usage**:
```julia
ν = 0.8
z = 1.2
lhs = meijerg((), (), (ν/2, -ν/2), (), z)
rhs = 2 * besselk(ν, 2*sqrt(z))
lhs ≈ rhs  # true
```

## Category 4: Logarithmic Functions (Confluent Poles)

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

**Why this matters**: 
- Direct computation avoids integrating the residue series, which has slow convergence near z=0
- Parameter perturbation (the full evaluator's approach) requires 4 independent G-function evalutions + Lagrange extrapolation
- Julia's `log1p` provides superior accuracy for small z, avoiding cancellation loss

**Usage**:
```julia
z = 0.3
result = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z)
# ≈ log(1.3) / 0.3  ≈ 0.8960

# High precision
z_big = BigFloat("0.3")
result_big = meijerg((BigFloat(1), BigFloat(1)), (BigFloat(1), BigFloat(0)), 1, 2, z_big)
# Preserves BigFloat precision throughout
```

## Fallback to Full Evaluation

If an input does not match any reduction rule, `meijerg` silently delegates to `meijerg_slater`:

```julia
# No reduction rule matches this input
a = (0.3, 0.8)
b = (0.2, 1.1)
z = 0.7
result = meijerg(a, b, 1, 1, z)
# Equivalent to:
result = meijerg_slater(a, b, 1, 1, z)
```


## What about the many G → hypergeometric reductions?

The Meijer G-function has a well-known relationship with the generalized hypergeometric
function `{}_pF_q`: many classical special functions can be written as both a Meijer G
and as a specific `{}_pF_q` or linear combinations thereof. However, the residue expansion used in this package already is the hypergeometric reduction.

To see why, consider the lower expansion for a simple case like
$G_{0,2}^{1,0}(z \mid -;\, \nu/2, -\nu/2)$ (related to Bessel J). With $m=1$, the
residue loop executes once ($k=1$) and produces:

```math
	ext{result} = A_1 \cdot z^{b_1} \cdot {}_{0}F_1\!\left(\begin{array}{c}-\\ \beta_1\end{array}\;\middle|\; (-1)^{p-m-n}z\right)
```

That single `pFq(α, β, argument)` call **is** the hypergeometric representation found in many formula collections, derived
analytically from the pole structure.

### Build sin via Meijer G

```julia
using MeijerG
using SpecialFunctions

y = 1.3
# Using reduction
result_reduced = sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2/4)

# Verify
result_direct = sin(y)

@assert result_reduced ≈ result_direct atol=1e-12
```

### Detect order cancellation

```julia
using MeijerG

# This should be recognized and cancelled
a = (2.0, 0.5)     # a = (a_left, a_right) = ((2.0), (0.5))
b = (0.5, 2.0)     # b = (b_left, b_right) = ((0.5), (2.0))
m, n = 1, 1

# a_left[1] = b_right[1] = 2.0 → cancel
# a_right[1] = b_left[1] = 0.5 → cancel
# Result: empty G-function = 1

result = meijerg(a, b, m, n, 0.7)  # Returns 1.0
```

### Mix reductions and fallback

```julia
using MeijerG

# First, reduce order if possible
a = (1.0, 0.3, 0.8)
b = (0.2, 1.0, 1.1)
m, n = 2, 2

# The (1.0, 1.0) pair will cancel
# Remaining: G_{2,2}^{1,1}(z | 0.3, 0.8 ; 0.2, 1.1) → falls back to meijerg
result = meijerg(a, b, m, n, 0.7)
```

## Extending the Reduction System

To add new reduction rules, edit `_reduce_special_case` in `src/reduce.jl`.

Each rule has the form:

```julia
# Rule: Description of mathematical identity
if condition_on_a_b_m_n_z
    z_promoted = float(z)  # Preserve type
    return computed_result
end
```

**Best practices**:
1. Preserve numeric types via `float(z)` or similar promotion
2. Use meaningful comments citing the mathematical reference (DLMF sections, papers)
3. Test with both Float64 and BigFloat inputs
4. Add test cases to `test/test_reductions.jl`

**Example**: Adding a hyperbolic sine reduction

```julia
# G_{0,2}^{1,0}(z | - ; ν/2, -ν/2) with ν=1 → sinh
if a === () && m == 1 && n == 0 && _tuple_isequal(b, (0.5, -0.5))
    z_promoted = float(z)
    return sinh(z_promoted)
end
```

## References

- **DLMF**: https://dlmf.nist.gov/16.17 (Meijer G-function, special cases)
- **Slater, L.J.** (1966). *Generalized Hypergeometric Functions*. Cambridge University Press.
- **Springer & Oldham** (1987). *An Atlas of Functions*. Hemisphere Publishing.
- **Gradshteyn & Ryzhik** (2015). *Table of Integrals, Series, and Products* (8th ed.). Academic Press.
