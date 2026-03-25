# Reduction System

`MeijerG.jl` includes a reduction system accessible via `meijerg_reduce(...)`.
This function attempts to recognize special cases and reduce them to elementary or well-known
special functions, then falls back to the full residue-expansion evaluator if no match is found.

## Overview

The Meijer G-function is extremely general and encodes many classical special functions
as special cases. When these cases are recognized, direct evaluation via existing functions
(exponential, trigonometric, Bessel, etc.) is typically faster and more numerically accurate
than generic residue summation.

The reduction system currently handles three categories:

1. **Order cancellation**: Remove redundant parameters
2. **Elementary functions**: Exponential, sine, cosine
3. **Bessel functions**: Modified Bessel function of the second kind

## API: `meijerg_reduce`

```julia
result = meijerg_reduce(a, b, m, n, z)
result = meijerg_reduce(a_left, a_right, b_left, b_right, z)
```

If a reduction rule matches, the result is returned directly.
Otherwise, `meijerg_reduce` delegates to `meijerg` for the full evaluation.

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

result = meijerg_reduce((1.0, 0.8), (0.25, 1.0), 1, 1, z)
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
meijerg_reduce((), (), (0,), (), -0.5)  # ≈ exp(0.5)
```

### Sine: $\sin(x)$

The Meijer G-function

```math
G_{0,2}^{1,0}\!\left(z\;\middle|\;\begin{matrix} - \\ 1/2, 0 \end{matrix}\right) = \frac{\sin(2\sqrt{z})}{\sqrt{\pi}}
```

**Mapping**: `meijerg_reduce((), (), (0.5,), (0,), z^2/4)` with scaling → `sin(z)`

Concretely:
$$\sin(y) = \sqrt{\pi} \cdot G_{0,2}^{1,0}\!\left(\frac{y^2}{4}\;\middle|\;\begin{matrix} - \\ 1/2, 0 \end{matrix}\right)$$

**Usage**:
```julia
y = 0.7
sqrt(pi) * meijerg_reduce((), (), (0.5,), (0,), y^2/4)  # ≈ sin(0.7)
```

### Cosine: $\cos(x)$

Similarly, the Meijer G-function

```math
G_{0,2}^{1,0}\!\left(z\;\middle|\;\begin{matrix} - \\ 0, 1/2 \end{matrix}\right) = \frac{\cos(2\sqrt{z})}{\sqrt{\pi}}
```

**Mapping**: `meijerg_reduce((), (), (0,), (0.5,), z^2/4)` with scaling → `cos(z)`

$$\cos(y) = \sqrt{\pi} \cdot G_{0,2}^{1,0}\!\left(\frac{y^2}{4}\;\middle|\;\begin{matrix} - \\ 0, 1/2 \end{matrix}\right)$$

**Usage**:
```julia
y = 0.7
sqrt(pi) * meijerg_reduce((), (), (0,), (0.5,), y^2/4)  # ≈ cos(0.7)
```

## Category 3: Bessel Functions

### Modified Bessel K: $K_\nu(x)$

The Meijer G-function

```math
G_{0,2}^{2,0}\!\left(z\;\middle|\;\begin{matrix} - \\ \nu/2, -\nu/2 \end{matrix}\right) = 2K_\nu(2\sqrt{z})
```

This is one of the most common Meijer G-function identities, arising in many applied contexts
(probability theory, physics, signal processing).

**Mapping**: `meijerg_reduce((), (), (nu/2, -nu/2), (), z)` → `2*besselk(nu, 2*sqrt(z))`

**Usage**:
```julia
ν = 0.8
z = 1.2
lhs = meijerg_reduce((), (), (ν/2, -ν/2), (), z)
rhs = 2 * besselk(ν, 2*sqrt(z))
lhs ≈ rhs  # true
```

## Fallback to Full Evaluation

If an input does not match any reduction rule, `meijerg_reduce` silently delegates to `meijerg`:

```julia
# No reduction rule matches this input
a = (0.3, 0.8)
b = (0.2, 1.1)
z = 0.7
result = meijerg_reduce(a, b, 1, 1, z)
# Equivalent to:
result = meijerg(a, b, 1, 1, z)
```

## Comparison: `meijerg` vs `meijerg_reduce`

| Aspect | `meijerg` | `meijerg_reduce` |
|--------|-----------|------------------|
| **Always reduces?** | No | Yes (when possible) |
| **Speed for special cases** | Slower | Faster (direct formula) |

## Examples

### Build sin via Meijer G

```julia
using MeijerG
using SpecialFunctions

y = 1.3
# Using reduction
result_reduced = sqrt(pi) * meijerg_reduce((), (), (0.5,), (0,), y^2/4)

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

result = meijerg_reduce(a, b, m, n, 0.7)  # Returns 1.0
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
result = meijerg_reduce(a, b, m, n, 0.7)
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
