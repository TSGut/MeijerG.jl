# MeijerG.jl

`MeijerG.jl` computes the Meijer G-function

```math
G_{p,q}^{m,n}\left(z\;\middle|\;\begin{matrix}a_1,\ldots,a_p \\ b_1,\ldots,b_q\end{matrix}\right)
```

using residue-based expansions through `HypergeometricFunctions.jl` and
`SpecialFunctions.jl`.

## Highlights

- Supports `Float64`, `BigFloat`, and complex arguments.
- Supports both full-index and split-parameter calling conventions.
- Uses automatic reductions for common special cases and perturbation for confluent poles.
- **Two evaluation modes**:
  - `meijerg(...)`: Recommended API with reductions, perturbation handling, and Slater fallback
  - `meijerg_slater(...)`: Pure Slater residue evaluation without reductions or perturbation

## Quick usage

```@example
using MeijerG

x = 0.3
v = meijerg((), (), (0,), (), -x)
v
```

The exponential is recognized directly by `meijerg`:

```@example
using MeijerG

x = 0.3
v = meijerg((), (), (0,), (), -x)  # Returns exp(x)
exp(x)
```

Trigonometric functions via Meijer G:

```@example
using MeijerG

y = 0.7
sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2/4)
sin(y)
```