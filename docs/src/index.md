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
- **Two exported evaluation modes**:
  - `meijerg(...)`: Recommended polyalgorithm API with reductions, perturbation handling, and Slater fallback
  - `meijerg_slater(...)`: Pure Slater residue evaluation without reductions or perturbation

## Quick usage

```@example
using MeijerG

x = 0.3
v = meijerg((), (), (0,), (), -x)  # Returns exp(x)
v-exp(x)
```

```@example
using MeijerG

y = 0.7
v = sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2/4)
v-sin(y)
```