# MeijerG.jl

`MeijerG.jl` computes the Meijer G-function

```math
G_{p,q}^{m,n}\!\left(z\;\middle|\;\begin{matrix}a_1,\ldots,a_p \\ b_1,\ldots,b_q\end{matrix}\right)
```

using residue-based expansions through `HypergeometricFunctions.jl` and
`SpecialFunctions.jl`.

## Highlights

- Supports `Float64`, `BigFloat`, and complex arguments.
- Supports both full-index and split-parameter calling conventions.
- Includes explicit guardrails for unsupported confluent-pole configurations.
- **Two evaluation modes**:
  - `meijerg(...)`: Always uses true mathematical definition via residue expansions
  - `meijerg_reduce(...)`: Attempts reduction to elementary/special functions, then falls back

## Quick usage

```@example
using MeijerG

x = 0.3
v = meijerg((), (), (0,), (), -x)
v
```

The exponential is also recognized by `meijerg_reduce`:
```@example
using MeijerG

x = 0.3
v = meijerg_reduce((), (), (0,), (), -x)  # Returns exp(x)
exp(x)
```

Trigonometric functions via Meijer G:
```@example
using MeijerG

y = 0.7
sqrt(pi) * meijerg_reduce((), (), (0.5,), (0,), y^2/4)
sin(y)
```

## Where to go next

- **[API Reference](api.md)**: Calling conventions and function signatures
- **[Reduction System](reductions.md)**: Special cases, mathematical identities, and extensibility
- **[Mathematical Notes](math.md)**: Slater expansions and algorithm details
