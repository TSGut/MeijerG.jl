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

## Quick usage

```@example
using MeijerG

x = 0.3
v = meijerg((), (), (0,), (), -x)
v
```

```@example
using MeijerG

y = 0.7
sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2/4)
```

See the API and mathematical notes pages for details.
