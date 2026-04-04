# MeijerG.jl

`MeijerG.jl` computes the Meijer G-function

```math
G_{p,q}^{m,n}\left(z\;\middle|\;\begin{matrix}a_1,\ldots,a_p \\ b_1,\ldots,b_q\end{matrix}\right)
```

using residue-based expansions through `HypergeometricFunctions.jl` and
`SpecialFunctions.jl`.

## Features

- Supports `Float64`, `BigFloat`, and complex versions.
- Supports both full-index and split-parameter calling conventions.
- Includes a formatter macro to convert Julia `meijerg()` calls to plaintext, LaTeX and Wolfram Mathematica compatible forms.
- Uses automatic reductions for common special cases and perturbation for confluent poles.
- Three exported evaluation modes:
  - `meijerg(...)`: Recommended polyalgorithm API with reductions, perturbation handling, and Slater fallback
  - `meijerg_slater(...)`: Pure Slater residue evaluation without reductions or perturbation
  - `meijerg_contour(...)`: Numerical contour-integral evaluation of the Mellin-Barnes definition

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

We support two input modes, one aligning with the mathematical convention and one which aligns with how other software such as Wolfram Mathematica handles Meijer-G function input:

```@example
using MeijerG

# Split form: meijerg((), (), (0,), (), -x)
# Full form:  meijerg(a, b, m, n, z) with a=(), b=(0,), m=1, n=0
x = 0.3
v_split = meijerg((), (), (0,), (), -x)
v_full  = meijerg((), (0,), 1, 0, -x)
v_split - v_full
```

## Real vs. complex inputs

Following Julia conventions like `sqrt` and `HypergeometricFunctions.jl`, this package does not auto-promote real inputs to complex branches.

```@example
using HypergeometricFunctions

try
  HypergeometricFunctions._₂F₁(1/2, 1/2, 1, 2.0)
catch err
  typeof(err)
end
```

```@example
using MeijerG

try
  meijerg((), (1, 0), 2, 0, -1.0)
catch err
  typeof(err)
end
```

For Meijer G, pass complex input explicitly when you want the complex branch:

```@example
using MeijerG

meijerg((), (1, 0), 2, 0, -1.0 + 0im)
```

## Logo

```@example
using MeijerG, ComplexPhasePortrait, ColorSchemes, Images
viridis_colormap = [RGB(c.r, c.g, c.b) for c in ColorSchemes.viridis.colors];
xs = range(-5.0, 5.0; length=1000)
ys = range(-5.0, 5.0; length=1000)
Z  = [x + im*y for y in ys, x in xs]
f(z) = meijerg((1/2,1/2,-1/4,3.0), (0.0,1/2,5,-3/2), 2, 1, z)
img = Images.clamp01nan.(portrait(f.(Z), PTcgrid, colormap=viridis_colormap))
h, w = size(img)
cy, cx = (h + 1) / 2, (w + 1) / 2
r2 = (min(h, w) / 2)^2
disk_logo = [
  (x - cx)^2 + (y - cy)^2 <= r2 ? RGBA(img[y, x], 1) : RGBA(1, 1, 1, 0)
  for y in 1:h, x in 1:w
]
save("assets/logo.png", disk_logo);
```
