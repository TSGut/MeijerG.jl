# API

## Main Functions

### Recommended API

```@docs
meijerg
```

### Formatter Macro

For producing human- and tool-readable representations of Meijer G-functions (plaintext, LaTeX, and Mathematica), use the `@formatter` macro. The macro takes a `meijerg(...)` call and returns a 3-tuple `(plaintext, latex, mathematica)`.

Example:

```julia
plain, latex, math = @formatter meijerg((1,2), (3,4), 1, 1, 5)
```

You can also optionally pass a `verbose` value to disable the automatic printing.

Example:

```julia
plain, latex, math = @formatter meijerg((1,2), (3,4), 1, 1, 5) verbose=false
```

```@docs
@formatter
meijerg_to_latex_and_math
```

### Pure Slater Evaluation

```@docs
meijerg_slater
```

### Log-Ratio Jet Evaluation

```@docs
meijerg_jet
```

### Exact-Confluence Laurent Evaluation

```@docs
meijerg_laurent
```

### Contour Evaluation

```@docs
meijerg_contour
```

## Calling Conventions

`meijerg`, `meijerg_slater`, `meijerg_jet`, `meijerg_laurent`, and `meijerg_contour` support
two calling styles:

### Full-parameter form

```julia
result = meijerg(a, b, m, n, z)
```

where `a` and `b` are tuples or vectors of parameters.

### Split-parameter form

```julia
result = meijerg(a_left, a_right, b_left, b_right, z)
```

Equivalent to `meijerg((a_left..., a_right...), (b_left..., b_right...),
                        length(b_left), length(a_left), z)`

## Working precision

Every residue sum `meijerg` evaluates, the confluent jet series as well as the separated Slater
expansion, carries a summation condition number `κ_sum = Σ|T_n| / |Σ T_n|`. Roughly
`log10(κ_sum)` digits are lost to cancellation, which happens off the convergent
regime of a branch (`|z| > 1` for the lower, `|z| < 1` for the upper expansion).
A sum that is off its convergent regime or has `κ_sum > 1e2` is
evaluated internally in doubled working precision: `Float32` inputs on the `Float64`
path, `Float64` inputs in an internal double-double kernel (about five times the cost
of the `Float64` kernel for the exact-pair series, more for the Slater expansion
and the finite-separation and higher-order jets), and `BigFloat` inputs with more
bits. Contour integration remains the fallback beyond `κ_sum ≈ 1e14`. The returned
type is always the float-promoted input type.

```julia
meijerg((0.3,), (0.2, 1.2), 2, 1, 10.0)                            # double-double internally
```

`meijerg_jet` exposes this choice for the jet series through its `working_precision`
keyword: `:auto` as in `meijerg`, `:input` for the uncertified input-precision sum, and
`:extended` to always sum in doubled precision.
