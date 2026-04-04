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

### Contour Evaluation

```@docs
meijerg_contour
```

## Calling Conventions

`meijerg`, `meijerg_slater`, and `meijerg_contour` support two calling styles:

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

