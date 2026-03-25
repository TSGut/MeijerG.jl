# API

## Main Functions

### Pure Evaluation

```@docs
meijerg
```

### Reduction with Special Cases

```@docs
meijerg_reduce
```

## Calling Conventions

Both `meijerg` and `meijerg_reduce` support two calling styles:

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
