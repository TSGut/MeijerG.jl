# API

## Main Functions

### Recommended API

```@docs
meijerg
```

### Pure Slater Evaluation

```@docs
meijerg_slater
```

## Calling Conventions

Both `meijerg` and `meijerg_slater` support two calling styles:

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

`meijerg` is the default user-facing function. It applies explicit reductions when available,
uses perturbation for confluent-pole cases, and otherwise delegates to `meijerg_slater`.
