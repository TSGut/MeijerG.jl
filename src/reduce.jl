@doc raw"""
    meijerg(a, b, m, n, z)
    meijerg(a_left, a_right, b_left, b_right, z)

Compute the Meijer G-function with automatic reduction to simpler functions and confluent-pole handling.

# Overview

This is the **primary public API** and recommended entry point. It combines three strategies:
1. **Special-case reduction**: Recognizes and reduces to elementary or special functions 
    (exponential, trigonometric, Bessel J/K, incomplete gamma, logarithm, etc.) for known patterns
2. **Confluent-pole handling**: Detects integer-separated parameters that create logarithmic singularities
   and uses parameter perturbation + Lagrange extrapolation for accurate evaluation
3. **Slater expansion fallback**: For unrecognized inputs, uses pure Slater residue expansions

# Arguments
- `a`, `b`: complete upper and lower parameter collections (tuple/vector form)
- `m`, `n`: Meijer G indices with `0 <= m <= length(b)` and `0 <= n <= length(a)`
- `z`: evaluation point (any Julia number type: Float64, BigFloat, Complex, etc.)

# Alternative signature
The split-parameter form `meijerg(a_left, a_right, b_left, b_right, z)` is equivalent to
`meijerg((a_left..., a_right...), (b_left..., b_right...), length(b_left), length(a_left), z)`.
This convention follows the standard notation where `n = |a_left|` and `m = |b_left|`.

# Returns
The computed G-function value, with return type matching the promoted input types.
Supports arbitrary precision via `BigFloat` and complex arguments.

# Examples
```julia
# Exponential: G_{0,1}^{1,0}(z | - ; 0) = exp(-z)
meijerg((), (), (0,), (), -0.5)  # ≈ exp(0.5)

# Sine: √π · G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2√z)
sqrt(pi) * meijerg((), (), (0.5,), (0,), 0.1^2/4)

# Bessel K: G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2√z)
meijerg((), (), (0.4, -0.4), (), 1.0)  # ≈ 2K_{0.8}(2)

# Bessel J: G_{0,2}^{1,0}(z | - ; ν/2, -ν/2) = J_ν(2√z)
meijerg((), (), (0.4,), (-0.4,), 1.0)  # ≈ J_{0.8}(2)

# Lower incomplete gamma: G_{1,2}^{1,1}(z | 1 ; a, 0) = γ(a,z)
meijerg((1,), (), (0.7,), (0,), 1.3)  # ≈ γ(0.7, 1.3)

# Upper incomplete gamma: G_{1,2}^{2,1}(z | 1 ; a, 0) = Γ(a,z)
meijerg((1,), (), (0.7, 0), (), 1.3)  # ≈ Γ(0.7, 1.3)

# Logarithm (1+z): G_{2,2}^{1,2}(z | 1, 1 ; 1, 0) = log(1+z)/z
meijerg((1, 1), (1, 0), 1, 2, 0.5)  # ≈ log(1.5)/0.5

# Confluent-pole case: automatic perturbation handling
meijerg((), (), (1.5, 0.5), (), z)  # Integer difference detected, uses perturbation

# Generic input: falls back to Slater expansion
meijerg((0.3, 0.8), (0.2, 1.1), 1, 1, 0.7)
```

# See Also
- `meijerg_slater`: Pure Slater residue expansion (power users, no reductions and unsafe for confluent poles)
"""
function meijerg(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z)
    meijerg(_totuple(a), _totuple(b), m, n, z)
end

function meijerg(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg is implemented for nonzero z only"))
    
    # Step 1: Check for special-case reductions
    reduced = _reduce_special_case(a, b, m, n, z)
    if reduced !== nothing
        return reduced
    end
    
    # Step 2: Check for confluent poles and use perturbation if needed
    mode = _expansion_mode(p, q, z)
    if _has_confluent_poles(a, b, m, n, mode)
        return _meijerg_perturb(a, b, m, n, z)
    end
    
    # Step 3: Fall back to pure Slater expansion
    return meijerg_slater(a, b, m, n, z)
end

function meijerg(a_left::ParameterInput, a_right::ParameterInput,
                 b_left::ParameterInput, b_right::ParameterInput, z)
    a_left_tuple = _totuple(a_left)
    a_right_tuple = _totuple(a_right)
    b_left_tuple = _totuple(b_left)
    b_right_tuple = _totuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg(a, b, length(b_left_tuple), length(a_left_tuple), z)
end

"""
    _has_confluent_poles(a::Tuple, b::Tuple, m::Integer, n::Integer, mode::Symbol)

Detect if the Meijer G-function input exhibits confluent poles.

Confluent poles arise when active parameters have integer differences, creating logarithmic
singularities in the residue expansion. Detection depends on expansion mode:
- **Lower expansion** (p < q): Check for integer differences in active lower parameters `b[1:m]`
- **Upper expansion** (p > q): Check for integer differences in active upper parameters `a[1:n]`

# Returns
`true` if confluent poles detected, `false` otherwise.
"""
function _has_confluent_poles(a::Tuple, b::Tuple, m::Integer, n::Integer, mode::Symbol)
    if mode === :lower
        @inbounds for j in 1:m, k in j+1:m
            _isintegerlike(b[j] - b[k]) && return true
        end
    else
        @inbounds for j in 1:n, k in j+1:n
            _isintegerlike(a[j] - a[k]) && return true
        end
    end
    return false
end

"""
    _reduce_special_case(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Attempt to recognize special cases and reduce to simpler functions.

This is the core dispatcher for the reduction system. It checks for reduction rules in order:
1. **Order cancellation**: Parameters appearing in both `a_left` and `b_right` (or vice versa)
   cancel out, reducing problem size recursively.
2. **Exponential**: `G_{0,1}^{1,0}(z | - ; 0) = exp(-z)`
3. **Sine**: `G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2√z) / √π`
4. **Cosine**: `G_{0,2}^{1,0}(z | - ; 0, 1/2) = cos(2√z) / √π`
5. **Bessel J**: `G_{0,2}^{1,0}(z | - ; ν/2, -ν/2) = J_ν(2√z)`
6. **Bessel K**: `G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2√z)`
7. **Lower incomplete gamma**: `G_{1,2}^{1,1}(z | 1 ; a, 0) = γ(a,z)`
8. **Upper incomplete gamma**: `G_{1,2}^{2,1}(z | 1 ; a, 0) = Γ(a,z)`
9. **Logarithm (1+z)**: `G_{2,2}^{1,2}(z | 1, 1 ; 1, 0) = log(1+z) / z` (direct evaluation via `log1p`)

Returns `nothing` if no rule matches (caller will fall back to full `meijerg` evaluation).
"""
function _reduce_special_case(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    # Cancellation/order reduction: eliminate exact parameter matches across
    # left/right partitions, then continue reducing recursively.
    reduced_orders = _reduce_orders(a, b, m, n)
    if reduced_orders !== nothing
        a_reduced, b_reduced, m_reduced, n_reduced = reduced_orders
        return meijerg(a_reduced, b_reduced, m_reduced, n_reduced, z)
    end

    # G_{0,1}^{1,0}(z | - ; 0) = exp(-z)
    if a === () && m == 1 && n == 0 && _tuple_isequal(b, (0,))
        z_promoted = float(z)
        return exp(-z_promoted)
    end

    # G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2sqrt(z))/sqrt(pi)
    if a === () && m == 1 && n == 0 && _tuple_isequal(b, (0.5, 0))
        z_promoted = float(z)
        return sin(2 * sqrt(z_promoted)) * _inv_sqrt_pi(z_promoted)
    end

    # G_{0,2}^{1,0}(z | - ; 0, 1/2) = cos(2sqrt(z))/sqrt(pi)
    if a === () && m == 1 && n == 0 && _tuple_isequal(b, (0, 0.5))
        z_promoted = float(z)
        return cos(2 * sqrt(z_promoted)) * _inv_sqrt_pi(z_promoted)
    end

    # G_{0,2}^{1,0}(z | - ; ν/2, -ν/2) = J_ν(2sqrt(z))
    if a === () && m == 1 && n == 0 && length(b) == 2 && isequal(b[2], -b[1])
        z_promoted = float(z)
        ν = float(b[1] - b[2])
        return besselj(ν, 2 * sqrt(z_promoted))
    end

    # G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2sqrt(z))
    if a === () && m == 2 && n == 0 && length(b) == 2 && isequal(b[2], -b[1])
        z_promoted = float(z)
        ν = float(b[1] - b[2])
        return 2 * besselk(ν, 2 * sqrt(z_promoted))
    end

    # G_{1,2}^{1,1}(z | 1 ; a, 0) = γ(a,z)
    if n == 1 && m == 1 && _tuple_isequal(a, (1,)) && length(b) == 2 && isequal(b[2], zero(b[2]))
        z_promoted = float(z)
        α = float(b[1])
        return gamma(α) - gamma(α, z_promoted)
    end

    # G_{1,2}^{2,1}(z | 1 ; a, 0) = Γ(a,z)
    if n == 1 && m == 2 && _tuple_isequal(a, (1,)) && length(b) == 2 && isequal(b[2], zero(b[2]))
        z_promoted = float(z)
        α = float(b[1])
        return gamma(α, z_promoted)
    end

    # Logarithmic confluent poles: G_{2,2}^{1,2}(z | 1, 1 ; 1, 0)
    # Maps to log(1+z) / z via parameter perturbation limit (Gradshteyn & Ryzhik 9.353, DLMF 15.8.2)
    # Direct computation avoids expensive 4× evaluation + Lagrange extrapolation
    if length(a) == 2 && length(b) == 2 && m == 1 && n == 2 && 
       _tuple_isequal(a, (1, 1)) && _tuple_isequal(b, (1, 0))
        z_promoted = float(z)
        # Avoid division by zero; use log1p for accuracy
        if iszero(z_promoted)
            return z_promoted
        else
            return log1p(z_promoted) / z_promoted
        end
    end

    return nothing
end

"""
    _reduce_orders(a::Tuple, b::Tuple, m::Integer, n::Integer)

Eliminate parameter cancellations to reduce Meijer G-function order.

When a parameter appears in both `a_left = a[1:n]` and `b_right = b[m+1:end]`,
the two cancel out identically (both contribute a factor to the Mellin–Barnes integral
that annihilate). Similarly for parameters in `a_right = a[n+1:end]` and `b_left = b[1:m]`.

This function repeatedly removes such pairs until no more cancellations remain.

# Returns
- `nothing` if no reductions occur (i.e., input is already reduced)
- `(a_reduced, b_reduced, m_reduced, n_reduced)` tuple with cancelled parameters removed

# Example
```
a = (1.0, 0.25)
b = (0.25, 1.0)
m = 1, n = 1

# a_left = (1.0), a_right = (0.25), b_left = (0.25), b_right = (1.0)
# Match: a_left[1] == b_right[1] (both 1.0)
# Match: a_right[1] == b_left[1] (both 0.25)
# Result: a_reduced=(), b_reduced=(), m_reduced=0, n_reduced=0
```
"""
function _reduce_orders(a::Tuple, b::Tuple, m::Integer, n::Integer)
    a_left = collect(a[1:n])
    a_right = collect(a[n+1:end])
    b_left = collect(b[1:m])
    b_right = collect(b[m+1:end])

    changed = true
    while changed
        changed = false

        for i in eachindex(a_left)
            j = findfirst(x -> isequal(x, a_left[i]), b_right)
            if j !== nothing
                deleteat!(a_left, i)
                deleteat!(b_right, j)
                n -= 1
                changed = true
                break
            end
        end
        changed && continue

        for i in eachindex(a_right)
            j = findfirst(x -> isequal(x, a_right[i]), b_left)
            if j !== nothing
                deleteat!(a_right, i)
                deleteat!(b_left, j)
                m -= 1
                changed = true
                break
            end
        end
    end

    reduced = (Tuple(vcat(a_left, a_right)), Tuple(vcat(b_left, b_right)), m, n)
    return reduced == (a, b, m, n) ? nothing : reduced
end

"""
    _tuple_isequal(left::Tuple, right::Tuple) -> Bool

Check tuple equality by element-wise `isequal`.

Used for robust pattern matching of parameter tuples, which may contain
floating-point values with NaN or unusual representations. Regular `==` may
behave differently for floating-point comparisons.

# Example
```
_tuple_isequal((0.5, 0), (0.5, 0))  # true
_tuple_isequal((0.5,), (0,))        # false (different length and values)
```
"""
function _tuple_isequal(left::Tuple, right::Tuple)
    length(left) == length(right) || return false
    for i in eachindex(left)
        isequal(left[i], right[i]) || return false
    end
    return true
end

"""
    _inv_sqrt_pi(z) -> Number

Compute `1 / √π` preserving the type of `z`.

For type stability, promotes `π` to the type of `z` before division.
Ensures BigFloat inputs produce BigFloat outputs with correct precision.

# Example
```
_inv_sqrt_pi(Float64(1))      # Float64: ≈ 0.5641895835477563
_inv_sqrt_pi(BigFloat(1))     # BigFloat with current precision
```
"""
_inv_sqrt_pi(z) = inv(sqrt(convert(typeof(z), pi)))
