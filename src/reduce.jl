@doc raw"""
    meijerg_reduce(a, b, m, n, z)
    meijerg_reduce(a_left, a_right, b_left, b_right, z)

Compute the Meijer G-function with explicit reductions to simpler functions for known special cases.

# Overview

This function attempts to recognize and reduce the Meijer G-function to elementary or special functions
(exponential, trigonometric, Bessel K, etc.). If no reduction rule applies, it delegates to the
full `meijerg(...)` evaluator.

The reduction system handles:
- **Order cancellation**: Annihilate matching parameters across left/right partitions
- **Elementary identities**: Exponential, sine, cosine, Bessel-K functions
- **Fallback**: Full Meijer G-function for unrecognized inputs

# Arguments
- `a`, `b`: complete upper and lower parameter collections (tuple/vector form)
- `m`, `n`: Meijer G indices with `0 <= m <= length(b)` and `0 <= n <= length(a)`
- `z`: evaluation point (any Julia number type: Float64, BigFloat, Complex, etc.)

# Alternative signature
The split-parameter form `meijerg_reduce(a_left, a_right, b_left, b_right, z)` is equivalent to
`meijerg_reduce((a_left..., a_right...), (b_left..., b_right...), length(b_left), length(a_left), z)`.
This convention follows the standard notation where `n = |a_left|` and `m = |b_left|`.

# Returns
The computed G-function value, with return type matching the promoted input types.
Supports arbitrary precision via `BigFloat` and complex arguments.

# Examples
```julia
# Exponential: G_{0,1}^{1,0}(z | - ; 0) = exp(-z)
meijerg_reduce((), (), (0,), (), -0.5)  # ≈ exp(0.5)

# Sine: √π · G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2√z)
sqrt(pi) * meijerg_reduce((), (), (0.5,), (0,), 0.1^2/4)

# Bessel K: G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2√z)
meijerg_reduce((), (), (0.4, -0.4), (), 1.0)  # ≈ 2K_{0.8}(2)

# Fallback to full evaluation (no special case)
meijerg_reduce((0.3, 0.8), (0.2, 1.1), 1, 1, 0.7)
```

# See Also
- `meijerg`: Always uses the full residue-expansion evaluator
"""
function meijerg_reduce(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z)
    meijerg_reduce(_totuple(a), _totuple(b), m, n, z)
end

function meijerg_reduce(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    reduced = _reduce_special_case(a, b, m, n, z)
    reduced === nothing && return meijerg(a, b, m, n, z)
    return reduced
end

function meijerg_reduce(a_left::ParameterInput, a_right::ParameterInput,
                        b_left::ParameterInput, b_right::ParameterInput, z)
    a_left_tuple = _totuple(a_left)
    a_right_tuple = _totuple(a_right)
    b_left_tuple = _totuple(b_left)
    b_right_tuple = _totuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg_reduce(a, b, length(b_left_tuple), length(a_left_tuple), z)
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
5. **Bessel K**: `G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2√z)`

Returns `nothing` if no rule matches (caller will fall back to full `meijerg` evaluation).
"""
function _reduce_special_case(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    # Cancellation/order reduction: eliminate exact parameter matches across
    # left/right partitions, then continue reducing recursively.
    reduced_orders = _reduce_orders(a, b, m, n)
    if reduced_orders !== nothing
        a_reduced, b_reduced, m_reduced, n_reduced = reduced_orders
        return meijerg_reduce(a_reduced, b_reduced, m_reduced, n_reduced, z)
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

    # G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2sqrt(z))
    if a === () && m == 2 && n == 0 && length(b) == 2 && isequal(b[2], -b[1])
        z_promoted = float(z)
        ν = float(b[1] - b[2])
        return 2 * besselk(ν, 2 * sqrt(z_promoted))
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
