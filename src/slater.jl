const ParameterInput = Union{Tuple, AbstractVector}

"""
    _totuple(xs) -> Tuple

Convert a `Tuple` or `AbstractVector` to a `Tuple`. Identity for tuples.
"""
_totuple(xs::Tuple) = xs
_totuple(xs::AbstractVector) = Tuple(xs)

@doc raw"""
    meijerg_slater(a, b, m, n, z)
    meijerg_slater(a_left, a_right, b_left, b_right, z)

Compute the Meijer G-function using Slater's residue expansions (DLMF §16.17).

This is the **pure Slater evaluator** for advanced users who want direct residue-based computation
without automatic reduction to special functions or perturbation handling for confluent poles.

# Overview

The Meijer G-function is defined via a Mellin–Barnes contour integral whose residues yield
finite sums of generalized hypergeometric functions (pFq). This implementation:
- Uses the lower expansion (Slater, DLMF 16.17.2) for p < q or p = q, |z| ≤ 1
- Uses the upper expansion (Slater, DLMF 16.17.3) for p > q or p = q, |z| > 1
- Supports arbitrary-precision arithmetic via BigFloat and complex numbers
- **Does NOT** apply reductions to special functions
- **Does NOT** use parameter perturbation for confluent poles
# Arguments
- `a`, `b`: complete upper and lower parameter collections (tuple or vector)
  - `p = length(a)` are the total upper parameters
  - `q = length(b)` are the total lower parameters
- `m` (0 ≤ m ≤ q): number of active lower parameters (contributing residues)
- `n` (0 ≤ n ≤ p): number of active upper parameters
- `z`: evaluation point (Float64, BigFloat, Complex variants, nonzero)

# Split-parameter form
`meijerg_slater(a_left, a_right, b_left, b_right, z)` is equivalent to
```
meijerg_slater((a_left..., a_right...), (b_left..., b_right...),
               length(b_left), length(a_left), z)
```

This follows the standard indexing convention where active top parameters come first,
active bottom parameters come first, etc.

# See Also
- `meijerg`: Main public API with special-case reductions and confluent-pole handling
"""
function meijerg_slater(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z)
    meijerg_slater(_totuple(a), _totuple(b), m, n, z)
end

function meijerg_slater(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg_slater is implemented for nonzero z only"))

    return _meijerg_slater_impl(a, b, m, n, z)
end

function meijerg_slater(a_left::ParameterInput, a_right::ParameterInput,
                        b_left::ParameterInput, b_right::ParameterInput, z)
    a_left_tuple = _totuple(a_left)
    a_right_tuple = _totuple(a_right)
    b_left_tuple = _totuple(b_left)
    b_right_tuple = _totuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg_slater(a, b, length(b_left_tuple), length(a_left_tuple), z)
end

"""
    _meijerg_slater_impl(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Inner implementation of Slater's residue expansion, called after input validation.

Selects the expansion mode, checks for pairing violations, and dispatches to
`_lower_expansion` or `_upper_expansion`.
"""
function _meijerg_slater_impl(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    mode = _expansion_mode(p, q, z)
    _validate_pairing(a, b, m, n)
    return _meijerg_simple(a, b, m, n, z, mode)
end

"""
    _meijerg_simple(a::Tuple, b::Tuple, m::Integer, n::Integer, z, mode::Symbol)

Dispatch to the appropriate Slater expansion based on `mode` (`:lower` or `:upper`).
"""
function _meijerg_simple(a::Tuple, b::Tuple, m::Integer, n::Integer, z, mode::Symbol)
    if mode === :lower
        return _lower_expansion(a, b, m, n, z)
    else
        return _upper_expansion(a, b, m, n, z)
    end
end

"""
    _expansion_mode(p::Integer, q::Integer, z) -> Symbol

Determine the Slater expansion direction for a Meijer G-function with `p` upper and `q` lower parameters.

Returns `:lower` (DLMF 16.17.2) when `p < q`, or when `p == q` and `|z| ≤ 1`.
Returns `:upper` (DLMF 16.17.3) when `p > q`, or when `p == q` and `|z| > 1`.
"""
_expansion_mode(p::Integer, q::Integer, z) = p < q ? :lower : p > q ? :upper : abs(z) <= one(abs(z)) ? :lower : :upper

"""
    _promote_inputs(a::Tuple, b::Tuple, z) -> (a_promoted, b_promoted, z_promoted)

Promote all parameters and `z` to a common floating-point type.

Converts `a`, `b`, and `z` to float and broadcasts `promote` across all values
to ensure consistent arithmetic types throughout the expansion.
"""
function _promote_inputs(a::Tuple, b::Tuple, z)
    a_float = map(float, a)
    b_float = map(float, b)
    promoted = promote(float(z), a_float..., b_float...)
    p = length(a)
    q = length(b)
    z_promoted = promoted[1]
    a_promoted = ntuple(i -> promoted[i + 1], p)
    b_promoted = ntuple(i -> promoted[p + i + 1], q)
    return a_promoted, b_promoted, z_promoted
end

"""
    _lower_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Compute the Slater lower expansion (DLMF 16.17.2).

Sums `m` residue terms, each a generalized hypergeometric function `pFq` evaluated
at `±z`. Used when `p < q`, or when `p == q` and `|z| ≤ 1`.
"""
function _lower_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    m == 0 && return zero(z)
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(z_promoted, p - m - n)
    total = zero(z_promoted)
    for k in 1:m
        @inbounds begin
            b_k = b_promoted[k]
            α = ntuple(j -> one(b_k) + b_k - a_promoted[j], p)
            β = Tuple(one(b_k) + b_k - b_promoted[j] for j in 1:q if j != k)
            term = _log_gamma_ratio(
                Iterators.flatten(((b_promoted[j] - b_k for j in 1:m if j != k),
                                   (one(b_k) + b_k - a_promoted[j] for j in 1:n))),
                Iterators.flatten(((one(b_k) + b_k - b_promoted[j] for j in m+1:q),
                                   (a_promoted[j] - b_k for j in n+1:p))),
                b_k)
            term *= _pow_meijerg(z_promoted, b_k)
            term *= pFq(α, β, argument)
            total += term
        end
    end
    return something(total, zero(z_promoted))
end

"""
    _upper_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Compute the Slater upper expansion (DLMF 16.17.3).

Sums `n` residue terms, each a generalized hypergeometric function `pFq` evaluated
at `±1/z`. Used when `p > q`, or when `p == q` and `|z| > 1`.
"""
function _upper_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    n == 0 && return zero(z)
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(inv(z_promoted), q - m - n)
    total = zero(z_promoted)
    for h in 1:n
        @inbounds begin
            a_h = a_promoted[h]
            α = ntuple(j -> one(a_h) - a_h + b_promoted[j], q)
            β = Tuple(one(a_h) - a_h + a_promoted[j] for j in 1:p if j != h)
            term = _log_gamma_ratio(
                Iterators.flatten(((a_h - a_promoted[j] for j in 1:n if j != h),
                                   (one(a_h) - a_h + b_promoted[j] for j in 1:m))),
                Iterators.flatten(((one(a_h) - a_h + a_promoted[j] for j in n+1:p),
                                   (a_h - b_promoted[j] for j in m+1:q))),
                a_h)
            term *= _pow_meijerg(z_promoted, a_h - one(a_h))
            term *= pFq(α, β, argument)
            total += term
        end
    end
    return something(total, zero(z_promoted))
end

"""
    _signed_argument(z, exponent::Integer) -> Number

Return `-z` if `exponent` is odd, otherwise `z`.

Applies the sign convention for the hypergeometric argument in the Slater expansions
(DLMF 16.17.2–3).
"""
_signed_argument(z, exponent::Integer) = isodd(exponent) ? -z : z

"""
    _pow_meijerg(z, exponent) -> Number

Compute `z^exponent`, promoting real negative `z` to complex for non-integer exponents.

For integer-like exponents, rounds to the nearest integer to avoid floating-point drift.
For real negative `z` with a non-integer exponent, wraps in `complex(z)` to return a
well-defined complex power.
"""
function _pow_meijerg(z::Real, exponent)
    if _isintegerlike(exponent)
        return z^Int(round(_realpart(exponent)))
    elseif z >= zero(z)
        return z^exponent
    else
        return complex(z)^exponent
    end
end
_pow_meijerg(z, exponent) = z^exponent

"""
    _log_gamma_ratio(nums, denoms, ref_val) -> Number

Compute `∏Γ(nᵢ) / ∏Γ(dⱼ)` in log-space to avoid intermediate overflow.

For real `ref_val`, accumulates `log|Γ|` and sign using `logabsgamma`, then returns
`sign * exp(log_ratio)`. For complex `ref_val`, accumulates `loggamma` (which encodes
the sign in its imaginary part) and returns `exp(log_ratio)`. Returns `1` when both
iterables are empty.
"""
function _log_gamma_ratio(nums, denoms, ref_val::Real)
    T = typeof(float(ref_val))
    log_ratio = zero(T)
    sign = one(T)
    for v in nums
        (lv, sv) = logabsgamma(v)
        log_ratio += lv
        sign *= sv
    end
    for v in denoms
        (lv, sv) = logabsgamma(v)
        log_ratio -= lv
        sign *= sv
    end
    return sign * exp(log_ratio)
end

function _log_gamma_ratio(nums, denoms, ref_val)
    log_ratio = zero(loggamma(ref_val))
    for v in nums
        log_ratio += loggamma(v)
    end
    for v in denoms
        log_ratio -= loggamma(v)
    end
    return exp(log_ratio)
end

"""
    _validate_pairing(a::Tuple, b::Tuple, m::Integer, n::Integer)

Validate that no active upper–lower parameter pair creates a higher-order pole.

Throws `DomainError` if `a[j] - b[k]` is a positive integer for any `j ≤ n` and `k ≤ m`,
which would produce a higher-order pole in the Mellin–Barnes integrand and invalidate
the simple residue expansion.
"""
function _validate_pairing(a::Tuple, b::Tuple, m::Integer, n::Integer)
    @inbounds for j in 1:n, k in 1:m
        difference = a[j] - b[k]
        _ispositiveintegerlike(difference) && throw(DomainError((a[j], b[k]), "a[j] - b[k] must not be a positive integer for j <= n and k <= m"))
    end
end

"""
    _realpart(x) -> Real

Extract the real part of `x`. Returns `real(x)` for complex values, identity otherwise.
"""
_realpart(x::Complex) = real(x)
_realpart(x) = x

"""
    _isintegerlike(x) -> Bool

Return `true` if `x` represents an integer value, regardless of its numeric type.

Handles `Integer`, `Rational`, `Real`, `Complex`, and unknown types (returns `false`).
"""
_isintegerlike(x::Integer) = true
_isintegerlike(x::Rational) = denominator(x) == 1
_isintegerlike(x::Real) = isfinite(x) && isinteger(x)
_isintegerlike(x::Complex) = iszero(imag(x)) && _isintegerlike(real(x))
_isintegerlike(::Any) = false

"""
    _ispositiveintegerlike(x) -> Bool

Return `true` if `x` is both integer-like (via `_isintegerlike`) and strictly positive.
"""
function _ispositiveintegerlike(x)
    return _isintegerlike(x) && _realpart(x) > zero(_realpart(x))
end

