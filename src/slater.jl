const ParameterInput = Union{Tuple, AbstractVector}
const _SLATER_BOUNDARY_BIGFLOAT_RADIUS = 0.02

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
- Promotes internally to `BigFloat` near `|z| = 1` for balanced cases to improve accuracy
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
    meijerg_slater(Tuple(a), Tuple(b), m, n, z)
end

function meijerg_slater(a_left::ParameterInput, a_right::ParameterInput,
                        b_left::ParameterInput, b_right::ParameterInput, z)
    a_left_tuple = Tuple(a_left)
    a_right_tuple = Tuple(a_right)
    b_left_tuple = Tuple(b_left)
    b_right_tuple = Tuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg_slater(a, b, length(b_left_tuple), length(a_left_tuple), z)
end

function meijerg_slater(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    # Validate input parameters
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg_slater is implemented for nonzero z only"))
    _validate_pairing(a, b, m, n)

    # Use internal BigFloat arithmetic near the balanced branch boundary abs(z) = 1.
    if p == q && abs(abs(z) - one(abs(z))) <= _SLATER_BOUNDARY_BIGFLOAT_RADIUS
        return _meijerg_slater_boundary_bigfloat(a, b, m, n, z)
    end

    # Determine expansion mode and dispatch
    mode = _expansion_mode(p, q, z)
    if mode === :lower
        return _lower_expansion(a, b, m, n, z)
    else # mode === :upper
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

Helper function to promote all parameters and `z` to a common numeric type and return appropriate tuples.
"""
function _promote_inputs(a::Tuple, b::Tuple, z)
    promoted = promote(z, a..., b...)
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
    T = eltype(b_promoted)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(z_promoted, p - m - n)
    total = zero(z_promoted)
    for k in 1:m
        @inbounds begin
            b_k = b_promoted[k]
            α = ntuple(j -> one(T) + b_k - a_promoted[j], p)
            β = Tuple(one(T) + b_k - b_promoted[j] for j in 1:q if j != k)
            term = _gamma_prod_ratio(
                Iterators.flatten(((b_promoted[j] - b_k for j in 1:m if j != k),
                                   (one(T) + b_k - a_promoted[j] for j in 1:n))),
                Iterators.flatten(((one(T) + b_k - b_promoted[j] for j in m+1:q),
                                   (a_promoted[j] - b_k for j in n+1:p))), 
                                   T)

            # Skip evaluating the hypergeometric if the gamma-ratio is zero to
            # avoid 0 * NaN when `pFq` is undefined due to nonpositive integer
            # denominator parameters.
            if iszero(term)
                continue
            end

            powv = z_promoted ^ b_k
            if iszero(powv)
                continue
            end

            term *= powv * pFq(α, β, argument)
            total += term
        end
    end
    return total
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
    T = eltype(a_promoted)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(inv(z_promoted), q - m - n)
    total = zero(z_promoted)
    for h in 1:n
        @inbounds begin
            a_h = a_promoted[h]
            α = ntuple(j -> one(T) - a_h + b_promoted[j], q)
            β = Tuple(one(T) - a_h + a_promoted[j] for j in 1:p if j != h)
            term = _gamma_prod_ratio(
                Iterators.flatten(((a_h - a_promoted[j] for j in 1:n if j != h),
                                   (one(T) - a_h + b_promoted[j] for j in 1:m))),
                Iterators.flatten(((one(T) - a_h + a_promoted[j] for j in n+1:p),
                                   (a_h - b_promoted[j] for j in m+1:q))), 
                                   T)

            # Skip evaluating the hypergeometric if the gamma-ratio is zero to
            # avoid 0 * NaN when `pFq` is undefined due to nonpositive integer
            # denominator parameters.
            if iszero(term)
                continue
            end

            powv = z_promoted^(a_h - one(T))
            if iszero(powv)
                continue
            end

            term *= powv * pFq(α, β, argument)
            total += term
        end
    end
    return total
end

"""
    _signed_argument(z, exponent::Integer) -> Number

Return `-z` if `exponent` is odd, otherwise `z`.

Applies the sign convention for the hypergeometric argument in the Slater expansions
(DLMF 16.17.2–3).
"""
_signed_argument(z, exponent::Integer) = isodd(exponent) ? -z : z

"""
    _gamma_prod_ratio(nums, denoms, T::Type) -> Number

Compute `∏Γ(nᵢ) / ∏Γ(dⱼ)` in log-space to avoid intermediate overflow.

When the accumulator type `T` is a real floating type (`T <: Real`) this method
accumulates `log|Γ|` and a separate sign using `logabsgamma`, then returns
`sign * exp(log_ratio)`. For other accumulator types (e.g. complex) the routine
uses `loggamma` to preserve phase information and returns `exp(log_ratio)`.
Returns `1` when both iterables are empty.
"""
function _gamma_prod_ratio(nums, denoms, ::Type{T}) where {T<:Real}
    log_ratio = zero(T)
    sign = one(T)
    num_inf = false
    denom_inf = false
    for v in nums
        (lv, sv) = logabsgamma(v)
        if !isfinite(lv)
            num_inf = true
        end
        log_ratio += lv
        sign *= sv
    end
    for v in denoms
        (lv, sv) = logabsgamma(v)
        if !isfinite(lv)
            denom_inf = true
        end
        log_ratio -= lv
        sign *= sv
    end
    if num_inf && denom_inf
        throw(ErrorException("Both the numerator and the denominator of the computed gamma product contain poles. Pole cancellation / analytic limiting behaviour not yet implemented."))
    end
    return sign * exp(log_ratio)
end

function _gamma_prod_ratio(nums, denoms, ::Type{T}) where {T} # Complex or other non-real accumulator types
    log_ratio = zero(T)
    sign = one(T)
    num_inf = false
    denom_inf = false
    for v in nums
        if imag(v) == 0 && real(v) <= 0
            lv, sv = logabsgamma(real(v))
            if !isfinite(lv)
                num_inf = true
            end
            log_ratio += lv
            sign *= sv
        else
            log_ratio += loggamma(v)
        end
    end
    for v in denoms
        if imag(v) == 0 && real(v) <= 0
            lv, sv = logabsgamma(real(v))
            if !isfinite(lv)
                denom_inf = true
            end
            log_ratio -= lv
            sign *= sv
        else
            log_ratio -= loggamma(v)
        end
    end
    if num_inf && denom_inf
        throw(ErrorException("Both the numerator and the denominator of the computed gamma product contain poles. Pole cancellation / analytic limiting behaviour not yet implemented."))
    end
    return sign * exp(log_ratio)
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
    _isintegerlike(x) -> Bool

Return `true` if `x` represents an integer value, regardless of its numeric type.
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
    return _isintegerlike(x) && real(x) > zero(real(x))
end

"""
    _meijerg_slater_boundary_bigfloat(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Evaluate the Slater expansion using internal `BigFloat` arithmetic and convert the
result back to the promoted output type of the original inputs.

This path is used only for balanced cases near `|z| = 1`, where the underlying
hypergeometric evaluations are noticeably more accurate in higher precision.
"""
function _meijerg_slater_boundary_bigfloat(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    promoted = promote(z, a..., b...)
    target_type = eltype(promoted)
    a_big = big.(a)
    b_big = big.(b)
    z_big = big(z)
    mode = _expansion_mode(p, q, z_big)
    value_big = mode === :lower ? _lower_expansion(a_big, b_big, m, n, z_big) : _upper_expansion(a_big, b_big, m, n, z_big)
    return convert(target_type, value_big)
end

