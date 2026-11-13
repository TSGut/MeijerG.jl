const ParameterInput = Union{Tuple, AbstractVector}
const _SLATER_BOUNDARY_BIGFLOAT_RADIUS = 0.02

@doc raw"""
    meijerg_slater(a, b, m, n, z)
    meijerg_slater(a_left, a_right, b_left, b_right, z)

Compute the Meijer G-function using Slater's residue expansions (DLMF §16.17).

This is the **pure Slater evaluator** for advanced users who want direct residue-based computation
without automatic reduction to special functions or confluence handling.

# Overview

The Meijer G-function is defined via a Mellin–Barnes contour integral whose residues yield
finite sums of generalized hypergeometric functions (pFq). This implementation:
- Uses the lower expansion (Slater, DLMF 16.17.2) for p < q or p = q, |z| ≤ 1
- Uses the upper expansion (Slater, DLMF 16.17.3) for p > q or p = q, |z| > 1
- Promotes internally to `BigFloat` near `|z| = 1` for balanced cases to improve accuracy
- Supports arbitrary-precision arithmetic via BigFloat and complex numbers
- **Does NOT** apply reductions to special functions
- **Does NOT** use log-ratio jets or Laurent residues for confluent poles
- **Does NOT** certify the cancellation between its residue terms: off the convergent
  regime of a branch (`|z| > 1` lower, `|z| < 1` upper) the returned value can be
  inaccurate by the factor `κ_sum = Σ|T| / |Σ T|`. `meijerg` applies that certificate
  and a doubled-precision re-evaluation on top of this expansion
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
    # Otherwise determine expansion mode and dispatch
    value = if p == q && abs(abs(z) - one(abs(z))) <= _SLATER_BOUNDARY_BIGFLOAT_RADIUS
        _meijerg_slater_boundary_bigfloat(a, b, m, n, z)
    elseif _expansion_mode(p, q, z) === :lower
        _lower_expansion(a, b, m, n, z)
    else
        _upper_expansion(a, b, m, n, z)
    end

    # Float64 powers of negative reals throw, BigFloat ones return NaN. Make both throw
    z isa Real && z < zero(z) && isnan(value) &&
        throw(DomainError(z, "meijerg_slater on the negative real axis requires an explicit complex argument; pass z + 0im."))
    return _output_like(value, a, b, z)
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
    _output_like(value, a, b, z)

Convert an evaluator result to the type implied by the inputs: the promoted
floating-point type of `a`, `b` and `z`, real when all of them are real. Paths that
work in complex arithmetic (log-ratio jets, contour quadrature) would otherwise
return a complex value with a round-off imaginary part for real input, and
integer input could produce an integer zero.
"""
function _output_like(value, a::Tuple, b::Tuple, z)
    T = float(promote_type(typeof(z), map(typeof, a)..., map(typeof, b)...))
    return T <: Real ? convert(T, real(value)) : convert(T, value)
end

"""
    _lower_term_coeff(ap, bp, m, n, k, T::Type)

Gamma product-ratio coefficient of the `k`-th lower-expansion residue term,
computed with explicit indexed loops in the same evaluation order, sign and
pole handling as the equivalent `_gamma_prod_ratio` call over flattened
generators (which it replaces for performance on this hot path).
"""
function _lower_term_coeff(ap::NTuple{P,T}, bp::NTuple{Q,T}, m::Integer, n::Integer,
                           k::Integer, ::Type{T}, regularize_beta_poles::Bool=false) where {P,Q,T}
    log_ratio = zero(T)
    sign = one(T)
    num_inf = false
    denom_inf = false
    b_k = bp[k]
    @inbounds for j in 1:m
        j == k && continue
        lv, sv, bad = _loggamma_factor_eval(bp[j] - b_k, T)
        log_ratio += lv; sign *= sv; num_inf |= bad
    end
    @inbounds for j in 1:n
        lv, sv, bad = _loggamma_factor_eval(one(T) + b_k - ap[j], T)
        log_ratio += lv; sign *= sv; num_inf |= bad
    end
    @inbounds for j in (m + 1):Q
        v = one(T) + b_k - bp[j]
        regularize_beta_poles && _isnonpositiveintegerlike(v) && continue
        lv, sv, bad = _loggamma_factor_eval(v, T)
        log_ratio -= lv; sign *= sv; denom_inf |= bad
    end
    @inbounds for j in (n + 1):P
        lv, sv, bad = _loggamma_factor_eval(ap[j] - b_k, T)
        log_ratio -= lv; sign *= sv; denom_inf |= bad
    end
    num_inf && denom_inf && _gamma_prod_pole_error()
    return sign * exp(log_ratio)
end

"""
    _upper_term_coeff(ap, bp, m, n, h, T::Type)

Upper-expansion analogue of `_lower_term_coeff` for the `h`-th residue term.
"""
function _upper_term_coeff(ap::NTuple{P,T}, bp::NTuple{Q,T}, m::Integer, n::Integer,
                           h::Integer, ::Type{T}, regularize_beta_poles::Bool=false) where {P,Q,T}
    log_ratio = zero(T)
    sign = one(T)
    num_inf = false
    denom_inf = false
    a_h = ap[h]
    @inbounds for j in 1:n
        j == h && continue
        lv, sv, bad = _loggamma_factor_eval(a_h - ap[j], T)
        log_ratio += lv; sign *= sv; num_inf |= bad
    end
    @inbounds for j in 1:m
        lv, sv, bad = _loggamma_factor_eval(one(T) - a_h + bp[j], T)
        log_ratio += lv; sign *= sv; num_inf |= bad
    end
    @inbounds for j in (n + 1):P
        v = one(T) - a_h + ap[j]
        regularize_beta_poles && _isnonpositiveintegerlike(v) && continue
        lv, sv, bad = _loggamma_factor_eval(v, T)
        log_ratio -= lv; sign *= sv; denom_inf |= bad
    end
    @inbounds for j in (m + 1):Q
        lv, sv, bad = _loggamma_factor_eval(a_h - bp[j], T)
        log_ratio -= lv; sign *= sv; denom_inf |= bad
    end
    num_inf && denom_inf && _gamma_prod_pole_error()
    return sign * exp(log_ratio)
end

"""
    _regularized_pFq_at_poles(alpha, beta, x, regularized)

Evaluate `pFq(alpha, beta, x) / prod(gamma(beta[j]))` over the denominator
indices marked by `regularized`, assuming those parameters are nonpositive
integers. If `K = 1 + max(-beta[j])` over the marked indices, the defining
series starts at `K`. Shifting that series gives an ordinary hypergeometric
function with finite denominator parameters. This evaluates the exact Slater
residue limit without perturbing any Meijer-G parameter.
"""
function _regularized_pFq_at_poles(alpha::Tuple, beta::Tuple, x,
                                    regularized::Tuple)
    length(beta) == length(regularized) ||
        throw(ArgumentError("regularization mask must match denominator parameters"))

    K = 0
    @inbounds for j in eachindex(beta)
        if regularized[j]
            _isnonpositiveintegerlike(beta[j]) ||
                throw(ArgumentError("regularized denominator parameters must be nonpositive integers"))
            K = max(K, 1 - Int(round(real(beta[j]))))
        end
    end
    K > 0 || return pFq(alpha, beta, x)

    # Coefficient of x^K in the regularized series. For a marked beta[j],
    # 1 / (Gamma(beta[j]) * (beta[j])_K) is 1 / Gamma(beta[j] + K).
    coefficient = one(x)
    @inbounds for r in 0:(K - 1)
        for a_i in alpha
            coefficient *= a_i + r
        end
        for j in eachindex(beta)
            regularized[j] || (coefficient /= beta[j] + r)
        end
        coefficient *= x / (r + 1)
    end
    @inbounds for j in eachindex(beta)
        regularized[j] && (coefficient /= gamma(beta[j] + K))
    end
    iszero(coefficient) && return coefficient

    shifted_alpha = (one(x), ((a_i + K) for a_i in alpha)...)
    shifted_beta = (one(x) + K, ((b_j + K) for b_j in beta)...)
    return coefficient * pFq(shifted_alpha, shifted_beta, x)
end

"""
    _pfq_absolute_sum(α, β, x)

Return `Σ_k |(α)_k / (β)_k| |x|^k / k!`, the absolute-value majorant of the
Maclaurin series of `pFq(α; β; x)`, in the real type of `x`.  It feeds the
summation-conditioning certificate of the Slater expansion, where cancellation
inside a family's hypergeometric series is not visible from the family's value.
"""
function _pfq_absolute_sum(α::Tuple, β::Tuple, x; maxterms::Integer=100_000)
    ax = abs(x)
    t = one(ax)
    s = one(ax)
    for k in 0:(maxterms - 1)
        r = ax / (k + 1)
        for aj in α
            r *= abs(aj + k)
        end
        for bj in β
            r /= abs(bj + k)
        end
        t *= r
        s += t
        (iszero(t) || !isfinite(s) || (r < 1 && t <= 1e-17 * s)) && break
    end
    return s
end

"""
    _lower_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Compute the Slater lower expansion (DLMF 16.17.2).

Sums `m` residue terms, each a generalized hypergeometric function `pFq` evaluated
at `±z`. Used when `p < q`, or when `p == q` and `|z| ≤ 1`.
`_lower_expansion_certified` also returns the absolute sum `Σ|T|` over all
residue terms (families times their hypergeometric majorants), from which
`meijerg` forms the summation conditioning `κ_sum = Σ|T| / |Σ T|`.
"""
_lower_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z) =
    _lower_expansion_certified(a, b, m, n, z).value

function _lower_expansion_certified(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    m == 0 && return (value=zero(z), absolute_sum=zero(abs(z)))
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    T = eltype(b_promoted)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(z_promoted, p - m - n)
    total = zero(z_promoted)
    absum = zero(abs(total))
    for k in 1:m
        @inbounds begin
            b_k = b_promoted[k]
            α = ntuple(j -> one(T) + b_k - a_promoted[j], p)
            β = ntuple(q - 1) do idx
                j = idx < k ? idx : idx + 1
                one(T) + b_k - b_promoted[j]
            end
            regularized = ntuple(q - 1) do idx
                j = idx < k ? idx : idx + 1
                j > m && _isnonpositiveintegerlike(β[idx])
            end
            has_regularized = any(regularized)
            term = _lower_term_coeff(a_promoted, b_promoted, m, n, k, T,
                                     has_regularized)

            # A zero unrelated to a hypergeometric denominator pole is a true
            # zero residue. Marked denominator poles are regularized below.
            if iszero(term)
                continue
            end

            powv = z_promoted ^ b_k
            if iszero(powv)
                continue
            end

            hypergeom = has_regularized ?
                _regularized_pFq_at_poles(α, β, argument, regularized) :
                pFq(α, β, argument)
            scale = abs(term * powv)
            term *= powv * hypergeom
            total += term
            absum += has_regularized ? abs(term) : scale * _pfq_absolute_sum(α, β, argument)
        end
    end
    return (value=total, absolute_sum=absum)
end

"""
    _upper_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Compute the Slater upper expansion (DLMF 16.17.3).

Sums `n` residue terms, each a generalized hypergeometric function `pFq` evaluated
at `±1/z`. Used when `p > q`, or when `p == q` and `|z| > 1`.
`_upper_expansion_certified` also returns the absolute sum, as for the lower
expansion.
"""
_upper_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z) =
    _upper_expansion_certified(a, b, m, n, z).value

function _upper_expansion_certified(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    n == 0 && return (value=zero(z), absolute_sum=zero(abs(z)))
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    T = eltype(a_promoted)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(inv(z_promoted), q - m - n)
    total = zero(z_promoted)
    absum = zero(abs(total))
    for h in 1:n
        @inbounds begin
            a_h = a_promoted[h]
            α = ntuple(j -> one(T) - a_h + b_promoted[j], q)
            β = ntuple(p - 1) do idx
                j = idx < h ? idx : idx + 1
                one(T) - a_h + a_promoted[j]
            end
            regularized = ntuple(p - 1) do idx
                j = idx < h ? idx : idx + 1
                j > n && _isnonpositiveintegerlike(β[idx])
            end
            has_regularized = any(regularized)
            term = _upper_term_coeff(a_promoted, b_promoted, m, n, h, T,
                                     has_regularized)

            # A zero unrelated to a hypergeometric denominator pole is a true
            # zero residue. Marked denominator poles are regularized below.
            if iszero(term)
                continue
            end

            powv = z_promoted^(a_h - one(T))
            if iszero(powv)
                continue
            end

            hypergeom = has_regularized ?
                _regularized_pFq_at_poles(α, β, argument, regularized) :
                pFq(α, β, argument)
            scale = abs(term * powv)
            term *= powv * hypergeom
            total += term
            absum += has_regularized ? abs(term) : scale * _pfq_absolute_sum(α, β, argument)
        end
    end
    return (value=total, absolute_sum=absum)
end

"""
    _meijerg_slater_certified(a, b, m, n, z)

Separated Slater expansion of validated inputs, returning `(value, kappa)` with
the summation conditioning `κ_sum = Σ|T| / |Σ T|` of the residue sum.  The
balanced boundary path (`p == q`, `|z| ≈ 1`) works in BigFloat internally and
reports `κ_sum = 1`.
"""
function _meijerg_slater_certified(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    if p == q && abs(abs(z) - one(abs(z))) <= _SLATER_BOUNDARY_BIGFLOAT_RADIUS
        return (value=_meijerg_slater_boundary_bigfloat(a, b, m, n, z), kappa=1.0)
    end
    r = _expansion_mode(p, q, z) === :lower ? _lower_expansion_certified(a, b, m, n, z) :
                                               _upper_expansion_certified(a, b, m, n, z)
    kappa = iszero(r.value) ? (iszero(r.absolute_sum) ? 1.0 : Inf) :
                              Float64(r.absolute_sum / abs(r.value))
    return (value=r.value, kappa=kappa)
end

"""
    _signed_argument(z, exponent::Integer) -> Number

Return `-z` if `exponent` is odd, otherwise `z`.

Applies the sign convention for the hypergeometric argument in the Slater expansions
(DLMF 16.17.2–3).
"""
_signed_argument(z, exponent::Integer) = isodd(exponent) ? -z : z

"""
    _loggamma_factor_eval(v, T::Type) -> (lv, sv, ispole)

Evaluate one gamma factor in log space, mirroring the two `_gamma_prod_ratio`
methods: `logabsgamma` with explicit sign for real accumulator types and for
nonpositive real arguments under non-real accumulator types, `loggamma`
elsewhere.  Shared by the Slater residue coefficients and the Mellin-Barnes
contour kernel.  The caller accumulates `log_ratio ± lv` and `sign * sv`
directly. No intermediate sum may be introduced because
`loggamma(::Complex{BigFloat})` carries extended internal precision and an
intermediate `zero(T) + lv` would double-round.
"""
@inline function _loggamma_factor_eval(v, ::Type{T}) where {T<:Real}
    lv, sv = logabsgamma(v)
    return lv, sv, !isfinite(lv)
end

@inline function _loggamma_factor_eval(v, ::Type{T}) where {T}
    if imag(v) == 0 && real(v) <= 0
        lv, sv = logabsgamma(real(v))
        return lv, sv, !isfinite(lv)
    end
    return loggamma(v), true, false
end

@inline function _gamma_prod_pole_error()
    throw(ErrorException("Both the numerator and the denominator of the computed gamma product contain poles. Pole cancellation / analytic limiting behaviour not yet implemented."))
end

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
    num_inf && denom_inf && _gamma_prod_pole_error()
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
    num_inf && denom_inf && _gamma_prod_pole_error()
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

_isnonpositiveintegerlike(x) = _isintegerlike(x) && real(x) <= zero(real(x))

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
