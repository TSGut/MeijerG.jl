# ConfluentJets: type-generic truncated Taylor jets and stable confluent divided differences.
# A self-contained submodule depending only on LinearAlgebra and Gamma.

module ConfluentJets

using LinearAlgebra
import Gamma
import Gamma: gamma, loggamma, logabsgamma, digamma

export Jet, jet_order, coefficient, coefficients,
       constant_jet, variable_jet, affine_jet, affine_log_jet, compose_scale,
       gamma, loggamma, digamma, polygamma, reciprocal_gamma,
       complete_homogeneous, taylor_matrix,
       exp_divided_difference, exp_product_divided_difference,
       default_opitz_pad,
       exprel, safe_log1p, scaled_exprel, finite_number,
       real_float_type, integer_as,
       logpochhammer, loggamma_divdiff_limit,
       logpochhammer_divdiff_limit, loggamma_divdiff,
       loggamma_divdiff_pair, logpochhammer_divdiff



# ======================================================================
# ---- jet_algebra.jl ----
# ======================================================================

"""
    Jet(coefficients)

A fixed-order truncated power series
`sum(coefficients[k+1] * t^k)`, computed modulo the first omitted power.
"""
struct Jet{T<:Number}
    coeffs::Vector{T}
end

Jet(coeffs::AbstractVector{T}) where {T<:Number} = Jet{T}(collect(coeffs))

"""Return a jet's truncation order."""
jet_order(a::Jet) = length(a.coeffs) - 1

"""Return a copy of the stored Taylor coefficients."""
coefficients(a::Jet) = copy(a.coeffs)

"""Return the coefficient of `t^k`, or zero outside the stored order."""
coefficient(a::Jet, k::Integer) =
    0 <= k <= jet_order(a) ? a.coeffs[k + 1] : zero(eltype(a.coeffs))

"""Construct a constant jet of order `N`."""
function constant_jet(c::T, N::Integer) where {T<:Number}
    values = zeros(T, N + 1)
    values[1] = c
    return Jet(values)
end

"""Construct the variable jet `t` with coefficient type `T` and order `N`."""
function variable_jet(::Type{T}, N::Integer) where {T<:Number}
    values = zeros(T, N + 1)
    N >= 1 && (values[2] = one(T))
    return Jet(values)
end

"""Construct the affine jet `constant + slope*t` to order `N`."""
function affine_jet(constant::T, slope::T, N::Integer) where {T<:Number}
    # Same coefficients as `constant_jet(constant, N) + slope * variable_jet(T, N)`
    # (every omitted operation adds or scales an exact zero), without temporaries.
    values = zeros(T, N + 1)
    values[1] = constant
    N >= 1 && (values[2] = slope)
    return Jet(values)
end

"""
    affine_log_jet(constant, slope, N)

Return the Taylor jet of `log(constant + slope*t)` through order `N`.
This closed form preserves the operation order used by the extracted
confluent divided-difference kernel.
"""
function affine_log_jet(constant::T, slope::T,
                        N::Integer) where {T<:Number}
    iszero(constant) &&
        error("affine_log_jet: zero constant term")
    values = Vector{T}(undef, N + 1)
    values[1] = log(constant)
    for p in 1:N
        values[p + 1] = -(-slope / constant)^p / p
    end
    return Jet(values)
end

function _promote_jets(a::Jet, b::Jet)
    return (N=min(jet_order(a), jet_order(b)),
            T=promote_type(eltype(a.coeffs), eltype(b.coeffs)))
end

function Base.:+(a::Jet, b::Jet)
    nt = _promote_jets(a, b)
    return Jet([a.coeffs[k + 1] + b.coeffs[k + 1] for k in 0:nt.N])
end

Base.:-(a::Jet) = Jet([-x for x in a.coeffs])
Base.:-(a::Jet, b::Jet) = a + (-b)
Base.:*(s::Number, a::Jet) = Jet([s * x for x in a.coeffs])
Base.:*(a::Jet, s::Number) = s * a
Base.:/(a::Jet, s::Number) = Jet([x / s for x in a.coeffs])

function Base.:+(a::Jet, s::Number)
    values = copy(a.coeffs)
    values[1] += s
    return Jet(values)
end

Base.:+(s::Number, a::Jet) = a + s
Base.:-(a::Jet, s::Number) = a + (-s)
Base.:-(s::Number, a::Jet) = s + (-a)

function Base.:*(a::Jet, b::Jet)
    nt = _promote_jets(a, b)
    values = zeros(nt.T, nt.N + 1)
    @inbounds for i in 0:nt.N, j in 0:(nt.N - i)
        values[i + j + 1] += a.coeffs[i + 1] * b.coeffs[j + 1]
    end
    return Jet(values)
end

function Base.inv(b::Jet)
    b0 = b.coeffs[1]
    iszero(b0) && error("inv(::Jet): zero constant term")
    N = jet_order(b)
    T = promote_type(eltype(b.coeffs), typeof(inv(b0)))
    values = zeros(T, N + 1)
    values[1] = inv(b0)
    for n in 1:N
        total = zero(T)
        @inbounds for k in 1:n
            total += b.coeffs[k + 1] * values[n - k + 1]
        end
        values[n + 1] = -total / b0
    end
    return Jet(values)
end

Base.:/(a::Jet, b::Jet) = a * inv(b)

function Base.:^(a::Jet, n::Integer)
    n < 0 && return inv(a^(-n))
    result = constant_jet(one(eltype(a.coeffs)), jet_order(a))
    base = a
    exponent = n
    while exponent > 0
        isodd(exponent) && (result = result * base)
        exponent >>= 1
        exponent > 0 && (base = base * base)
    end
    return result
end

Base.zero(a::Jet) = constant_jet(zero(eltype(a.coeffs)), jet_order(a))
Base.one(a::Jet) = constant_jet(one(eltype(a.coeffs)), jet_order(a))
Base.iszero(a::Jet) = all(iszero, a.coeffs)

function Base.exp(a::Jet)
    N = jet_order(a)
    T = promote_type(eltype(a.coeffs), typeof(exp(a.coeffs[1])))
    values = zeros(T, N + 1)
    values[1] = exp(a.coeffs[1])
    for n in 1:N
        total = zero(T)
        @inbounds for k in 1:n
            total += k * a.coeffs[k + 1] * values[n - k + 1]
        end
        values[n + 1] = total / n
    end
    return Jet(values)
end

# Log of the affine jet `c0 + slope*t` into `values`: the general recurrence below restricted
# to its nonzero terms, so the result is bit-identical at O(N) instead of O(N^2).
function _log_affine_into!(values::AbstractVector{T}, c0, slope, N::Integer) where {T}
    values[1] = log(c0)
    N >= 1 || return values
    values[2] = (1 * slope) / (1 * c0)
    @inbounds for n in 2:N
        total = n * zero(slope)
        total -= slope * (n - 1) * values[n]
        values[n + 1] = total / (n * c0)
    end
    return values
end

# Orders 1..N of `_log_affine_into!` without the order-0 logarithm (the recurrence only
# reads `values[2:N]`), for callers that discard the constant term.
function _log_affine_derivatives_into!(values::AbstractVector{T}, c0, slope, N::Integer) where {T}
    values[1] = zero(T)
    N >= 1 || return values
    values[2] = (1 * slope) / (1 * c0)
    @inbounds for n in 2:N
        total = n * zero(slope)
        total -= slope * (n - 1) * values[n]
        values[n + 1] = total / (n * c0)
    end
    return values
end

function Base.log(a::Jet)
    a0 = a.coeffs[1]
    iszero(a0) && error("log(::Jet): zero constant term")
    N = jet_order(a)
    T = promote_type(eltype(a.coeffs), typeof(log(a0)))
    values = zeros(T, N + 1)
    affine = true
    @inbounds for j in 3:length(a.coeffs)
        affine &= iszero(a.coeffs[j])
    end
    if affine
        _log_affine_into!(values, a0, N >= 1 ? a.coeffs[2] : zero(a0), N)
        return Jet(values)
    end
    values[1] = log(a0)
    for n in 1:N
        total = n * a.coeffs[n + 1]
        @inbounds for k in 1:(n - 1)
            total -= a.coeffs[k + 1] * (n - k) *
                     values[n - k + 1]
        end
        values[n + 1] = total / (n * a0)
    end
    return Jet(values)
end

Base.log1p(a::Jet) = log(one(a) + a)
Base.expm1(a::Jet) = exp(a) - one(a)

"""Substitute `t -> scale*t` in a jet."""
compose_scale(a::Jet, scale::Number) =
    Jet([a.coeffs[k + 1] * scale^k for k in 0:jet_order(a)])

# ======================================================================
# ---- gamma_extensions.jl ----
# ======================================================================

function _bernoulli_even(count::Integer)
    count >= 1 || return Rational{BigInt}[]
    work = Vector{Rational{BigInt}}(undef, 2count + 1)
    result = Rational{BigInt}[]
    for m in 0:2count
        work[m + 1] = BigInt(1) // BigInt(m + 1)
        for j in m:-1:1
            work[j] = j * (work[j] - work[j + 1])
        end
        iseven(m) && m > 0 && push!(result, work[1])
    end
    return result
end

const _LOGGAMMA_SHIFT_THRESHOLD = 6
const _LOGGAMMA_STIRLING_TERMS = 16
const _LOGGAMMA_BERNOULLI = Ref(_bernoulli_even(_LOGGAMMA_STIRLING_TERMS))
const _LOGGAMMA_BERNOULLI_LOCK = ReentrantLock()

function _loggamma_bernoulli(count::Integer)
    length(_LOGGAMMA_BERNOULLI[]) >= count &&
        return _LOGGAMMA_BERNOULLI[]
    lock(_LOGGAMMA_BERNOULLI_LOCK) do
        if length(_LOGGAMMA_BERNOULLI[]) < count
            _LOGGAMMA_BERNOULLI[] = _bernoulli_even(count)
        end
        return _LOGGAMMA_BERNOULLI[]
    end
end

function _loggamma_parameters(::Type{R}) where {R<:Real}
    if R === BigFloat
        bits = precision(BigFloat)
        threshold = ceil(Int, bits * log(2) / (2pi) + 24)
        return threshold, ceil(Int, pi * threshold)
    elseif R <: Union{Float16,Float32,Float64}
        return _LOGGAMMA_SHIFT_THRESHOLD, _LOGGAMMA_STIRLING_TERMS
    end
    # Other fixed-precision real types (e.g. a double-double): shift until the
    # Stirling remainder after `terms` terms is below the type's unit roundoff.
    bits = try
        precision(R)
    catch
        53
    end
    threshold = ceil(Int, bits * log(2) / (2pi) + 12)
    return threshold, ceil(Int, 0.7 * threshold)
end

function _rational_as(::Type{T}, value::Rational{BigInt}) where {T<:Number}
    R = typeof(real(zero(T)))
    return T(convert(R, numerator(value)) / convert(R, denominator(value)))
end

function _affine_negative_power(
        constant::T, slope::T, m::Integer,
        N::Integer) where {T<:Number}
    m >= 1 || throw(DomainError(m))
    values = Vector{T}(undef, N + 1)
    values[1] = constant^(-m)
    ratio = -slope / constant
    for k in 1:N
        values[k + 1] =
            values[k] * (T(m + k - 1) / T(k)) * ratio
    end
    return Jet(values)
end

"""
    _jet_fast_generic(::Type{R})

Trait: `true` for fixed-precision real types above Float64 (such as a
double-double) whose log-gamma jets take the allocation-light affine path
`_loggamma_jet_affine_fast`.  Float64 and BigFloat keep the reference path.
"""
_jet_fast_generic(::Type) = false

const _LOGGAMMA_SCALAR_CACHE = Dict{Tuple{DataType,Int},Vector}()

# Stirling scalars B_{2k} / (2k(2k-1)), k = 1..terms, converted once per type.
function _loggamma_stirling_scalars(::Type{T}, terms::Integer) where {T}
    key = (T, Int(terms))
    lock(_LOGGAMMA_BERNOULLI_LOCK) do
        get!(_LOGGAMMA_SCALAR_CACHE, key) do
            bernoulli = _loggamma_bernoulli(terms)
            T[_rational_as(T, bernoulli[k]) / T(2k * (2k - 1)) for k in 1:terms]
        end
    end::Vector{T}
end

# The same scalars for the hardware types at the fixed hardware term count,
# built once at load so the hardware jet path needs no lock.
_hardware_stirling_table(::Type{T}) where {T} =
    T[_rational_as(T, _LOGGAMMA_BERNOULLI[][k]) / T(2k * (2k - 1))
      for k in 1:_LOGGAMMA_STIRLING_TERMS]
const _LOGGAMMA_STIRLING_F64 = _hardware_stirling_table(Float64)
const _LOGGAMMA_STIRLING_CF64 = _hardware_stirling_table(ComplexF64)
const _LOGGAMMA_STIRLING_F32 = _hardware_stirling_table(Float32)
const _LOGGAMMA_STIRLING_CF32 = _hardware_stirling_table(ComplexF32)
_loggamma_stirling_hardware(::Type{Float64}) = _LOGGAMMA_STIRLING_F64
_loggamma_stirling_hardware(::Type{ComplexF64}) = _LOGGAMMA_STIRLING_CF64
_loggamma_stirling_hardware(::Type{Float32}) = _LOGGAMMA_STIRLING_F32
_loggamma_stirling_hardware(::Type{ComplexF32}) = _LOGGAMMA_STIRLING_CF32

const _HardwareJetTypes = Union{Float32,Float64,ComplexF32,ComplexF64}

# Argument-independent ratios `T(m + j - 1) / T(j)` of the Stirling term `w^-m`, `m = 2k - 1`,
# tabulated per type as `[k][j]` for `j <= 64` (same expression, evaluated once).
const _LOGGAMMA_RATIO_ORDERS = 64
_hardware_ratio_table(::Type{T}) where {T} =
    [T[T(2k - 1 + j - 1) / T(j) for j in 1:_LOGGAMMA_RATIO_ORDERS]
     for k in 1:_LOGGAMMA_STIRLING_TERMS]
const _LOGGAMMA_RATIOS_F64 = _hardware_ratio_table(Float64)
const _LOGGAMMA_RATIOS_CF64 = _hardware_ratio_table(ComplexF64)
const _LOGGAMMA_RATIOS_F32 = _hardware_ratio_table(Float32)
const _LOGGAMMA_RATIOS_CF32 = _hardware_ratio_table(ComplexF32)
_loggamma_ratios_hardware(::Type{Float64}) = _LOGGAMMA_RATIOS_F64
_loggamma_ratios_hardware(::Type{ComplexF64}) = _LOGGAMMA_RATIOS_CF64
_loggamma_ratios_hardware(::Type{Float32}) = _LOGGAMMA_RATIOS_F32
_loggamma_ratios_hardware(::Type{ComplexF32}) = _LOGGAMMA_RATIOS_CF32

# Log-gamma jet of an affine input for the hardware float types.  Orders 0 and
# 1 come from the scalar log-gamma and digamma.
function _loggamma_jet_affine_hardware(z::Jet{T}, threshold_integer::Integer,
                                       terms::Integer) where {T<:_HardwareJetTypes}
    N = jet_order(z)
    R = typeof(real(zero(T)))
    z0 = z.coeffs[1]
    values = Vector{T}(undef, N + 1)
    if N >= 2
        slope = z.coeffs[2]
        threshold = convert(R, threshold_integer)
        shift = max(0, ceil(Int, threshold - real(z0)))
        w0 = z0 + T(shift)
        half = T(one(R) / 2)
        a0 = w0 - half
        # (w - 1/2) log w, orders >= 2, with log w by the affine recurrence.
        logterm = Vector{T}(undef, N + 1)
        _log_affine_derivatives_into!(logterm, w0, slope, N)
        @inbounds for n in 2:N
            values[n + 1] = a0 * logterm[n + 1] + slope * logterm[n]
        end
        compensation = zeros(T, N + 1)
        scalars = terms == _LOGGAMMA_STIRLING_TERMS ?
            _loggamma_stirling_hardware(T) : _loggamma_stirling_scalars(T, terms)
        power = Vector{T}(undef, N + 1)
        ratio = -slope / w0
        inverse = inv(w0)
        tabulated = terms == _LOGGAMMA_STIRLING_TERMS && N <= _LOGGAMMA_RATIO_ORDERS
        ratios = _loggamma_ratios_hardware(T)
        for k in 1:terms
            m = 2k - 1
            scalar = scalars[k]
            # Base evaluates complex `w0^(-m)` as `inv(w0)^m`, real floats use a compensated power.
            power[1] = T <: Complex ? inverse^m : w0^(-m)
            if tabulated
                row = ratios[k]
                @inbounds for j in 1:N
                    power[j + 1] = power[j] * row[j] * ratio
                end
            else
                @inbounds for j in 1:N
                    power[j + 1] = power[j] * (T(m + j - 1) / T(j)) * ratio
                end
            end
            @inbounds for j in 3:(N + 1)
                y = scalar * power[j] - compensation[j]
                updated = values[j] + y
                compensation[j] = (updated - values[j]) - y
                values[j] = updated
            end
        end
        for k in 0:(shift - 1)
            _log_affine_derivatives_into!(logterm, z0 + T(k), slope, N)
            @inbounds for j in 3:(N + 1)
                y = -logterm[j] - compensation[j]
                updated = values[j] + y
                compensation[j] = (updated - values[j]) - y
                values[j] = updated
            end
        end
    end
    values[1] = convert(T, loggamma(z0))
    N >= 1 && (values[2] = convert(T, digamma(z0)) * z.coeffs[2])
    return Jet(values)
end

# Log-gamma jet of an affine input for fast generic types: the reference shift-plus-Stirling
# evaluation with closed-form affine logs, one log of the shift product and no intermediate jets.
function _loggamma_jet_affine_fast(z::Jet{T}, threshold::Integer, terms::Integer) where {T}
    N = jet_order(z)
    R = typeof(real(zero(T)))
    z0 = z.coeffs[1]
    slope = N >= 1 ? z.coeffs[2] : zero(T)
    shift = max(0, ceil(Int, convert(R, threshold) - real(z0)))
    w0 = z0 + T(shift)
    half = T(one(R) / 2)
    values = zeros(T, N + 1)
    # (w - 1/2) log w - w + log(2π)/2 with log w = log w0 + Σ_p -(-slope/w0)^p / p t^p.
    r = -slope / w0
    logw0 = log(w0)
    values[1] = (w0 - half) * logw0 - w0 + T(log(convert(R, 2) * convert(R, pi)) / 2)
    power = one(T)
    prev = logw0                       # coefficient p-1 of log w
    @inbounds for p in 1:N
        power *= r
        lp = -power / R(p)
        values[p + 1] = (w0 - half) * lp + slope * prev - (p == 1 ? slope : zero(T))
        prev = lp
    end
    # Stirling series Σ_k c_k w^{-(2k-1)}. Coefficient j of w^{-m} is
    # w0^{-m} r^j (m)_j / j!, built by the same recurrence as `_affine_negative_power`.
    scalars = _loggamma_stirling_scalars(T, terms)
    cinv = inv(w0)
    cinv2 = cinv * cinv
    cpow = cinv
    @inbounds for k in 1:terms
        m = 2k - 1
        coeff = scalars[k] * cpow
        values[1] += coeff
        for j in 1:N
            coeff = coeff * r * R(m + j - 1) / R(j)
            values[j + 1] += coeff
        end
        cpow = cpow * cinv2
    end
    # Shift: subtract Σ_{k<shift} log(z0 + k + slope t).
    if shift > 0
        prod = one(T)
        @inbounds for k in 0:(shift - 1)
            c0 = z0 + T(k)
            prod *= c0
            if (k + 1) % 32 == 0
                values[1] -= log(prod)
                prod = one(T)
            end
            rk = -slope / c0
            pk = one(T)
            for p in 1:N
                pk *= rk
                values[p + 1] += pk / R(p)
            end
        end
        values[1] -= log(prod)
    end
    N >= 1 && (values[2] = convert(T, digamma(z0)) * slope)
    return Jet(values)
end

function Gamma.loggamma(z::Jet)
    N = jet_order(z)
    T = eltype(z.coeffs)
    R = typeof(real(zero(T)))
    z0 = z.coeffs[1]
    threshold_integer, terms = _loggamma_parameters(R)
    if _jet_fast_generic(R)
        affine = true
        @inbounds for j in 3:length(z.coeffs)
            affine &= iszero(z.coeffs[j])
        end
        affine && return _loggamma_jet_affine_fast(z, threshold_integer, terms)
    end
    if T <: _HardwareJetTypes
        affine = true
        @inbounds for j in 3:length(z.coeffs)
            affine &= iszero(z.coeffs[j])
        end
        affine && return _loggamma_jet_affine_hardware(z, threshold_integer, terms)
    end
    threshold = convert(R, threshold_integer)
    shift = max(0, ceil(Int, threshold - real(z0)))
    w = z + T(shift)
    half = T(one(R) / 2)
    logtwopi_half = T(log(convert(R, 2) * convert(R, pi)) / 2)
    base = (w - half) * log(w) - w + constant_jet(logtwopi_half, N)
    values = copy(base.coeffs)
    compensation = zeros(T, N + 1)
    bernoulli = _loggamma_bernoulli(terms)
    # Hardware floats replace orders 0 and 1 by the scalar log-gamma and digamma below, so
    # their accumulation and the shift logarithms (which only feed order 0) are skipped.
    hardware = T === Float32 || T === Float64 ||
               T === ComplexF32 || T === ComplexF64
    jstart = hardware ? 3 : 1

    affine_input = true
    @inbounds for j in 3:length(w.coeffs)
        affine_input &= iszero(w.coeffs[j])
    end
    if affine_input
        slope = N >= 1 ? w.coeffs[2] : zero(T)
        # Coefficients of w^{-m} formed in a reused buffer by the recurrence of
        # `_affine_negative_power`, with the same operations (no temporary jets).
        power = Vector{T}(undef, N + 1)
        constant = w.coeffs[1]
        ratio = -slope / constant
        for k in 1:terms
            m = 2k - 1
            scalar = _rational_as(T, bernoulli[k]) / T(2k * (2k - 1))
            power[1] = constant^(-m)
            @inbounds for j in 1:N
                power[j + 1] = power[j] * (T(m + j - 1) / T(j)) * ratio
            end
            @inbounds for j in jstart:length(values)
                y = scalar * power[j] - compensation[j]
                updated = values[j] + y
                compensation[j] = (updated - values[j]) - y
                values[j] = updated
            end
        end
    else
        inverse = inv(w)
        inverse2 = inverse * inverse
        power = inverse
        for k in 1:terms
            scalar = _rational_as(T, bernoulli[k]) / T(2k * (2k - 1))
            term = scalar * power
            @inbounds for j in jstart:length(values)
                y = term.coeffs[j] - compensation[j]
                updated = values[j] + y
                compensation[j] = (updated - values[j]) - y
                values[j] = updated
            end
            power = power * inverse2
        end
    end
    if shift > 0 && length(values) >= jstart
        logterm = Vector{T}(undef, N + 1)
        slope = N >= 1 ? z.coeffs[2] : zero(T)
        for k in 0:(shift - 1)
            if affine_input && hardware
                _log_affine_derivatives_into!(logterm, z0 + T(k), slope, N)
            elseif affine_input
                _log_affine_into!(logterm, z0 + T(k), slope, N)
            else
                copyto!(logterm, log(z + T(k)).coeffs)
            end
            @inbounds for j in jstart:length(values)
                y = -logterm[j] - compensation[j]
                updated = values[j] + y
                compensation[j] = (updated - values[j]) - y
                values[j] = updated
            end
        end
    end
    if hardware
        values[1] = convert(T, loggamma(z0))
    end
    N >= 1 && (values[2] = convert(T, digamma(z0)) * z.coeffs[2])
    return Jet(values)
end

Gamma.gamma(z::Jet) = exp(loggamma(z))

function _compose_taylor(values, delta::Jet)
    result = constant_jet(last(values), jet_order(delta))
    for k in (length(values) - 1):-1:1
        result = result * delta + values[k]
    end
    return result
end

function polygamma(m::Integer, z::Jet)
    m >= 0 || throw(DomainError(m, "polygamma order must be nonnegative"))
    N = jet_order(z)
    T = eltype(z.coeffs)
    center = z.coeffs[1]
    source = loggamma(affine_jet(center, one(T), N + m + 1)).coeffs
    values = Vector{eltype(source)}(undef, N + 1)
    for k in 0:N
        factor = one(real(one(eltype(source))))
        for j in 1:(m + 1)
            factor *= k + j
        end
        values[k + 1] = source[k + m + 2] * factor
    end
    delta = z - center
    return _compose_taylor(values, delta)
end

Gamma.digamma(z::Jet) = polygamma(0, z)

"""
    reciprocal_gamma(jet; radius=0.1)

Return a pole-aware reciprocal-gamma jet near a nonpositive integer, or
`nothing` when the constant term is outside `radius`.
"""
function reciprocal_gamma(z::Jet; radius=0.1)
    z0 = z.coeffs[1]
    pole_order = max(0, round(Int, -real(z0)))
    u0 = z0 + pole_order
    abs(u0) <= radius || return nothing
    T = eltype(z.coeffs)
    u = z + T(pole_order)
    numerator_jet = one(z)
    for k in 0:pole_order
        numerator_jet = numerator_jet * (u - T(k))
    end
    return numerator_jet * exp(-loggamma(one(z) + u))
end

# ======================================================================
# ---- opitz.jl ----
# ======================================================================

"""Return the complete homogeneous symmetric polynomial `h_p(vars)`."""
function complete_homogeneous(p::Integer, vars::AbstractVector{T}) where {T}
    p == 0 && return one(T)
    isempty(vars) && return zero(T)
    x = vars[1]
    rest = @view vars[2:end]
    value = zero(T)
    xpower = one(T)
    for j in 0:p
        value += xpower * complete_homogeneous(p - j, rest)
        xpower *= x
    end
    return value
end

"""Evaluate a Taylor jet on the upper-triangular Opitz node matrix."""
# `h[p+1, i, j]` is the complete homogeneous polynomial `h_p(deltas[i:j])`, built from
# `h_p(x_i..x_j) = h_p(x_{i+1}..x_j) + x_i h_{p-1}(x_i..x_j)`.
function taylor_matrix(values::AbstractVector{T},
                       deltas::AbstractVector) where {T}
    r = length(deltas)
    length(values) >= r || error("taylor_matrix: insufficient jet order")
    converted = collect(T, deltas)
    maximum_order = length(values) - 1
    h = zeros(T, maximum_order + 1, r, r)
    @inbounds for i in 1:r, j in i:r
        h[1, i, j] = one(T)
    end
    @inbounds for p in 1:maximum_order, j in 1:r, i in j:-1:1
        h[p + 1, i, j] = (i < j ? h[p + 1, i + 1, j] : zero(T)) + converted[i] * h[p, i, j]
    end
    matrix = zeros(T, r, r)
    @inbounds for i in 1:r, j in i:r
        degree = j - i
        entry = zero(T)
        for m in degree:maximum_order
            entry += values[m + 1] * h[m - degree + 1, i, j]
        end
        matrix[i, j] = entry
    end
    return matrix
end

taylor_matrix(jet::Jet, deltas::AbstractVector) =
    taylor_matrix(jet.coeffs, deltas)

# Coefficient types with a LinearAlgebra matrix exponential. Other types take the jet route,
# since `exp(taylor_matrix(L))` equals `taylor_matrix(exp(L))` up to the truncation order.
const _MatrixExpTypes = Union{Float32,Float64,ComplexF32,ComplexF64}

"""Return the highest divided difference of `exp(logjet)` at node offsets."""
function exp_divided_difference(logjet::Jet{T}, deltas::AbstractVector) where {T}
    T <: _MatrixExpTypes ||
        return taylor_matrix(exp(logjet), deltas)[1, length(deltas)]
    matrix = taylor_matrix(logjet, deltas)
    return exp(matrix)[1, length(deltas)]
end

"""
    exp_product_divided_difference(factorjet, logjet, deltas)

Return the highest divided difference of `factorjet * exp(logjet)`.
"""
function exp_product_divided_difference(factorjet::Jet, logjet::Jet{T},
                                        deltas::AbstractVector) where {T}
    T <: _MatrixExpTypes ||
        return taylor_matrix(factorjet * exp(logjet), deltas)[1, length(deltas)]
    logmatrix = taylor_matrix(logjet, deltas)
    factormatrix = taylor_matrix(factorjet, deltas)
    return (factormatrix * exp(logmatrix))[1, length(deltas)]
end

"""Return the existing diameter-adaptive default Opitz jet pad."""
function default_opitz_pad(deltas; maximum_pad::Integer=16)
    diameter = 0.0
    @inbounds for i in eachindex(deltas), j in (i + 1):lastindex(deltas)
        diameter = max(diameter, abs(deltas[i] - deltas[j]))
    end
    diameter <= 1e-12 && return 6
    return clamp(ceil(Int, 2.5 * log10(diameter) + 21), 6, maximum_pad)
end

# ======================================================================
# ---- scalar_helpers.jl ----
# ======================================================================

"""Compute `(exp(x)-1)/x` with the existing small-argument expansion."""
function exprel(x)
    ax = abs(x)
    ax < 1e-8 && return one(x) + x / 2 + x^2 / 6 + x^3 / 24
    return expm1(x) / x
end

"""Return `log1p(x)`, or NaN for a real argument below `-1`."""
safe_log1p(x::Real) = x >= -1 ? log1p(x) : oftype(float(x), NaN)
safe_log1p(x) = log1p(x)

"""Compute `M*exprel(a*M)` without division by `M`."""
function scaled_exprel(a, M)
    a == zero(a) && return M
    return expm1(a * M) / a
end

"""Return whether both components of a number are finite."""
finite_number(value) =
    value isa Number && isfinite(real(value)) && isfinite(imag(value))

"""Return the floating real component type associated with a number."""
real_float_type(value) = typeof(float(real(value)))

"""Convert an integer to the floating real component type of `value`."""
integer_as(value, n::Integer) = convert(real_float_type(value), n)

"""Return `log((a)_n)` using the scalar log-gamma backend."""
logpochhammer(a, n::Integer) =
    n == 0 ? zero(a + 0im) : loggamma(a + n) - loggamma(a)

"""Return the zero-separation log-gamma divided-difference limit."""
loggamma_divdiff_limit(a, sigma::Integer) = sigma * digamma(a)

"""Return the zero-separation log-Pochhammer divided-difference limit."""
function logpochhammer_divdiff_limit(a, n::Integer, sigma::Integer)
    n == 0 && return zero(a + 0im)
    total = zero(a + 0im)
    for k in 0:(n - 1)
        total += sigma / (a + k)
    end
    return total
end

"""
    loggamma_divdiff(a, epsilon, sigma; tolerance=nothing)

Evaluate `(loggamma(a + sigma*epsilon)-loggamma(a))/epsilon` from a
log-gamma jet, including the native BigFloat order and tolerance policy.
"""
function loggamma_divdiff(a, epsilon, sigma::Integer; tolerance=nothing)
    R = real_float_type(a + epsilon)
    hardware = R <: Union{Float16,Float32,Float64}
    if tolerance === nothing
        tolerance = hardware ? 1e-16 : convert(R, 10) * eps(R)
    end
    Nmax = R === BigFloat ? 19 : 13
    T = a isa Real && a > zero(a) ? typeof(float(a)) :
        typeof(complex(float(a)))
    aT = convert(T, a)
    sT = convert(T, sigma)
    N = _loggamma_divdiff_degree(aT, epsilon, tolerance, Nmax)
    if N < Nmax
        result, converged = _loggamma_divdiff_sum(aT, sT, epsilon, N, tolerance)
        converged && return result
    end
    result, _ = _loggamma_divdiff_sum(aT, sT, epsilon, Nmax, tolerance)
    return result
end

# Degree estimate: log-gamma Taylor coefficients at `a` decay like `d^-r / r`, `d` the distance
# to the nearest pole, so about `log(tolerance) / log(|epsilon| / d)` terms are needed (+2 margin).
function _loggamma_divdiff_degree(a, epsilon, tolerance, Nmax::Integer)
    e = Float64(abs(epsilon))
    (isfinite(e) && e > 0) || return 2
    tol = Float64(tolerance)
    tol > 0 || return Nmax
    re = Float64(real(a))
    k = (-1e9 < re < 0) ? -round(Int, re) : 0
    d = min(1.0, Float64(abs(a + k)))
    ratio = e / d
    (isfinite(ratio) && ratio < 1) || return Nmax
    estimate = log(tol) / log(ratio) + 2
    estimate < Nmax || return Nmax
    return max(2, ceil(Int, estimate))
end

function _loggamma_divdiff_sum(aT, sT, epsilon, N::Integer, tolerance)
    jet = loggamma(affine_jet(aT, sT, N))
    return _loggamma_divdiff_sum_coeffs(jet.coeffs, epsilon, N, tolerance, 1)
end

# Divided-difference sum from the slope-`+1` jet. The slope-`-1` jet has the same coefficients
# with every odd order negated (each jet operation is odd in the slope).
function _loggamma_divdiff_sum_coeffs(coeffs, epsilon, N::Integer, tolerance, sigma::Integer)
    result = sigma == 1 ? coeffs[2] : -coeffs[2]
    epsilon_power = one(epsilon)
    converged = false
    for r in 1:(N - 1)
        epsilon_power *= epsilon
        term = coeffs[r + 2] * epsilon_power
        (sigma == 1 || isodd(r)) || (term = -term)
        result += term
        if abs(term) < tolerance * max(abs(result), one(abs(result)))
            converged = true
            break
        end
    end
    return result, converged
end

"""
    loggamma_divdiff_pair(a, epsilon; tolerance=nothing)

Return `(loggamma_divdiff(a, epsilon, +1), loggamma_divdiff(a, epsilon, -1))`.
For hardware floats both are formed from one log-gamma jet (the slope `-1`
jet is the slope `+1` jet with alternating signs), with the same degree
policy and fallback as the single-sign evaluation, so the results are
identical to two separate calls.
"""
function loggamma_divdiff_pair(a, epsilon; tolerance=nothing)
    R = real_float_type(a + epsilon)
    hardware = R <: Union{Float16,Float32,Float64}
    hardware || return (loggamma_divdiff(a, epsilon, 1; tolerance),
                        loggamma_divdiff(a, epsilon, -1; tolerance))
    if tolerance === nothing
        tolerance = 1e-16
    end
    Nmax = 13
    T = a isa Real && a > zero(a) ? typeof(float(a)) :
        typeof(complex(float(a)))
    aT = convert(T, a)
    sT = one(T)
    N = _loggamma_divdiff_degree(aT, epsilon, tolerance, Nmax)
    coeffs = loggamma(affine_jet(aT, sT, N)).coeffs
    plus, plus_ok = _loggamma_divdiff_sum_coeffs(coeffs, epsilon, N, tolerance, 1)
    minus, minus_ok = _loggamma_divdiff_sum_coeffs(coeffs, epsilon, N, tolerance, -1)
    if N < Nmax && !(plus_ok && minus_ok)
        coeffs = loggamma(affine_jet(aT, sT, Nmax)).coeffs
        plus_ok || ((plus, _) = _loggamma_divdiff_sum_coeffs(coeffs, epsilon, Nmax, tolerance, 1))
        minus_ok || ((minus, _) = _loggamma_divdiff_sum_coeffs(coeffs, epsilon, Nmax, tolerance, -1))
    end
    return plus, minus
end

"""Evaluate the signed log-Pochhammer divided difference with `log1p`."""
function logpochhammer_divdiff(a, n::Integer, epsilon, sigma::Integer)
    n == 0 && return zero(a + epsilon + 0im)
    total = zero(a + epsilon + 0im)
    for k in 0:(n - 1)
        total += log1p(sigma * epsilon / (a + k)) / epsilon
    end
    return total
end

end # module ConfluentJets
