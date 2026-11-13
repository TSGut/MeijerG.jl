# Minimal double-double arithmetic (`hi + lo`, about 106 bits) for the extended-precision lanes.
# `_DD <: Real` with complex values as `Complex{_DD}` (`_CDD`), so Base and the generic jets apply.
# Only the primitives the lanes need, with Dekker/Knuth error-free transformations and `fma`.
# `abs(::_CDD)` returns a Float64 and there is no atan, so complex powers are `exp(s * log(z))`.

struct _DD <: Real
    hi::Float64
    lo::Float64
end

_DD(x::Float64) = _DD(x, 0.0)
_DD(x::Integer) = _DD(Float64(x))
function _DD(x::BigFloat)
    hi = Float64(x)
    return _DD(hi, Float64(x - hi))
end
_DD(x::BigInt) = _DD(BigFloat(x))
_DD(x::Rational) = _DD(BigFloat(x))
_DD(x::AbstractIrrational) = _DD(BigFloat(x))
_DD(x::Real) = _DD(Float64(x))
_DD(x::_DD) = x

const _CDD = Complex{_DD}

Base.Float64(x::_DD) = x.hi + x.lo
Base.Float32(x::_DD) = Float32(x.hi + x.lo)
Base.BigFloat(x::_DD) = BigFloat(x.hi) + BigFloat(x.lo)
Base.float(::Type{_DD}) = _DD
Base.float(x::_DD) = x
Base.promote_rule(::Type{_DD}, ::Type{<:Union{Float16,Float32,Float64,Integer,Rational}}) = _DD
Base.precision(::Type{_DD}; base::Integer=2) = base == 2 ? 106 : floor(Int, 106 * log(2) / log(base))
Base.eps(::Type{_DD}) = _DD(2.0^-104)
Base.zero(::Type{_DD}) = _DD(0.0)
Base.one(::Type{_DD}) = _DD(1.0)
Base.iszero(a::_DD) = a.hi == 0.0
Base.isfinite(a::_DD) = isfinite(a.hi)
Base.isnan(a::_DD) = isnan(a.hi)
Base.isinf(a::_DD) = isinf(a.hi)
Base.signbit(a::_DD) = signbit(a.hi)
Base.:(==)(a::_DD, b::_DD) = a.hi == b.hi && a.lo == b.lo
Base.isless(a::_DD, b::_DD) = a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo)
Base.:<(a::_DD, b::_DD) = isless(a, b)
Base.:<=(a::_DD, b::_DD) = !isless(b, a)
Base.abs(a::_DD) = signbit(a.hi) ? -a : a
Base.round(a::_DD, r::RoundingMode) = _DD(round(a.hi + a.lo, r))
Base.round(::Type{T}, a::_DD) where {T<:Integer} = round(T, a.hi + a.lo)
Base.round(::Type{T}, a::_DD, r::RoundingMode) where {T<:Integer} = round(T, a.hi + a.lo, r)
# Julia 1.10 does not route ceil/floor/trunc(::Type{T}, ::Real) through `round`.
Base.trunc(::Type{T}, a::_DD) where {T<:Integer} = trunc(T, a.hi + a.lo)
Base.ceil(::Type{T}, a::_DD) where {T<:Integer} = ceil(T, a.hi + a.lo)
Base.floor(::Type{T}, a::_DD) where {T<:Integer} = floor(T, a.hi + a.lo)
Base.ldexp(a::_DD, k::Integer) = _DD(ldexp(a.hi, k), ldexp(a.lo, k))
Base.show(io::IO, a::_DD) = print(io, "_DD(", a.hi, ", ", a.lo, ")")

@inline function _dd_two_sum(a::Float64, b::Float64)
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)
end

@inline function _dd_fast_two_sum(a::Float64, b::Float64)
    s = a + b
    return s, b - (s - a)
end

@inline _dd_two_prod(a::Float64, b::Float64) = (p = a * b; (p, fma(a, b, -p)))

@inline function Base.:+(a::_DD, b::_DD)
    s, e = _dd_two_sum(a.hi, b.hi)
    t, f = _dd_two_sum(a.lo, b.lo)
    e += t
    s, e = _dd_fast_two_sum(s, e)
    e += f
    return _DD(_dd_fast_two_sum(s, e)...)
end

@inline function Base.:+(a::_DD, b::Float64)
    s, e = _dd_two_sum(a.hi, b)
    e += a.lo
    return _DD(_dd_fast_two_sum(s, e)...)
end

Base.:+(a::Float64, b::_DD) = b + a
Base.:-(a::_DD) = _DD(-a.hi, -a.lo)
Base.:-(a::_DD, b::_DD) = a + (-b)
Base.:-(a::_DD, b::Float64) = a + (-b)
Base.:-(a::Float64, b::_DD) = (-b) + a

@inline function Base.:*(a::_DD, b::_DD)
    p, e = _dd_two_prod(a.hi, b.hi)
    e += a.hi * b.lo + a.lo * b.hi
    return _DD(_dd_fast_two_sum(p, e)...)
end

@inline function Base.:*(a::_DD, b::Float64)
    p, e = _dd_two_prod(a.hi, b)
    e += a.lo * b
    return _DD(_dd_fast_two_sum(p, e)...)
end

Base.:*(a::Float64, b::_DD) = b * a

function Base.:/(a::_DD, b::_DD)
    q1 = a.hi / b.hi
    r = a - b * q1
    q2 = r.hi / b.hi
    r = r - b * q2
    q3 = r.hi / b.hi
    return _DD(_dd_fast_two_sum(q1, q2)...) + q3
end

Base.:/(a::_DD, b::Float64) = a / _DD(b)
Base.inv(a::_DD) = _DD(1.0) / a
Base.:^(a::_DD, n::Integer) = n >= 0 ? Base.power_by_squaring(a, n) : Base.power_by_squaring(inv(a), -n)
Base.:^(z::_CDD, n::Integer) = n >= 0 ? Base.power_by_squaring(z, n) : Base.power_by_squaring(inv(z), -n)
_dd_ldexp(a::_DD, k::Integer) = ldexp(a, k)

function Base.sqrt(a::_DD)                        # a >= 0
    iszero(a) && return a
    y = _DD(sqrt(a.hi))
    return y + (a - y * y) / (2.0 * y.hi)
end

# Constants rounded from 256-bit values.
const _DD_LN2 = _DD(log(big(2)))
const _DD_PI = _DD(big(pi))
const _DD_HALF_PI = _DD(big(pi) / 2)
const _DD_HALF_LOG_2PI = _DD(log(2 * big(pi)) / 2)
const _DD_INV_FACT = Tuple(_DD(1 / factorial(big(k))) for k in 0:30)

function _dd_exp(a::_DD)
    isnan(a.hi) && return _DD(NaN)
    a.hi > 709.0 && return _DD(Inf)
    a.hi < -745.0 && return _DD(0.0)
    k = round(a.hi / _DD_LN2.hi)
    r = ldexp(a - _DD_LN2 * k, -4)                # |r| <= 0.022
    s = _DD_INV_FACT[16]
    for i in 15:-1:1
        s = s * r + _DD_INV_FACT[i]
    end
    for _ in 1:4
        s = s * s
    end
    return ldexp(s, Int(k))
end

function _dd_log(a::_DD)                          # a > 0
    (isfinite(a.hi) && a.hi > 0.0) || return _DD(a.hi > 0.0 ? Inf : a.hi == 0.0 ? -Inf : NaN)
    y = _DD(log(a.hi))
    return y + (a * _dd_exp(-y) - 1.0)            # one Newton step
end

# sin and cos on |r| <= π/4 by Taylor series (terms below 1e-33).
function _dd_sincos_reduced(r::_DD)
    r2 = r * r
    s = _DD_INV_FACT[30]                          # 1/29!
    c = _DD_INV_FACT[29]                          # 1/28!
    for k in 13:-1:0
        s = s * r2 * (-1.0) + _DD_INV_FACT[2k + 2]
        c = c * r2 * (-1.0) + _DD_INV_FACT[2k + 1]
    end
    return s * r, c
end

function _dd_sincos(a::_DD)
    isfinite(a.hi) || return _DD(NaN), _DD(NaN)
    k = round(a.hi / _DD_HALF_PI.hi)
    s, c = _dd_sincos_reduced(a - _DD_HALF_PI * k)
    q = mod(Int(k), 4)
    q == 0 && return s, c
    q == 1 && return c, -s
    q == 2 && return -s, -c
    return -c, s
end

Base.exp(a::_DD) = _dd_exp(a)
Base.log(a::_DD) = _dd_log(a)
Base.sincos(a::_DD) = _dd_sincos(a)
Base.sin(a::_DD) = _dd_sincos(a)[1]
Base.cos(a::_DD) = _dd_sincos(a)[2]

# ---------------------------------------------------------------------------
# Complex double-double.
# ---------------------------------------------------------------------------

Base.abs(z::_CDD) = abs(ComplexF64(z))

function _dd_exp(z::_CDD)
    e = _dd_exp(real(z))
    s, c = _dd_sincos(imag(z))
    return Complex(e * c, e * s)
end

function _dd_log(z::_CDD)
    y = _CDD(log(ComplexF64(z)))
    w = z * _dd_exp(-y) - 1.0                     # |w| ~ 1e-16
    return y + w - w * w * 0.5
end

Base.exp(z::_CDD) = _dd_exp(z)
Base.log(z::_CDD) = _dd_log(z)

# log1p and expm1 by Taylor series for small arguments (the lanes use them
# with |x| far below 1), otherwise through log and exp.
function Base.log1p(z::_CDD)
    abs(z) < 0.25 || return _dd_log(z + 1.0)
    term = z
    total = z
    for k in 2:60
        term = term * z * (-1.0)
        v = term / k
        total = total + v
        abs(v) <= 1e-34 * abs(total) && break
    end
    return total
end
Base.log1p(a::_DD) = real(log1p(_CDD(a)))

function Base.expm1(z::_CDD)
    abs(z) < 0.5 || return _dd_exp(z) - 1.0
    term = z
    total = z
    for k in 2:40
        term = term * z / k
        total = total + term
        abs(term) <= 1e-34 * abs(total) && break
    end
    return total
end
Base.expm1(a::_DD) = real(expm1(_CDD(a)))

function _dd_sinpi(z::_CDD)                       # sin(πz)
    x = _DD_PI * real(z)
    y = _DD_PI * imag(z)
    s, c = _dd_sincos(x)
    ep = _dd_exp(y)
    em = inv(ep)
    return Complex(s * ((ep + em) * 0.5), c * ((ep - em) * 0.5))
end

function _dd_cospi(z::_CDD)                       # cos(πz)
    x = _DD_PI * real(z)
    y = _DD_PI * imag(z)
    s, c = _dd_sincos(x)
    ep = _dd_exp(y)
    em = inv(ep)
    return Complex(c * ((ep + em) * 0.5), -(s * ((ep - em) * 0.5)))
end

Base.sinpi(z::_CDD) = _dd_sinpi(z)
Base.cospi(z::_CDD) = _dd_cospi(z)
Base.sinpi(a::_DD) = real(_dd_sinpi(_CDD(a)))
Base.cospi(a::_DD) = real(_dd_cospi(_CDD(a)))

# log-gamma (mod 2πi) and digamma by upward shift to real part >= 24 and 16 Stirling terms.

function _dd_bernoulli_even(n)
    A = Vector{Rational{BigInt}}(undef, 2n + 1)
    B = Rational{BigInt}[]
    for mm in 0:2n
        A[mm + 1] = 1 // (mm + 1)
        for j in mm:-1:1
            A[j] = j * (A[j] - A[j + 1])
        end
        mm >= 2 && iseven(mm) && push!(B, A[1])
    end
    return B
end

const _DD_B2K = _dd_bernoulli_even(16)
const _DD_LOGGAMMA_COEFFS = Tuple(_DD(BigFloat(_DD_B2K[k] / (2k * (2k - 1)))) for k in 1:16)
const _DD_DIGAMMA_COEFFS = Tuple(_DD(BigFloat(_DD_B2K[k] / (2k))) for k in 1:16)
const _DD_SHIFT = 24.0

function _dd_loggamma(x::_CDD)
    logshift = zero(_CDD)
    prod = one(_CDD)
    w = x
    steps = 0
    while real(w).hi < _DD_SHIFT
        prod = prod * w
        w = w + 1
        steps += 1
        if steps % 32 == 0                      # keep the product in range
            logshift = logshift + _dd_log(prod)
            prod = one(_CDD)
        end
    end
    logshift = logshift + _dd_log(prod)
    lw = _dd_log(w)
    val = (w - 0.5) * lw - w + _DD_HALF_LOG_2PI
    p = inv(w)
    q = p * p
    for c in _DD_LOGGAMMA_COEFFS
        val = val + p * c
        p = p * q
    end
    return val - logshift
end

function _dd_digamma(x::_CDD)
    val = zero(_CDD)
    w = x
    while real(w).hi < _DD_SHIFT
        val = val - inv(w)
        w = w + 1
    end
    val = val + _dd_log(w) - inv(w) * 0.5
    q = inv(w * w)
    p = q
    for c in _DD_DIGAMMA_COEFFS
        val = val - p * c
        p = p * q
    end
    return val
end

ConfluentJets._jet_fast_generic(::Type{_DD}) = true

Gamma.loggamma(z::_CDD) = _dd_loggamma(z)
Gamma.digamma(z::_CDD) = _dd_digamma(z)
Gamma.gamma(z::_CDD) = _dd_exp(_dd_loggamma(z))
function Gamma.loggamma(a::_DD)
    signbit(a.hi) && throw(DomainError(a, "loggamma(::_DD) requires a positive argument"))
    return real(_dd_loggamma(_CDD(a)))
end
Gamma.digamma(a::_DD) = real(_dd_digamma(_CDD(a)))
Gamma.gamma(a::_DD) = real(_dd_exp(_dd_loggamma(_CDD(a))))
function Gamma.logabsgamma(a::_DD)
    lg = _dd_loggamma(_CDD(a))
    return real(lg), cos(Float64(imag(lg))) >= 0 ? one(_DD) : -one(_DD)
end
