# Stable sums of log-gamma and digamma values at a shared shifted argument. Numerator and
# denominator factors of `prod_i Γ(t + K + α_i) / prod_j Γ(t + K + β_j)` are paired and each
# pair evaluated as a ratio with the exact difference `α_i - β_j`, so the error scales with the
# result instead of `u * sum |loggamma|`. The `|p - q|` unpaired factors are evaluated directly.

# Shift threshold and Stirling coefficients B_{2k}/(2k(2k-1)) (log-gamma) and
# B_{2k}/(2k) (digamma), k = 1..8: truncation error below 1e-19 for real(w) >= 15.
const _GAMMARATIO_SHIFT = 15
const _GAMMARATIO_LOGGAMMA_COEFFS = (1/12, -1/360, 1/1260, -1/1680, 1/1188,
                                     -691/360360, 1/156, -3617/122400)
const _GAMMARATIO_DIGAMMA_COEFFS = (1/12, -1/120, 1/252, -1/240, 1/132,
                                    -691/32760, 1/12, -3617/8160)
const _GAMMARATIO_GUARD_BITS = 64

"""
    _loggamma_ratio(x, d)

Return `loggamma(x + d) - loggamma(x)` modulo `2πi` for `Complex{Float64}`
arguments, where `d` is supplied exactly (not as a difference of rounded values).
"""
function _loggamma_ratio(x::ComplexF64, d::ComplexF64)
    # log Γ(w+d) - log Γ(w) = [same](w+1) - log((w+d)/w). Collect the shift as one product.
    shift = one(x)
    w = x
    while real(w) < _GAMMARATIO_SHIFT
        shift *= 1 + d / w
        w += 1
    end
    x1 = w + d
    # (x1 - 1/2) log x1 - (w - 1/2) log w - d, rearranged without cancellation.
    value = d * log(w) + (x1 - 0.5) * log1p(d / w) - d
    p1 = inv(x1)
    p2 = inv(w)
    q1 = p1 * p1
    q2 = p2 * p2
    for c in _GAMMARATIO_LOGGAMMA_COEFFS
        value += c * (p1 - p2)
        p1 *= q1
        p2 *= q2
    end
    return value - log(shift)
end

"""
    _digamma_diff(x, d)

Return `digamma(x + d) - digamma(x)` for `Complex{Float64}` arguments, with `d`
supplied exactly.
"""
function _digamma_diff(x::ComplexF64, d::ComplexF64)
    value = zero(x)
    w = x
    while real(w) < _GAMMARATIO_SHIFT
        value += d / (w * (w + d))     # ψ(w+d) - ψ(w) = [same](w+1) + d/(w(w+d))
        w += 1
    end
    x1 = w + d
    value += log1p(d / w) - (inv(x1) - inv(w)) / 2
    q1 = inv(x1 * x1)
    q2 = inv(w * w)
    p1 = q1
    p2 = q2
    for c in _GAMMARATIO_DIGAMMA_COEFFS
        value -= c * (p1 - p2)
        p1 *= q1
        p2 *= q2
    end
    return value
end

"""
    _gammaratio_pairs(plus, minus)

Pair each entry of the shorter offset list with the nearest unused entry of the
longer one.  Returns `(pairs, rest_plus, rest_minus)` with index pairs.
"""
function _gammaratio_pairs(plus, minus)
    np = length(plus)
    nm = length(minus)
    usedp = falses(np)
    usedm = falses(nm)
    pairs = Tuple{Int,Int}[]
    for _ in 1:min(np, nm)
        best = (0, 0)
        bestd = Inf
        for i in 1:np, j in 1:nm
            (usedp[i] || usedm[j]) && continue
            dist = abs(plus[i] - minus[j])
            dist < bestd && (bestd = dist; best = (i, j))
        end
        usedp[best[1]] = true
        usedm[best[2]] = true
        push!(pairs, best)
    end
    return pairs, findall(!, usedp), findall(!, usedm)
end

"""
    _loggamma_sum_shifted(K, plus, minus)

Return `sum loggamma(plus_i + K) - sum loggamma(minus_j + K)` modulo `2πi`.
Offsets are small (formed from parameters). `K >= 0` is an integer shift.
"""
_gammaratio_real(plus, minus) = typeof(float(real(first(promote(plus..., minus...)))))

function _loggamma_sum_shifted(K::Integer, plus, minus)
    R = _gammaratio_real(plus, minus)
    R === BigFloat && return _loggamma_sum_shifted(Complex{BigFloat}, K, plus, minus)
    return convert(Complex{R}, _loggamma_sum_shifted(ComplexF64, K, plus, minus))
end

function _loggamma_sum_shifted(::Type{ComplexF64}, K::Integer, plus, minus)
    pl = ComplexF64[x for x in plus]
    mi = ComplexF64[x for x in minus]
    pairs, rp, rm = _gammaratio_pairs(pl, mi)
    out = zero(ComplexF64)
    for (i, j) in pairs
        out += _loggamma_ratio(mi[j] + K, pl[i] - mi[j])
    end
    for i in rp
        out += loggamma(pl[i] + K)
    end
    for j in rm
        out -= loggamma(mi[j] + K)
    end
    return out
end

# BigFloat: evaluate the direct sum with guard bits on exactly carried inputs.
function _loggamma_sum_shifted(::Type{Complex{BigFloat}}, K::Integer, plus, minus)
    return _gammaratio_guarded(K, plus, minus) do x
        loggamma(x)
    end
end

"""
    _digamma_sum_shifted(K, plus, minus)

Return `sum digamma(plus_i + K) - sum digamma(minus_j + K)`.
"""
function _digamma_sum_shifted(K::Integer, plus, minus)
    R = _gammaratio_real(plus, minus)
    R === BigFloat && return _digamma_sum_shifted(Complex{BigFloat}, K, plus, minus)
    return convert(Complex{R}, _digamma_sum_shifted(ComplexF64, K, plus, minus))
end

function _digamma_sum_shifted(::Type{ComplexF64}, K::Integer, plus, minus)
    pl = ComplexF64[x for x in plus]
    mi = ComplexF64[x for x in minus]
    pairs, rp, rm = _gammaratio_pairs(pl, mi)
    out = zero(ComplexF64)
    for (i, j) in pairs
        out += _digamma_diff(mi[j] + K, pl[i] - mi[j])
    end
    for i in rp
        out += digamma(pl[i] + K)
    end
    for j in rm
        out -= digamma(mi[j] + K)
    end
    return out
end

function _digamma_sum_shifted(::Type{Complex{BigFloat}}, K::Integer, plus, minus)
    return _gammaratio_guarded(K, plus, minus) do x
        digamma(x)
    end
end

function _gammaratio_guarded(f, K::Integer, plus, minus)
    prec = precision(BigFloat)
    hi = setprecision(BigFloat, prec + _GAMMARATIO_GUARD_BITS) do
        acc = zero(Complex{BigFloat})
        for x in plus
            acc += f(Complex{BigFloat}(x) + K)
        end
        for x in minus
            acc -= f(Complex{BigFloat}(x) + K)
        end
        acc
    end
    return complex(BigFloat(real(hi)), BigFloat(imag(hi)))
end
