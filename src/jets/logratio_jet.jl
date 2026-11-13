# Private exact and near-confluent log-ratio jet support for guarded lower and upper active-pole
# pairs, returning structured results so `meijerg` can fall back without changing raw Slater.

const _LOGRATIO_JET_MAX_ABS_EPSILON = 1.1e-3
const _LOGRATIO_JET_MAX_ABS_EPSILON_LOGZ = 2e-1
const _LOGRATIO_JET_MIN_LOWER_RAW_CANCELLATION = 10.0
const _LOGRATIO_JET_MAX_GROUPED_TERMS = 50
const _LOGRATIO_JET_PAIR_ATOL = 1.1e-2
const _LOGRATIO_JET_PAIR_RTOL = 1e-9
const _LOGRATIO_JET_TOL_FACTOR = 10
const _LOGRATIO_JET_REFRESH = 64
const _LOGRATIO_JET_EXACT_SNAP_ULPS = 4

# Series-level conditioning guard: kappa = sum_n |P_n| / |sum_n P_n| certifies the digits lost to
# cancellation between terms. It grows like exp(nu*|z|^(1/nu)) off the convergent regime.
const _LOGRATIO_JET_MAX_SUMMATION_KAPPA = 1e2
const _LOGRATIO_JET_ESCALATION_GUARD_BITS = 32
const _LOGRATIO_JET_MAX_ESCALATION_FACTOR = 4

_logratio_summation_kappa(absum, value) = absum / abs(value)
_logratio_summation_ok(kappa) = kappa <= _LOGRATIO_JET_MAX_SUMMATION_KAPPA

_logratio_reprecision(x::BigFloat) = BigFloat(x)
_logratio_reprecision(x::Complex{BigFloat}) = complex(BigFloat(real(x)), BigFloat(imag(x)))
_logratio_reprecision(x::Tuple) = map(_logratio_reprecision, x)
_logratio_reprecision(x) = x

"""
    _logratio_escalated(f, args, z, kappa; min_extra=0)

Rerun `f(args...)` at a working precision raised by `log2(kappa)` plus guard bits, and
by at least `min_extra` bits, with the (exact) inputs carried to that precision, and
round the result's `value` back.  Returns `nothing` when `z` is not `BigFloat`-based or
the required extra precision exceeds `_LOGRATIO_JET_MAX_ESCALATION_FACTOR` times the
current one.
"""
function _logratio_escalated(f, args::Tuple, z, kappa; min_extra::Integer=0)
    real(z) isa BigFloat || return nothing
    isfinite(kappa) || return nothing
    prec = precision(BigFloat)
    extra = max(min_extra, ceil(Int, log2(kappa)) + _LOGRATIO_JET_ESCALATION_GUARD_BITS)
    extra <= _LOGRATIO_JET_MAX_ESCALATION_FACTOR * prec || return nothing
    hi = setprecision(BigFloat, prec + extra) do
        f(map(_logratio_reprecision, args)...)
    end
    return merge(hi, (value=_logratio_reprecision(hi.value),))
end

_logratio_jet_enabled() = get(ENV, "MEIJERG_LOGRATIO_JET", "1") != "0"
_logratio_nan_value(z) = convert(Complex{Float64}, NaN + NaN * im)
function _logratio_exact_tol(z)
    return _logratio_series_tol(z)
end
function _logratio_exact_snap_tol(x, y)
    R = _logratio_real_float_type(x + y)
    scale = max(abs(x), abs(y), one(abs(x)))
    return convert(R, _LOGRATIO_JET_EXACT_SNAP_ULPS) * eps(R) * scale
end
function _logratio_series_tol(z)
    R = _logratio_real_float_type(z)
    return convert(R, _LOGRATIO_JET_TOL_FACTOR) * eps(R)
end

"""
    _logratio_refresh_every(branch, z)

Return the recurrence refresh interval for a selected near-confluent log-ratio
jet series.  The usual lower-small-`z` and upper-large-`z` regimes use the fast
recurrence, while the opposite magnitude regimes refresh from the log-gamma
product every term to suppress recurrence drift in cancellation-sensitive cases.
(The exact-limit series refreshes in ratio form and always uses
`_LOGRATIO_JET_REFRESH`.)
"""
function _logratio_refresh_every(branch::Symbol, z)
    az = abs(z)
    oneaz = one(az)
    off_regime = (branch === :lower && az > oneaz) ||
        (branch === :upper && az < oneaz)
    return off_regime ? 1 : _LOGRATIO_JET_REFRESH
end

# Compensated summation (Sum2 of Ogita, Rump & Oishi, 2005): `total` plus the accumulated
# rounding errors `corr`, with error at most u|S| + γ²_{n-1} Σ|T_n|.
@inline function _logratio_sum2(total, corr, v)
    t = total + v
    bb = t - total
    return t, corr + ((total - (t - bb)) + (v - bb))
end

function _logratio_result(status::Symbol, value, branch::Symbol, selected::Bool,
                          reason::Symbol, diagnostics=(;))
    return (status=status, value=value, branch=branch, selected=selected,
            reason=reason, diagnostics=diagnostics)
end

# Cancellation-free tail majorant for the grouped pair series. The prefactor ratio from index
# r-1 to r is rho * prod |r + c_num| / prod |r + c_den| with |r+c| <= r(1+|c|/R) and
# |r+c| >= r(1+min(real(c),0)/R) for r >= R, and D advances by log1p(epsilon/(r+c))/epsilon, so
# |Delta D| <= (1/r) / (1 + (min(real(c),0)-|epsilon|)/R) and |D_r| <= M_r grows logarithmically.
# With h(M) = M*exprel(|epsilon|M) and d(log h)/dM <= |epsilon|+1 for M >= 1 the envelope has a
# uniform geometric ratio once theta < 1, without a 1/|epsilon| separated-stream bound.
function _logratio_tail_constants(a, b, L::Integer, epsv, branch::Symbol)
    # Concretely typed vectors, fixed once from the argument types, so that the per-term tail
    # bound does not box every element.
    T = promote_type(typeof(epsv), eltype(a), eltype(b), Float64)
    if branch === :lower
        anchor = b[1] + L
        numerator = T[anchor - aj for aj in a]
        denominator = T[T(L), -epsv]
        append!(denominator, (anchor - b[j] for j in 3:length(b)))
        increments = T[anchor - aj for aj in a]
        append!(increments, (anchor - b[j] for j in 3:length(b)))
        push!(increments, zero(T), T(L))
    elseif branch === :upper
        anchor = L - a[1]
        numerator = T[anchor + bj for bj in b]
        denominator = T[T(L), epsv]
        append!(denominator, (anchor + a[j] for j in 3:length(a)))
        increments = T[anchor + bj for bj in b]
        append!(increments, (anchor + a[j] for j in 3:length(a)))
        push!(increments, T(L), zero(T))
    else
        throw(ArgumentError("tail constants require a lower or upper branch"))
    end
    return numerator, denominator, increments
end

"""
    _logratio_pair_tail_bound(P, D, a, b, z, L, epsilon, n, branch)

Bound the uncomputed grouped-pair tail after index n.  The returned
tail bound and theta are finite only after the rational recurrence
majorants are contractive.  The bound is mathematical in exact arithmetic.
The caller inflates it modestly before using it as a floating-point stopping
test.
"""
function _logratio_pair_tail_bound(P, D, a, b, z, L::Integer, epsv,
                                   nn::Integer, branch::Symbol;
                                   constants=nothing)
    R = nn + 1
    numerator, denominator, increments = constants === nothing ?
        _logratio_tail_constants(a, b, L, epsv, branch) : constants
    RT = typeof(float(abs(P) + abs(D) + abs(epsv) + abs(z)))
    r = convert(RT, R)
    epsabs = convert(RT, abs(epsv))
    rho = branch === :lower ? convert(RT, abs(z)) : convert(RT, abs(inv(z)))
    d = branch === :lower ? length(b) - length(a) : length(a) - length(b)
    d > 0 || return (tail_bound=convert(RT, Inf), theta=convert(RT, Inf))

    numscale = one(RT)
    for c in numerator
        numscale *= one(RT) + convert(RT, abs(c)) / r
    end
    denscale = one(RT)
    for c in denominator
        f = one(RT) + min(convert(RT, real(c)), zero(RT)) / r
        f > zero(RT) || return (tail_bound=convert(RT, Inf), theta=convert(RT, Inf))
        denscale *= f
    end
    qbar = rho * numscale / (denscale * r^d)

    cD = zero(RT)
    for c in increments
        f = one(RT) +
            (min(convert(RT, real(c)), zero(RT)) - epsabs) / r
        f > zero(RT) || return (tail_bound=convert(RT, Inf), theta=convert(RT, Inf))
        cD += inv(f)
    end
    deltaM = cD / r
    M = one(RT) + convert(RT, abs(D))
    theta = qbar * exp((epsabs + one(RT)) * deltaM)
    theta < one(RT) || return (tail_bound=convert(RT, Inf), theta=theta)

    Enext = convert(RT, abs(P)) * qbar *
        _logratio_tail_h(epsabs, M + deltaM)
    # Allow for ordinary roundoff in this positive floating-point majorant.
    safety = one(RT) + convert(RT, 64) * eps(RT)
    theta_safe = safety * theta
    theta_safe < one(RT) ||
        return (tail_bound=convert(RT, Inf), theta=theta_safe)
    return (tail_bound=safety * Enext / (one(RT) - theta_safe),
            theta=theta_safe)
end

"""
    _detect_exact_logratio_jet_pair(a, b, m, n, branch)

Detect and canonicalize the single active exact pair supported by the
log-ratio exact-limit series.  The formulas use canonical pair positions
`(1, 2)` with lower convention `b[2] = b[1] + L` and upper convention
`a[2] = a[1] - L`. Active parameters are therefore permuted internally when
the caller supplies the colliding families in the opposite order.
"""
function _detect_exact_logratio_jet_pair(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                         branch::Symbol)
    active = branch === :lower ? b[1:m] : branch === :upper ? a[1:n] : ()
    length(active) <= 1 && return (status=:not_near_confluent, pair=nothing, a=a, b=b)

    candidates = Any[]
    has_near = false
    for i in eachindex(active), j in (i + 1):length(active)
        δ = active[j] - active[i]
        signed_L = round(Int, real(δ))
        residual = δ - signed_L
        if abs(residual) <= _logratio_exact_snap_tol(δ, signed_L)
            L = abs(signed_L)
            if branch === :lower
                first, second = signed_L >= 0 ? (i, j) : (j, i)
            else
                first, second = signed_L <= 0 ? (i, j) : (j, i)
            end
            push!(candidates, (i=first, j=second, L=L, epsilon=zero(δ),
                               abs_epsilon=zero(abs(δ)), branch=branch))
        elseif abs(residual) <= _confluent_near_tol(δ, signed_L)
            has_near = true
        end
    end

    has_near && isempty(candidates) && return (status=:not_near_confluent, pair=nothing, a=a, b=b)
    isempty(candidates) && return (status=:not_near_confluent, pair=nothing, a=a, b=b)

    involved = Set{Int}()
    for c in candidates
        push!(involved, c.i); push!(involved, c.j)
    end
    length(involved) > 2 && return (status=:unsupported_higher_multiplicity, pair=nothing, a=a, b=b)
    length(candidates) > 1 && return (status=:unsupported_multiple_pairs, pair=nothing, a=a, b=b)

    pair = candidates[1]
    if branch === :lower
        active_order = [pair.i, pair.j, (k for k in 1:m if k != pair.i && k != pair.j)...]
        bcanon = ((b[k] for k in active_order)..., b[(m + 1):end]...)
        return (status=:ok, pair=merge(pair, (i=1, j=2)), a=a, b=bcanon)
    else
        active_order = [pair.i, pair.j, (k for k in 1:n if k != pair.i && k != pair.j)...]
        acanon = ((a[k] for k in active_order)..., a[(n + 1):end]...)
        return (status=:ok, pair=merge(pair, (i=1, j=2)), a=acanon, b=b)
    end
end

"""
    _detect_near_logratio_jet_pair(a, b, m, n, branch; atol, rtol)

Detect and canonicalize the single active near pair supported by the private
log-ratio jet candidate.  Lower pairs are tested with `b[j] - b[i] ≈ L`,
upper pairs with `a[j] - a[i] ≈ -L`.  Nonzero separations are kept in the
finite-epsilon log-ratio jet even when they are below the exact-confluence
tolerance. Only a literally zero separation is left for the exact-zero
log-ratio path.

As for the exact detector, the pair is moved to the canonical positions
`(1, 2)` with lower convention `b[2] = b[1] + L + epsilon` and upper convention
`a[2] = a[1] - L + epsilon`, `L >= 0`. The colliding families are swapped when
the caller supplies them in the opposite order.  Returns `(status, pair, a, b)`
with the canonical parameter tuples.

The scan tracks the single admissible candidate directly: two distinct
candidate pairs always involve at least three active indices, so more than
one candidate is already a higher-multiplicity rejection.
"""
function _detect_near_logratio_jet_pair(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                        branch::Symbol; atol=_LOGRATIO_JET_PAIR_ATOL,
                                        rtol=_LOGRATIO_JET_PAIR_RTOL)
    active = branch === :lower ? b[1:m] : branch === :upper ? a[1:n] : ()
    isempty(active) && return (status=:unsupported_branch, pair=nothing, a=a, b=b)
    length(active) == 1 && return (status=:not_near_confluent, pair=nothing, a=a, b=b)
    found = 0
    first = 0
    second = 0
    bestL = 0
    besteps = zero(active[2] - active[1])
    for i in eachindex(active), j in (i + 1):length(active)
        δ = active[j] - active[i]
        signed_L = round(Int, real(δ))
        epsv = δ - signed_L
        abs(epsv) <= atol + rtol * max(abs(δ), one(abs(δ))) || continue
        abs(epsv) <= _logratio_exact_snap_tol(δ, signed_L) && continue
        found += 1
        found == 1 || continue
        if branch === :lower
            first, second, bestL, besteps = signed_L >= 0 ? (i, j, signed_L, epsv) :
                                                            (j, i, -signed_L, -epsv)
        else
            first, second, bestL, besteps = signed_L <= 0 ? (i, j, -signed_L, epsv) :
                                                            (j, i, signed_L, -epsv)
        end
    end
    found == 0 && return (status=:not_near_confluent, pair=nothing, a=a, b=b)
    found > 1 && return (status=:unsupported_higher_multiplicity, pair=nothing, a=a, b=b)
    pair = (i=1, j=2, L=bestL, epsilon=besteps, abs_epsilon=abs(besteps), branch=branch)
    if branch === :lower
        active_order = [first, second, (k for k in 1:m if k != first && k != second)...]
        bcanon = ((b[k] for k in active_order)..., b[(m + 1):end]...)
        return (status=:ok, pair=pair, a=a, b=bcanon)
    end
    active_order = [first, second, (k for k in 1:n if k != first && k != second)...]
    acanon = ((a[k] for k in active_order)..., a[(n + 1):end]...)
    return (status=:ok, pair=pair, a=acanon, b=b)
end

"""
    _logratio_validated_z(branch, z)

Return whether `z` is admissible for the near-confluent pair jet.  Only a thin
neighborhood of the `log(z)` branch cut is excluded (`|arg z| > 0.95*pi`, where
the value is two-valued). A 40-family sweep showed the accuracy is flat or even
improves toward the cut.  The magnitude of `z` is not restricted: off the
convergent regime of a branch the summation-conditioning certificate
(`_logratio_summation_kappa`) and the double-double lane take over, as for the
exact-limit series.
"""
function _logratio_validated_z(branch::Symbol, z)
    branch in (:lower, :upper) || return false
    return abs(angle(z)) <= 0.95 * pi
end

function _logratio_lower_logH(a, b, m, nactive, s)
    out = zero(s + 0im)
    signH = one(real(s))
    T = typeof(s)
    for j in 3:m
        lv, sv, _ = _loggamma_factor_eval(b[j] - s, T)
        out += lv
        signH *= sv
    end
    for j in 1:nactive
        lv, sv, _ = _loggamma_factor_eval(1 - a[j] + s, T)
        out += lv
        signH *= sv
    end
    for j in (m + 1):length(b)
        lv, sv, _ = _loggamma_factor_eval(1 + s - b[j], T)
        out -= lv
        signH *= sv
    end
    for j in (nactive + 1):length(a)
        lv, sv, _ = _loggamma_factor_eval(a[j] - s, T)
        out -= lv
        signH *= sv
    end
    return out, signH
end

function _logratio_lower_residue_coeff(a, b, m, nactive, k::Integer, ell::Integer)
    bk = b[k]
    s = bk + ell
    coeff = (-one(s))^ell / gamma(_logratio_int_as(s, ell + 1))
    for j in 1:m
        j == k && continue
        coeff *= gamma(b[j] - s)
    end
    for j in 1:nactive
        coeff *= gamma(1 + s - a[j])
    end
    for j in (m + 1):length(b)
        coeff /= gamma(1 + s - b[j])
    end
    for j in (nactive + 1):length(a)
        coeff /= gamma(a[j] - s)
    end
    return coeff
end

function _logratio_lower_single_family_term(a, b, m, nactive, z, k::Integer)
    ap, bp, zp = _promote_inputs(Tuple(a), Tuple(b), z)
    T = eltype(bp)
    p = length(ap)
    q = length(bp)
    bk = bp[k]
    α = ntuple(j -> one(T) + bk - ap[j], p)
    β = ntuple(q - 1) do idx
        j = idx < k ? idx : idx + 1
        one(T) + bk - bp[j]
    end
    C = _logratio_family_prefactor(ap, bp, m, nactive, k, :lower)
    C === nothing && (C = _gamma_prod_ratio(
        Iterators.flatten(((bp[j] - bk for j in 1:m if j != k),
                           (one(T) + bk - ap[j] for j in 1:nactive))),
        Iterators.flatten(((one(T) + bk - bp[j] for j in m+1:q),
                           (ap[j] - bk for j in nactive+1:p))), T))
    return C * zp^bk * pFq(α, β, _signed_argument(zp, p - m - nactive))
end

# Lower convention: the near pair is `b[2] = b[1] + L + epsilon`.
# The grouped coefficient aligns family 1 at `K = L + n` with family 2 at `n`.
function _logratio_lower_D(a, b, m, nactive, z, L::Integer, epsv, nn::Integer)
    K = L + nn
    s = b[1] + K
    one_s = one(s)
    plus_one, minus_one = _logratio_loggamma_divdiff_pair(one_s, epsv)
    D = minus_one - plus_one
    for j in 3:m
        D += _logratio_loggamma_divdiff_signed(b[j] - s, epsv, -1)
    end
    for j in 1:nactive
        D += _logratio_loggamma_divdiff_signed(1 - a[j] + s, epsv, +1)
    end
    for j in (m + 1):length(b)
        D -= _logratio_loggamma_divdiff_signed(1 + s - b[j], epsv, +1)
    end
    for j in (nactive + 1):length(a)
        D -= _logratio_loggamma_divdiff_signed(a[j] - s, epsv, -1)
    end
    D += log(z)
    D += _logratio_logpoch_divdiff_signed(one_s, nn, epsv, -1)
    D -= _logratio_logpoch_divdiff_signed(one_s, K, epsv, +1)
    return D
end

function _logratio_lower_D_limit_direct(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    s = b[1] + K
    one_s = one(s)
    D = _logratio_loggamma_divdiff_limit(one_s, -1) -
        _logratio_loggamma_divdiff_limit(one_s, +1)
    for j in 3:m
        D += _logratio_loggamma_divdiff_limit(b[j] - s, -1)
    end
    for j in 1:nactive
        D += _logratio_loggamma_divdiff_limit(1 - a[j] + s, +1)
    end
    for j in (m + 1):length(b)
        D -= _logratio_loggamma_divdiff_limit(1 + s - b[j], +1)
    end
    for j in (nactive + 1):length(a)
        D -= _logratio_loggamma_divdiff_limit(a[j] - s, -1)
    end
    D += log(z)
    D += _logratio_logpoch_divdiff_limit(one_s, nn, -1)
    D -= _logratio_logpoch_divdiff_limit(one_s, K, +1)
    return D
end

# The grouped pair term is `P * D * exprel(epsilon * D)` with `P = (-1)^(L+1) z^s exp(logX)`.
# The sign is applied here so the recurrence only updates the local product.
function _logratio_lower_pair_state(a, b, m, nactive, z, L::Integer, epsv, nn::Integer)
    K = L + nn
    s = b[1] + K
    logH, signH = _logratio_lower_logH(a, b, m, nactive, s)
    logX = loggamma(1 + epsv) + logH -
           loggamma(_logratio_int_as(s, K + 1)) - _logratio_logpoch(1 - epsv, nn)
    D = _logratio_lower_D(a, b, m, nactive, z, L, epsv, nn)
    return signH * (-1)^(L + 1) * z^s * exp(logX), D
end

function _logratio_lower_pair_term(a, b, m, nactive, z, L::Integer, epsv, nn::Integer)
    P, D = _logratio_lower_pair_state(a, b, m, nactive, z, L, epsv, nn)
    return P * D * _logratio_exprel(epsv * D)
end

function _logratio_lower_pair_state_limit_direct(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    s = b[1] + K
    logH, signH = _logratio_lower_logH(a, b, m, nactive, s)
    logX = logH -
           loggamma(_logratio_int_as(s, K + 1)) - _logratio_logpoch(one(s), nn)
    D = _logratio_lower_D_limit_direct(a, b, m, nactive, z, L, nn)
    return signH * (-1)^(L + 1) * z^s * exp(logX), D
end

function _logratio_lower_pair_term_limit(a, b, m, nactive, z, L::Integer, nn::Integer)
    P, D = _logratio_lower_pair_state_limit(a, b, m, nactive, z, L, nn)
    return P * D
end

function _logratio_upper_logH(a, b, m, nactive, s; skip=(1, 2))
    out = zero(s + 0im)
    signH = one(real(s))
    T = typeof(s)
    for j in 1:m
        lv, sv, _ = _loggamma_factor_eval(b[j] - s, T)
        out += lv
        signH *= sv
    end
    for j in 1:nactive
        j in skip && continue
        lv, sv, _ = _loggamma_factor_eval(1 - a[j] + s, T)
        out += lv
        signH *= sv
    end
    for j in (m + 1):length(b)
        lv, sv, _ = _loggamma_factor_eval(1 - b[j] + s, T)
        out -= lv
        signH *= sv
    end
    for j in (nactive + 1):length(a)
        lv, sv, _ = _loggamma_factor_eval(a[j] - s, T)
        out -= lv
        signH *= sv
    end
    return out, signH
end

function _logratio_upper_unpaired_residue(a, b, m, nactive, z, h::Integer, ell::Integer)
    s = a[h] - 1 - ell
    coeff = (-one(s))^ell / gamma(_logratio_int_as(s, ell + 1))
    for j in 1:m
        coeff *= gamma(b[j] - s)
    end
    for j in 1:nactive
        j == h && continue
        coeff *= gamma(1 - a[j] + s)
    end
    for j in (m + 1):length(b)
        coeff /= gamma(1 - b[j] + s)
    end
    for j in (nactive + 1):length(a)
        coeff /= gamma(a[j] - s)
    end
    return coeff * z^s
end

function _logratio_upper_single_family_term(a, b, m, nactive, z, h::Integer)
    ap, bp, zp = _promote_inputs(Tuple(a), Tuple(b), z)
    T = eltype(ap)
    p = length(ap)
    q = length(bp)
    ah = ap[h]
    α = ntuple(j -> one(T) - ah + bp[j], q)
    β = ntuple(p - 1) do idx
        j = idx < h ? idx : idx + 1
        one(T) - ah + ap[j]
    end
    C = _logratio_family_prefactor(ap, bp, m, nactive, h, :upper)
    C === nothing && (C = _gamma_prod_ratio(
        Iterators.flatten(((ah - ap[j] for j in 1:nactive if j != h),
                           (one(T) - ah + bp[j] for j in 1:m))),
        Iterators.flatten(((one(T) - ah + ap[j] for j in nactive+1:p),
                           (ah - bp[j] for j in m+1:q))), T))
    return C * zp^(ah - one(T)) * pFq(α, β, _signed_argument(inv(zp), q - m - nactive))
end

# Upper convention: the near pair is `a[2] = a[1] - L + epsilon`.
# The grouped coefficient aligns family 1 at `K = L + n` with family 2 at `n`.
function _logratio_upper_D(a, b, m, nactive, z, L::Integer, epsv, nn::Integer)
    K = L + nn
    s = a[1] - 1 - K
    one_s = one(s)
    plus_one, minus_one = _logratio_loggamma_divdiff_pair(one_s, epsv)
    D = plus_one - minus_one
    for j in 1:m
        D += _logratio_loggamma_divdiff_signed(b[j] - s, epsv, -1)
    end
    for j in 1:nactive
        (j == 1 || j == 2) && continue
        D += _logratio_loggamma_divdiff_signed(1 - a[j] + s, epsv, +1)
    end
    for j in (m + 1):length(b)
        D -= _logratio_loggamma_divdiff_signed(1 - b[j] + s, epsv, +1)
    end
    for j in (nactive + 1):length(a)
        D -= _logratio_loggamma_divdiff_signed(a[j] - s, epsv, -1)
    end
    D += log(z)
    D -= _logratio_logpoch_divdiff_signed(one_s, K, epsv, -1)
    D += _logratio_logpoch_divdiff_signed(one_s, nn, epsv, +1)
    return D
end

function _logratio_upper_D_limit_direct(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    s = a[1] - 1 - K
    one_s = one(s)
    D = _logratio_loggamma_divdiff_limit(one_s, +1) -
        _logratio_loggamma_divdiff_limit(one_s, -1)
    for j in 1:m
        D += _logratio_loggamma_divdiff_limit(b[j] - s, -1)
    end
    for j in 1:nactive
        (j == 1 || j == 2) && continue
        D += _logratio_loggamma_divdiff_limit(1 - a[j] + s, +1)
    end
    for j in (m + 1):length(b)
        D -= _logratio_loggamma_divdiff_limit(1 - b[j] + s, +1)
    end
    for j in (nactive + 1):length(a)
        D -= _logratio_loggamma_divdiff_limit(a[j] - s, -1)
    end
    D += log(z)
    D -= _logratio_logpoch_divdiff_limit(one_s, K, -1)
    D += _logratio_logpoch_divdiff_limit(one_s, nn, +1)
    return D
end

function _logratio_upper_pair_state(a, b, m, nactive, z, L::Integer, epsv, nn::Integer)
    K = L + nn
    s = a[1] - 1 - K
    logH, signH = _logratio_upper_logH(a, b, m, nactive, s)
    logX = loggamma(1 - epsv) + logH -
           loggamma(_logratio_int_as(s, K + 1)) - _logratio_logpoch(1 + epsv, nn)
    D = _logratio_upper_D(a, b, m, nactive, z, L, epsv, nn)
    return signH * (-1)^L * z^s * exp(logX), D
end

function _logratio_upper_pair_term(a, b, m, nactive, z, L::Integer, epsv, nn::Integer)
    P, D = _logratio_upper_pair_state(a, b, m, nactive, z, L, epsv, nn)
    return P * D * _logratio_exprel(epsv * D)
end

function _logratio_upper_pair_state_limit_direct(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    s = a[1] - 1 - K
    logH, signH = _logratio_upper_logH(a, b, m, nactive, s)
    logX = logH -
           loggamma(_logratio_int_as(s, K + 1)) - _logratio_logpoch(one(s), nn)
    D = _logratio_upper_D_limit_direct(a, b, m, nactive, z, L, nn)
    return signH * (-1)^L * z^s * exp(logX), D
end

function _logratio_upper_pair_term_limit(a, b, m, nactive, z, L::Integer, nn::Integer)
    P, D = _logratio_upper_pair_state_limit(a, b, m, nactive, z, L, nn)
    return P * D
end

# Ratio forms of the exact-pair state and the ordinary-family prefactor. By reflection every
# gamma factor at s = b1 + K becomes Γ(s+1-a_j) or Γ(s+1-b_j) up to K-independent constants:
#   log X = Σ_j log Γ(s+1-a_j) - Σ_j log Γ(s+1-b_j) + Σ_{j=3}^{m} [log π - log sinpi(b_j - b1)]
#           - Σ_{j>n} [log π - log sinpi(a_j - b1)]   (mod 2πi),
#   D     = Σ_j ψ(s+1-a_j) - Σ_j ψ(s+1-b_j) + log z
#           + Σ_{j=3}^{m} π cot π(b_j - b1) - Σ_{j>n} π cot π(a_j - b1).
# The upper branch mirrors this with t = K + 1 - a1. The sums are evaluated pairwise by
# `_loggamma_sum_shifted` / `_digamma_sum_shifted`, so the rounding error scales with the result.
# Each form returns `nothing` at a singular reflection constant and the caller evaluates directly.

_logratio_finite_or_nothing(x) = _logratio_finite(x) ? x : nothing

# Restore the real type of real-parameter inputs after complex log-space work
# (the grouped pair states are complex, as in the direct evaluation).
_logratio_like(x, s) = s isa Real ? real(x) : x

_logratio_logpi(x) = log(oftype(float(real(x)), pi))

function _logratio_lower_pair_offsets(a, b, L)
    b1 = b[1]
    plus = map(aj -> 1 - aj + b1, a)
    minus = ntuple(length(b)) do j
        j == 1 ? one(b1) : j == 2 ? one(b1) - L : 1 - b[j] + b1
    end
    return plus, minus
end

function _logratio_lower_logX_ratio(a, b, m, nactive, L::Integer, K::Integer)
    b1 = b[1]
    plus, minus = _logratio_lower_pair_offsets(a, b, L)
    out = _loggamma_sum_shifted(K, plus, minus)
    logpi = _logratio_logpi(b1)
    flips = 0
    for j in 3:m
        out += logpi - log(complex(sinpi(b[j] - b1)))
        flips += 1
    end
    for j in (nactive + 1):length(a)
        out -= logpi - log(complex(sinpi(a[j] - b1)))
        flips += 1
    end
    isodd(K * flips) && (out += im * oftype(real(out), pi))
    return _logratio_finite_or_nothing(out)
end

function _logratio_lower_D_limit(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    b1 = b[1]
    plus, minus = _logratio_lower_pair_offsets(a, b, L)
    D = _digamma_sum_shifted(K, plus, minus) + log(z)
    piT = oftype(float(real(b1)), pi)
    for j in 3:m
        D += piT * cospi(b[j] - b1) / sinpi(b[j] - b1)
    end
    for j in (nactive + 1):length(a)
        D -= piT * cospi(a[j] - b1) / sinpi(a[j] - b1)
    end
    _logratio_finite(D) || return _logratio_lower_D_limit_direct(a, b, m, nactive, z, L, nn)
    return D
end

function _logratio_lower_pair_state_limit(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    logX = _logratio_lower_logX_ratio(a, b, m, nactive, L, K)
    logX === nothing && return _logratio_lower_pair_state_limit_direct(a, b, m, nactive, z, L, nn)
    s = b[1] + K
    P = (-1)^(L + 1) * z^s * exp(logX)
    return P, _logratio_lower_D_limit(a, b, m, nactive, z, L, nn)
end

function _logratio_upper_pair_offsets(a, b, L)
    a1 = a[1]
    plus = map(bj -> bj - a1 + 1, b)
    minus = ntuple(length(a)) do j
        j == 1 ? one(a1) : j == 2 ? one(a1) - L : a[j] - a1 + 1
    end
    return plus, minus
end

function _logratio_upper_logX_ratio(a, b, m, nactive, L::Integer, K::Integer)
    a1 = a[1]
    plus, minus = _logratio_upper_pair_offsets(a, b, L)
    out = _loggamma_sum_shifted(K, plus, minus)
    logpi = _logratio_logpi(a1)
    flips = 0
    for j in 3:nactive
        out += logpi - log(complex(sinpi(a[j] - a1)))
        flips += 1
    end
    for j in (m + 1):length(b)
        out -= logpi - log(complex(sinpi(b[j] - a1)))
        flips += 1
    end
    isodd((K + 1) * flips) && (out += im * oftype(real(out), pi))
    return _logratio_finite_or_nothing(out)
end

function _logratio_upper_D_limit(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    a1 = a[1]
    plus, minus = _logratio_upper_pair_offsets(a, b, L)
    # Σ_a ψ(t + a_j) - Σ_b ψ(t + b_j): the a-offsets are the positive set here.
    D = _digamma_sum_shifted(K, minus, plus) + log(z)
    piT = oftype(float(real(a1)), pi)
    for j in 3:nactive
        D += piT * cospi(a[j] - a1) / sinpi(a[j] - a1)
    end
    for j in (m + 1):length(b)
        D -= piT * cospi(b[j] - a1) / sinpi(b[j] - a1)
    end
    _logratio_finite(D) || return _logratio_upper_D_limit_direct(a, b, m, nactive, z, L, nn)
    return D
end

function _logratio_upper_pair_state_limit(a, b, m, nactive, z, L::Integer, nn::Integer)
    K = L + nn
    logX = _logratio_upper_logX_ratio(a, b, m, nactive, L, K)
    logX === nothing && return _logratio_upper_pair_state_limit_direct(a, b, m, nactive, z, L, nn)
    s = a[1] - 1 - K
    P = (-1)^L * z^s * exp(logX)
    return P, _logratio_upper_D_limit(a, b, m, nactive, z, L, nn)
end

"""
    _logratio_family_prefactor(a, b, m, nactive, k, branch)

Gamma prefactor of the ordinary Slater family `k` in ratio form, or `nothing`
when a reflection constant is singular.
"""
# Tuple `t` without its `k`-th entry, typed (a filtered generator collects through `Vector{Any}`).
_logratio_drop_index(t::Tuple, k::Integer) =
    ntuple(i -> t[i < k ? i : i + 1], length(t) - 1)

function _logratio_family_prefactor(a, b, m, nactive, k::Integer, branch::Symbol)
    p = length(a)
    q = length(b)
    if branch === :lower
        c = b[k]
        plus = map(aj -> 1 + c - aj, a)
        minus = _logratio_drop_index(map(bj -> 1 + c - bj, b), k)
        out = _loggamma_sum_shifted(0, plus, minus)
        logpi = _logratio_logpi(c)
        for j in 1:m
            j == k && continue
            out += logpi - log(complex(sinpi(b[j] - c)))
        end
        for j in (nactive + 1):p
            out -= logpi - log(complex(sinpi(a[j] - c)))
        end
    else
        c = a[k]
        plus = map(bj -> 1 - c + bj, b)
        minus = _logratio_drop_index(map(aj -> 1 - c + aj, a), k)
        out = _loggamma_sum_shifted(0, plus, minus)
        logpi = _logratio_logpi(c)
        for j in 1:nactive
            j == k && continue
            out += logpi - log(complex(sinpi(a[j] - c + 1)))
        end
        for j in (m + 1):q
            out -= logpi - log(complex(sinpi(b[j] - c + 1)))
        end
    end
    _logratio_finite(out) || return nothing
    return _logratio_like(exp(out), c)
end

function _logratio_lower_raw_cancellation(a, b, m, nactive, z, L)
    t1 = _logratio_lower_residue_coeff(a, b, m, nactive, 1, L) * z^(b[1] + L)
    t2 = _logratio_lower_residue_coeff(a, b, m, nactive, 2, 0) * z^b[2]
    return (abs(t1) + abs(t2)) / max(abs(t1 + t2), eps(Float64))
end

function _logratio_upper_raw_cancellation(a, b, m, nactive, z, L, epsv)
    t1 = _logratio_upper_unpaired_residue(a, b, m, nactive, z, 1, L)
    s2 = a[2] - 1
    logH, signH = _logratio_upper_logH(a, b, m, nactive, s2)
    t2 = gamma(-L + epsv) * signH * exp(logH) * z^s2
    return (abs(t1) + abs(t2)) / max(abs(t1 + t2), eps(Float64))
end

"""
    _lower_logratio_jet_exact_series(a, b, m, n, z, L; maxterms, tol)

Evaluate the lower exact-confluence log-ratio limit series for a canonical
active pair.  This is the epsilon-zero limit of `_lower_logratio_jet_series`
and is the exact-confluence path used by `meijerg`.

The grouped state `(P, D)` advances by exact one-step recurrences (gamma
recurrence for `P`, digamma recurrence for `D`) and is recomputed from
scratch every `refresh_every` terms to bound rounding drift.
"""
function _lower_logratio_jet_exact_series(a::Tuple, b::Tuple, m::Integer, nactive::Integer,
                                          z, L::Integer; maxterms=10_000,
                                          tol=_logratio_exact_tol(z),
                                          refresh_every::Integer=_LOGRATIO_JET_REFRESH)
    total = zero(z + 0im)
    corr = zero(total)
    absum = zero(abs(total))
    ordinary = 0
    for ell in 0:(L - 1)
        # `oftype` keeps the accumulator type inferable (pFq's return type currently is not).
        v = oftype(total, _logratio_lower_residue_coeff(a, b, m, nactive, 1, ell) * z^(b[1] + ell))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    for k in 3:m
        v = oftype(total, _logratio_lower_single_family_term(a, b, m, nactive, z, k))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    p = length(a)
    q = length(b)
    P, D = _logratio_lower_pair_state_limit(a, b, m, nactive, z, L, 0)
    tail_constants = _logratio_tail_constants(a, b, L, zero(D), :lower)
    status = :terms_exceeded
    tail = Inf
    tail_ratio = Inf
    used = 0
    for nn in 0:(maxterms - 1)
        K = L + nn
        s = b[1] + K
        if nn > 0
            if nn % refresh_every == 0
                P, D = _logratio_lower_pair_state_limit(a, b, m, nactive, z, L, nn)
            else
                sold = s - 1
                ratio = z / (K * nn)
                for j in 3:m
                    ratio /= (b[j] - s)
                    D += inv(b[j] - s)
                end
                for j in 1:nactive
                    ratio *= (1 - a[j] + sold)
                    D += inv(1 - a[j] + sold)
                end
                for j in (m + 1):q
                    ratio /= (1 + sold - b[j])
                    D -= inv(1 + sold - b[j])
                end
                for j in (nactive + 1):p
                    ratio *= (a[j] - s)
                    D -= inv(a[j] - s)
                end
                D -= inv(_logratio_int_as(z, nn))
                D -= inv(_logratio_int_as(z, K))
                P *= ratio
            end
        end
        term = P * D
        absum += abs(term)
        total, corr = _logratio_sum2(total, corr, term)
        used = nn + 1
        if !_logratio_finite(total) || !_logratio_finite(term)
            status = :nonfinite
            break
        end
        cert = _logratio_pair_tail_bound(P, D, a, b, z, L, zero(D), nn, :lower;
                                         constants=tail_constants)
        tail = cert.tail_bound
        tail_ratio = cert.theta
        if tail <= tol * max(abs(total), one(abs(total)))
            status = :ok
            break
        end
    end
    return (status=status, value=total + corr, terms_used=used, grouped_terms_used=used,
            ordinary_terms_used=ordinary, tail_estimate=tail,
            absolute_sum=absum,
            tail_ratio_bound=tail_ratio,
            tail_status=status === :ok ? :accepted : :not_accepted,
            convergence_status=status === :ok ? :ok : status,
            primitive_status=:ok)
end

"""
    _lower_logratio_jet_series(a, b, m, n, z, L, epsilon; maxterms, tol)

Evaluate the selected lower near-confluent grouped Slater stream.  Ordinary
unpaired active families are added directly, while the canonical near pair is
summed as log-ratio grouped coefficients with Kahan compensation and the
cancellation-free recurrence tail certificate.

The grouped state `(P, D)` advances by exact one-step recurrences (gamma
recurrence for `P`, `log1p(±epsilon/x)/epsilon` increments for `D`) and is
recomputed from scratch every `refresh_every` terms to bound rounding drift.
"""
function _lower_logratio_jet_series(a::Tuple, b::Tuple, m::Integer, nactive::Integer,
                                    z, L::Integer, epsv; maxterms=10_000,
                                    tol=_logratio_series_tol(z),
                                    refresh_every::Integer=_LOGRATIO_JET_REFRESH)
    total = zero(z + 0im)
    corr = zero(total)
    absum = zero(abs(total))
    ordinary = 0
    for ell in 0:(L - 1)
        v = oftype(total, _logratio_lower_residue_coeff(a, b, m, nactive, 1, ell) * z^(b[1] + ell))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    for k in 3:m
        v = oftype(total, _logratio_lower_single_family_term(a, b, m, nactive, z, k))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    p = length(a)
    q = length(b)
    inv_eps = inv(epsv)
    P, D = _logratio_lower_pair_state(a, b, m, nactive, z, L, epsv, 0)
    tail_constants = _logratio_tail_constants(a, b, L, epsv, :lower)
    status = :terms_exceeded
    tail = Inf
    tail_ratio = Inf
    used = 0
    for nn in 0:(maxterms - 1)
        K = L + nn
        s = b[1] + K
        if nn > 0
            if nn % refresh_every == 0
                P, D = _logratio_lower_pair_state(a, b, m, nactive, z, L, epsv, nn)
            else
                sold = s - 1
                ratio = z / (K * (nn - epsv))
                for j in 3:m
                    ratio /= (b[j] - s)
                    D -= _logratio_log1p(-epsv / (b[j] - s)) * inv_eps
                end
                for j in 1:nactive
                    ratio *= (1 - a[j] + sold)
                    D += _logratio_log1p(epsv / (1 - a[j] + sold)) * inv_eps
                end
                for j in (m + 1):q
                    ratio /= (1 + sold - b[j])
                    D -= _logratio_log1p(epsv / (1 + sold - b[j])) * inv_eps
                end
                for j in (nactive + 1):p
                    ratio *= (a[j] - s)
                    D += _logratio_log1p(-epsv / (a[j] - s)) * inv_eps
                end
                D += _logratio_log1p(-epsv / nn) * inv_eps
                D -= _logratio_log1p(epsv / K) * inv_eps
                P *= ratio
            end
        end
        term = P * D * _logratio_exprel(epsv * D)
        absum += abs(term)
        total, corr = _logratio_sum2(total, corr, term)
        used = nn + 1
        if !_logratio_finite(total) || !_logratio_finite(term)
            status = :nonfinite
            break
        end
        cert = _logratio_pair_tail_bound(P, D, a, b, z, L, epsv, nn, :lower;
                                         constants=tail_constants)
        tail = cert.tail_bound
        tail_ratio = cert.theta
        if tail <= tol * max(abs(total), one(abs(total)))
            status = :ok
            break
        end
    end
    return (status=status, value=total + corr, terms_used=used, grouped_terms_used=used,
            ordinary_terms_used=ordinary, tail_estimate=tail,
            absolute_sum=absum,
            tail_ratio_bound=tail_ratio,
            tail_status=status === :ok ? :accepted : :not_accepted,
            convergence_status=status === :ok ? :ok : status,
            primitive_status=:ok)
end

"""
    _upper_logratio_jet_exact_series(a, b, m, n, z, L; maxterms, tol)

Evaluate the upper exact-confluence log-ratio limit series for a canonical
active pair.  This is the infinity-branch analogue of
`_lower_logratio_jet_exact_series`, with the same recurrence-based grouped
state advance.
"""
function _upper_logratio_jet_exact_series(a::Tuple, b::Tuple, m::Integer, nactive::Integer,
                                          z, L::Integer; maxterms=10_000,
                                          tol=_logratio_exact_tol(z),
                                          refresh_every::Integer=_LOGRATIO_JET_REFRESH)
    total = zero(z + 0im)
    corr = zero(total)
    absum = zero(abs(total))
    ordinary = 0
    for ell in 0:(L - 1)
        v = oftype(total, _logratio_upper_unpaired_residue(a, b, m, nactive, z, 1, ell))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    for h in 3:nactive
        v = oftype(total, _logratio_upper_single_family_term(a, b, m, nactive, z, h))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    p = length(a)
    q = length(b)
    invz = inv(z)
    P, D = _logratio_upper_pair_state_limit(a, b, m, nactive, z, L, 0)
    tail_constants = _logratio_tail_constants(a, b, L, zero(D), :upper)
    status = :terms_exceeded
    tail = Inf
    tail_ratio = Inf
    used = 0
    for nn in 0:(maxterms - 1)
        K = L + nn
        s = a[1] - 1 - K
        if nn > 0
            if nn % refresh_every == 0
                P, D = _logratio_upper_pair_state_limit(a, b, m, nactive, z, L, nn)
            else
                sold = s + 1
                ratio = invz / (K * nn)
                for j in 1:m
                    ratio *= (b[j] - sold)
                    D -= inv(b[j] - sold)
                end
                for j in 3:nactive
                    ratio /= (1 - a[j] + s)
                    D -= inv(1 - a[j] + s)
                end
                for j in (m + 1):q
                    ratio *= (1 - b[j] + s)
                    D += inv(1 - b[j] + s)
                end
                for j in (nactive + 1):p
                    ratio /= (a[j] - sold)
                    D += inv(a[j] - sold)
                end
                D += inv(_logratio_int_as(z, K))
                D += inv(_logratio_int_as(z, nn))
                P *= ratio
            end
        end
        term = P * D
        absum += abs(term)
        total, corr = _logratio_sum2(total, corr, term)
        used = nn + 1
        if !_logratio_finite(total) || !_logratio_finite(term)
            status = :nonfinite
            break
        end
        cert = _logratio_pair_tail_bound(P, D, a, b, z, L, zero(D), nn, :upper;
                                         constants=tail_constants)
        tail = cert.tail_bound
        tail_ratio = cert.theta
        if tail <= tol * max(abs(total), one(abs(total)))
            status = :ok
            break
        end
    end
    return (status=status, value=total + corr, terms_used=used, grouped_terms_used=used,
            ordinary_terms_used=ordinary, tail_estimate=tail,
            absolute_sum=absum,
            tail_ratio_bound=tail_ratio,
            tail_status=status === :ok ? :accepted : :not_accepted,
            convergence_status=status === :ok ? :ok : status,
            primitive_status=:ok)
end

"""
    _upper_logratio_jet_series(a, b, m, n, z, L, epsilon; maxterms, tol)

Evaluate the selected upper near-confluent grouped Slater stream.  This is
the infinity-branch analogue of `_lower_logratio_jet_series`, with the upper
pole convention `s = a[1] - 1 - K` and the same recurrence-based grouped
state advance.
"""
function _upper_logratio_jet_series(a::Tuple, b::Tuple, m::Integer, nactive::Integer,
                                    z, L::Integer, epsv; maxterms=10_000,
                                    tol=_logratio_series_tol(z),
                                    refresh_every::Integer=_LOGRATIO_JET_REFRESH)
    total = zero(z + 0im)
    corr = zero(total)
    absum = zero(abs(total))
    ordinary = 0
    for ell in 0:(L - 1)
        v = oftype(total, _logratio_upper_unpaired_residue(a, b, m, nactive, z, 1, ell))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    for h in 3:nactive
        v = oftype(total, _logratio_upper_single_family_term(a, b, m, nactive, z, h))
        absum += abs(v)
        total, corr = _logratio_sum2(total, corr, v)
        ordinary += 1
    end
    p = length(a)
    q = length(b)
    inv_eps = inv(epsv)
    invz = inv(z)
    P, D = _logratio_upper_pair_state(a, b, m, nactive, z, L, epsv, 0)
    tail_constants = _logratio_tail_constants(a, b, L, epsv, :upper)
    status = :terms_exceeded
    tail = Inf
    tail_ratio = Inf
    used = 0
    for nn in 0:(maxterms - 1)
        K = L + nn
        s = a[1] - 1 - K
        if nn > 0
            if nn % refresh_every == 0
                P, D = _logratio_upper_pair_state(a, b, m, nactive, z, L, epsv, nn)
            else
                sold = s + 1
                ratio = invz / (K * (nn + epsv))
                for j in 1:m
                    ratio *= (b[j] - sold)
                    D += _logratio_log1p(-epsv / (b[j] - sold)) * inv_eps
                end
                for j in 3:nactive
                    ratio /= (1 - a[j] + s)
                    D -= _logratio_log1p(epsv / (1 - a[j] + s)) * inv_eps
                end
                for j in (m + 1):q
                    ratio *= (1 - b[j] + s)
                    D += _logratio_log1p(epsv / (1 - b[j] + s)) * inv_eps
                end
                for j in (nactive + 1):p
                    ratio /= (a[j] - sold)
                    D -= _logratio_log1p(-epsv / (a[j] - sold)) * inv_eps
                end
                D -= _logratio_log1p(-epsv / K) * inv_eps
                D += _logratio_log1p(epsv / nn) * inv_eps
                P *= ratio
            end
        end
        term = P * D * _logratio_exprel(epsv * D)
        absum += abs(term)
        total, corr = _logratio_sum2(total, corr, term)
        used = nn + 1
        if !_logratio_finite(total) || !_logratio_finite(term)
            status = :nonfinite
            break
        end
        cert = _logratio_pair_tail_bound(P, D, a, b, z, L, epsv, nn, :upper;
                                         constants=tail_constants)
        tail = cert.tail_bound
        tail_ratio = cert.theta
        if tail <= tol * max(abs(total), one(abs(total)))
            status = :ok
            break
        end
    end
    return (status=status, value=total + corr, terms_used=used, grouped_terms_used=used,
            ordinary_terms_used=ordinary, tail_estimate=tail,
            absolute_sum=absum,
            tail_ratio_bound=tail_ratio,
            tail_status=status === :ok ? :accepted : :not_accepted,
            convergence_status=status === :ok ? :ok : status,
            primitive_status=:ok)
end

"""
    _logratio_jet_policy(a, b, m, n, z, branch, pair)

Apply the conservative first-inclusion policy for the private near-confluent
candidate.  The policy accepts only lower/upper, finite non-balanced,
canonical active pairs with small `epsilon`, controlled `epsilon*log(z)`,
validated `z`, and branch-specific safety gates.
"""
function _logratio_jet_policy(a::Tuple, b::Tuple, m::Integer, n::Integer, z,
                              branch::Symbol, pair)
    p = length(a)
    q = length(b)
    p == q && return (selected=false, status=:unsupported_finite_p_eq_q,
                      reason=:unsupported_finite_p_eq_q, diagnostics=(pair=pair,))
    branch in (:lower, :upper) || return (selected=false, status=:unsupported_branch,
                                          reason=:unsupported_branch, diagnostics=(pair=pair,))
    _logratio_validated_z(branch, z) || return (selected=false, status=:not_near_confluent,
                                                reason=:outside_validated_z_mode,
                                                diagnostics=(pair=pair,))
    epsv = pair.epsilon
    abs(epsv) > 0 || return (selected=false, status=:exact_limit_sufficient,
                             reason=:exact_confluence, diagnostics=(pair=pair,))
    abs(epsv) <= _LOGRATIO_JET_MAX_ABS_EPSILON ||
        return (selected=false, status=:not_near_confluent, reason=:epsilon_too_large,
                diagnostics=(pair=pair,))
    abs(epsv * log(z)) <= _LOGRATIO_JET_MAX_ABS_EPSILON_LOGZ ||
        return (selected=false, status=:not_near_confluent, reason=:epsilon_logz_too_large,
                diagnostics=(pair=pair,))

    cancel = branch === :lower ?
        _logratio_lower_raw_cancellation(a, b, m, n, z, pair.L) :
        _logratio_upper_raw_cancellation(a, b, m, n, z, pair.L, epsv)
    if branch === :lower && !(isfinite(cancel) && cancel >= _LOGRATIO_JET_MIN_LOWER_RAW_CANCELLATION)
        return (selected=false, status=:raw_slater_safe, reason=:raw_slater_safe,
                diagnostics=(pair=pair, raw_cancellation_ratio=cancel))
    end
    return (selected=true, status=:ok, reason=:selected,
            diagnostics=(pair=pair, raw_cancellation_ratio=cancel))
end

_logratio_widen(x::Float32) = Float64(x)
_logratio_widen(x::Complex{Float32}) = ComplexF64(x)
_logratio_widen(x) = x

# Convergent regime of each branch: lower small |z|, upper large |z|.
_logratio_in_regime(branch::Symbol, z) = branch === :lower ? abs(z) <= 1 : abs(z) >= 1

# Working-precision modes of the jet entry points, see `meijerg_jet`. `:auto` is the
# behaviour of `meijerg`.
const _JET_WORKING_PRECISIONS = (:auto, :input, :extended)
function _check_working_precision(mode::Symbol)
    mode in _JET_WORKING_PRECISIONS ||
        throw(ArgumentError("working_precision must be :auto, :input or :extended"))
    return mode
end

# Re-evaluation of a native pair-jet series result `diag` with conditioning `kappa`.
# Returns the result to keep, or `nothing` when the cancellation cannot be absorbed.
# `:extended` always re-evaluates BigFloat inputs at least at doubled precision.
function _logratio_reevaluate(run, args::Tuple, z, diag, kappa, mode::Symbol)
    force = mode === :extended
    (force || !_logratio_summation_ok(kappa)) || return diag
    hi = mode === :input ? nothing :
         _logratio_escalated(run, args, z, kappa;
                             min_extra=force ? precision(BigFloat) : 0)
    hi === nothing && return _logratio_summation_ok(kappa) ? diag : nothing
    return hi
end


# Statuses of a detected but declined confluent geometry. `meijerg` skips raw Slater on these
# (the separated expansion is exactly what cancels there) and continues to contour integration.
const _LOGRATIO_DECLINED_STATUSES = (:summation_cancellation, :tail_not_accepted,
                                     :terms_exceeded, :nonfinite, :internal_error)

# Acceptance of a double-double lane result `ext = (status, value, kappa)`.
function _logratio_lane_result(ext, pair, br, return_diagnostics)
    diagnostics = return_diagnostics ?
        (status=ext.status, value=ext.value, pair=pair, lane=:extended,
         summation_kappa=ext.kappa) : (;)
    ext.status === :ok ||
        return _logratio_result(ext.status, ext.value, br, false, ext.status, diagnostics)
    ext.kappa <= _LOGRATIO_JET_EXTENDED_MAX_KAPPA ||
        return _logratio_result(:summation_cancellation, ext.value, br, false,
                                :summation_cancellation, diagnostics)
    _logratio_finite(ext.value) ||
        return _logratio_result(:nonfinite, ext.value, br, false, :nonfinite, diagnostics)
    return _logratio_result(:ok, ext.value, br, true, :selected, diagnostics)
end

_logratio_real_type(z, a, b) = typeof(float(real(first(promote(z, a..., b...)))))

"""
    _try_meijerg_slater_logratio_jet_exact(a, b, m, n, z; branch=:auto)

Try the private exact-confluence log-ratio jet evaluator.  It selects only a
single active exact pair, canonicalizes the pair order internally, and returns
a structured result.  Unsupported exact configurations are non-selected so
`meijerg` can continue through Slater and contour fallback.

With `working_precision=:auto`, a series whose summation conditioning exceeds
`_LOGRATIO_JET_MAX_SUMMATION_KAPPA` is re-evaluated in doubled working precision:
Float32 inputs on the Float64 path, Float64 inputs in the double-double lane
(`_logratio_jet_exact_extended`, taken directly off the convergent regime), and
BigFloat inputs with raised precision.  With `:input` such a series is declined, and
`:extended` always takes the doubled precision (see `meijerg_jet`).
"""
function _try_meijerg_slater_logratio_jet_exact(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
        branch::Symbol=:auto, maxterms::Integer=10_000, return_diagnostics::Bool=false,
        working_precision::Symbol=:auto)
    extended = working_precision !== :input
    !_logratio_jet_enabled() &&
        return _logratio_result(:not_near_confluent, _logratio_nan_value(z), :unknown,
                                false, :disabled, (;))
    try
        p = length(a)
        q = length(b)
        p == q && return _logratio_result(:unsupported_finite_p_eq_q, _logratio_nan_value(z),
                                          :unknown, false, :unsupported_finite_p_eq_q, (;))
        br = branch === :auto ? (p < q ? :lower : :upper) : branch
        br in (:lower, :upper) || return _logratio_result(:unsupported_branch,
                                                          _logratio_nan_value(z), br,
                                                          false, :unsupported_branch, (;))
        det = _detect_exact_logratio_jet_pair(a, b, m, n, br)
        det.status === :ok || return _logratio_result(det.status, _logratio_nan_value(z),
                                                      br, false, det.status,
                                                      return_diagnostics ? (pair=det.pair,) : (;))
        R = _logratio_real_type(z, a, b)
        # Float32 uses the Float64 path. Widen the canonical parameters with the pair kept exactly
        # confluent, since widening the raw inputs would turn an exact pair into a near pair.
        if extended && R === Float32
            aw = map(_logratio_widen, det.a)
            bw = map(_logratio_widen, det.b)
            L = det.pair.L
            br === :lower ? (bw = (bw[1], bw[1] + L, bw[3:end]...)) :
                            (aw = (aw[1], aw[1] - L, aw[3:end]...))
            r = _try_meijerg_slater_logratio_jet_exact(aw, bw, m, n, _logratio_widen(z);
                    branch=br, maxterms, return_diagnostics, working_precision=:auto)
            return merge(r, (value=convert(Complex{Float32}, r.value),))
        end
        lane = extended && R === Float64
        extended_lane() = _logratio_lane_result(
            _logratio_jet_exact_extended(det.a, det.b, m, n, z, det.pair.L, br; maxterms),
            det.pair, br, return_diagnostics)
        # Off the convergent regime the series cancels (κ ~ exp(ν w^(1/ν))), so Float64
        # inputs go straight to the double-double lane.
        lane && (working_precision === :extended || !_logratio_in_regime(br, z)) &&
            return extended_lane()
        # The exact states are refreshed in ratio form (`_logratio_*_pair_state_limit`),
        # accurate enough that the gamma recurrence needs no per-term refresh.
        refresh_every = _LOGRATIO_JET_REFRESH
        series = br === :lower ? _lower_logratio_jet_exact_series :
                                 _upper_logratio_jet_exact_series
        run(aa, bb, zz) = series(aa, bb, m, n, zz, det.pair.L;
                                 maxterms=maxterms, refresh_every=refresh_every)
        diag = run(det.a, det.b, z)
        kappa = _logratio_summation_kappa(diag.absolute_sum, diag.value)
        if diag.status === :ok
            # Cancellation between terms: re-evaluate in doubled working precision
            # where allowed, otherwise decline so `meijerg` falls back to contour.
            lane && !_logratio_summation_ok(kappa) && return extended_lane()
            hi = _logratio_reevaluate(run, (det.a, det.b, z), z, diag, kappa, working_precision)
            hi === nothing &&
                return _logratio_result(:summation_cancellation, diag.value, br, false,
                                        :summation_cancellation,
                                        return_diagnostics ? merge(diag, (pair=det.pair,
                                            refresh_every=refresh_every,
                                            summation_kappa=kappa)) : (;))
            diag = hi
        end
        diagnostics = return_diagnostics ? merge(diag, (pair=det.pair,
                                                        refresh_every=refresh_every,
                                                        summation_kappa=kappa)) : (;)
        diag.status === :ok ||
            return _logratio_result(diag.status, diag.value, br, false, diag.status, diagnostics)
        diag.tail_status === :accepted ||
            return _logratio_result(:tail_not_accepted, diag.value, br, false, :tail_not_accepted, diagnostics)
        _logratio_finite(diag.value) ||
            return _logratio_result(:nonfinite, diag.value, br, false, :nonfinite, diagnostics)
        return _logratio_result(:ok, diag.value, br, true, :selected, diagnostics)
    catch e
        return _logratio_result(:internal_error, _logratio_nan_value(z),
                                branch === :auto ? :unknown : branch,
                                false, :internal_error, (error=sprint(showerror, e),))
    end
end

"""
    _try_meijerg_slater_logratio_jet(a, b, m, n, z; branch=:auto)

Try the private log-ratio jet evaluator for active-pole confluence.  This is
the single internal entry point used by `meijerg`.  It evaluates nonzero
separations with the finite-epsilon log-ratio form and does not snap near
points to an exact-confluence limit.  Literal exact-zero confluence is handled
by the exact-zero limit of the same log-ratio representation. Unsupported,
ambiguous, or failed cases return non-selected statuses so the caller can
continue through the existing fallback sequence. A status in
`_LOGRATIO_DECLINED_STATUSES` reports a detected confluent geometry whose
series was declined, so that the caller skips the separated Slater expansion.
"""
function _try_meijerg_slater_logratio_jet(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
        branch::Symbol=:auto, maxterms::Integer=10_000, return_diagnostics::Bool=false,
        working_precision::Symbol=:auto)
    near = _try_meijerg_slater_logratio_jet_near(a, b, m, n, z;
        branch=branch, maxterms=maxterms, return_diagnostics=return_diagnostics,
        working_precision)
    near.status === :ok && near.selected && return near
    near.reason === :disabled && return near

    exact = _try_meijerg_slater_logratio_jet_exact(a, b, m, n, z;
        branch=branch, maxterms=maxterms, return_diagnostics=return_diagnostics,
        working_precision)
    exact.status === :ok && exact.selected && return exact

    higher = _try_meijerg_slater_higher_logratio_jet(a, b, m, n, z;
        branch=branch, return_diagnostics=return_diagnostics,
        working_precision)
    higher.status === :ok && higher.selected && return higher
    for r in (near, exact, higher)
        r.status in _LOGRATIO_DECLINED_STATUSES && return r
    end
    return near
end

"""
    _try_meijerg_slater_logratio_jet_near(a, b, m, n, z; branch=:auto)

Try the private near-confluent log-ratio jet evaluator.  It returns a
structured internal result and selects only validated punctured neighborhoods.
Unsupported, exact, ambiguous, or failed cases return non-selected statuses so
`meijerg` can continue through the existing fallback sequence.

The summation-conditioning policy is that of
`_try_meijerg_slater_logratio_jet_exact`: with `working_precision=:auto` a
Float64 series off the convergent regime, or one whose `κ_sum` exceeds
`_LOGRATIO_JET_MAX_SUMMATION_KAPPA`, is evaluated in the double-double lane
`_logratio_jet_near_extended`. Float32 inputs use the Float64 path and BigFloat
inputs raised precision. `:extended` always takes the doubled precision.
"""
function _try_meijerg_slater_logratio_jet_near(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
        branch::Symbol=:auto, maxterms::Integer=10_000, return_diagnostics::Bool=false,
        working_precision::Symbol=:auto)
    extended = working_precision !== :input
    !_logratio_jet_enabled() &&
        return _logratio_result(:not_near_confluent, _logratio_nan_value(z), :unknown,
                                false, :disabled, (;))
    try
        p = length(a)
        q = length(b)
        p == q && return _logratio_result(:unsupported_finite_p_eq_q, _logratio_nan_value(z),
                                          :unknown, false, :unsupported_finite_p_eq_q, (;))
        br = branch === :auto ? (p < q ? :lower : :upper) : branch
        br in (:lower, :upper) || return _logratio_result(:unsupported_branch,
                                                          _logratio_nan_value(z), br,
                                                          false, :unsupported_branch, (;))
        det = _detect_near_logratio_jet_pair(a, b, m, n, br)
        det.status === :ok || return _logratio_result(det.status, _logratio_nan_value(z),
                                                      br, false, det.status,
                                                      return_diagnostics ? (pair=det.pair,) : (;))
        policy = _logratio_jet_policy(det.a, det.b, m, n, z, br, det.pair)
        policy.selected || return _logratio_result(policy.status, _logratio_nan_value(z),
                                                   br, false, policy.reason,
                                                   return_diagnostics ? policy.diagnostics : (;))
        R = _logratio_real_type(z, a, b)
        # Float32: the Float64 path is the doubled working precision.  Widening is exact
        # and the canonical pair re-detects with the same (L, ε).
        if extended && R === Float32
            r = _try_meijerg_slater_logratio_jet_near(map(_logratio_widen, det.a),
                    map(_logratio_widen, det.b), m, n, _logratio_widen(z);
                    branch=br, maxterms, return_diagnostics, working_precision=:auto)
            return merge(r, (value=convert(Complex{Float32}, r.value),))
        end
        lane = extended && R === Float64
        extended_lane() = _logratio_lane_result(
            _logratio_jet_near_extended(det.a, det.b, m, n, z, det.pair.L, det.pair.epsilon,
                                        br; maxterms),
            det.pair, br, return_diagnostics)
        lane && (working_precision === :extended || !_logratio_in_regime(br, z)) &&
            return extended_lane()
        refresh_every = _logratio_refresh_every(br, z)
        series = br === :lower ? _lower_logratio_jet_series : _upper_logratio_jet_series
        epsv = det.pair.epsilon
        run(aa, bb, zz) = series(aa, bb, m, n, zz, det.pair.L, epsv;
                                 maxterms=maxterms, refresh_every=refresh_every)
        diag = run(det.a, det.b, z)
        kappa = _logratio_summation_kappa(diag.absolute_sum, diag.value)
        if diag.status === :ok
            lane && !_logratio_summation_ok(kappa) && return extended_lane()
            hi = _logratio_reevaluate(run, (det.a, det.b, z), z, diag, kappa, working_precision)
            hi === nothing &&
                return _logratio_result(:summation_cancellation, diag.value, br, false,
                                        :summation_cancellation,
                                        return_diagnostics ? merge(diag, policy.diagnostics,
                                            (refresh_every=refresh_every,
                                             summation_kappa=kappa)) : (;))
            diag = hi
        end
        diagnostics = return_diagnostics ? merge(diag, policy.diagnostics,
                                                 (refresh_every=refresh_every,
                                                  summation_kappa=kappa)) : (;)
        diag.terms_used <= _LOGRATIO_JET_MAX_GROUPED_TERMS ||
            return _logratio_result(:terms_exceeded, diag.value, br, false, :terms_exceeded, diagnostics)
        diag.status === :ok ||
            return _logratio_result(diag.status, diag.value, br, false, diag.status, diagnostics)
        diag.tail_status === :accepted ||
            return _logratio_result(:tail_not_accepted, diag.value, br, false, :tail_not_accepted, diagnostics)
        _logratio_finite(diag.value) ||
            return _logratio_result(:nonfinite, diag.value, br, false, :nonfinite, diagnostics)
        return _logratio_result(:ok, diag.value, br, true, :selected, diagnostics)
    catch e
        return _logratio_result(:internal_error, _logratio_nan_value(z),
                                branch === :auto ? :unknown : branch,
                                false, :internal_error, (error=sprint(showerror, e),))
    end
end
