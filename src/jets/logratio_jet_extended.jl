# Double-double lanes for the pair log-ratio jets. Off the convergent regime the grouped series
# cancels by κ = Σ|T_n|/|ΣT_n|, so the whole series is re-evaluated in double-double with direct
# log-gamma sums and no refreshes. `_logratio_jet_exact_extended` (ε = 0) and
# `_logratio_jet_near_extended` return `(status, value::ComplexF64, kappa)`.

const _LOGRATIO_JET_EXTENDED_MAX_KAPPA = 1e14
const _LOGRATIO_JET_EXTENDED_TOL = 2.0^-60
# A family may be κ times larger than the final value, so each family series is summed to full
# double-double precision. The pair series instead stops on its tail certificate.
const _DD_FAMILY_SERIES_TOL = 2.0^-100

_dd_lg(x) = _dd_loggamma(_CDD(x))
_dd_lg(x::_CDD) = _dd_loggamma(x)
_dd_psi(x::_CDD) = _dd_digamma(x)

# Ordinary family series C * w^c * pFq(α; β; ±w) in double-double, adding every
# term to the absolute sum so that cancellation inside the family is measured.
function _dd_family_series(logC::_CDD, logw::_CDD, c::_CDD, α, β, x::_CDD; maxterms=10_000)
    pref = _dd_exp(logC + c * logw)
    # A prefactor at a gamma pole (unresolved confluence) is non-finite in any precision.
    isfinite(abs(pref)) || return _CDD(_DD(NaN), _DD(NaN)), Inf, true
    t = one(_CDD)
    total = zero(_CDD)
    absum = 0.0
    small = 0
    for i in 0:(maxterms - 1)
        v = pref * t
        isfinite(abs(v)) || return _CDD(_DD(NaN), _DD(NaN)), Inf, true
        total = total + v
        absum += abs(v)
        if abs(v) <= _DD_FAMILY_SERIES_TOL * abs(total)
            small += 1
            small >= 2 && return total, absum, true
        else
            small = 0
        end
        ratio = x / (i + 1)
        for aj in α
            ratio = ratio * (aj + i)
        end
        for bj in β
            ratio = ratio / (bj + i)
        end
        t = t * ratio
    end
    return total, absum, false
end

function _dd_lower_family_term(A, B, m, nactive, Z::_CDD, logZ::_CDD, k::Integer; maxterms=10_000)
    p = length(A)
    q = length(B)
    c = B[k]
    logC = zero(_CDD)
    for j in 1:m
        j == k && continue
        logC = logC + _dd_lg(B[j] - c)
    end
    for j in 1:nactive
        logC = logC + _dd_lg(1 + c - A[j])
    end
    for j in (m + 1):q
        logC = logC - _dd_lg(1 + c - B[j])
    end
    for j in (nactive + 1):p
        logC = logC - _dd_lg(A[j] - c)
    end
    α = ntuple(j -> 1 + c - A[j], length(A))
    β = ntuple(idx -> (j = idx < k ? idx : idx + 1; 1 + c - B[j]), length(B) - 1)
    x = isodd(p - m - nactive) ? -Z : Z
    return _dd_family_series(logC, logZ, c, α, β, x; maxterms)
end

function _dd_upper_family_term(A, B, m, nactive, Z::_CDD, logZ::_CDD, h::Integer; maxterms=10_000)
    p = length(A)
    q = length(B)
    c = A[h]
    logC = zero(_CDD)
    for j in 1:nactive
        j == h && continue
        logC = logC + _dd_lg(c - A[j])
    end
    for j in 1:m
        logC = logC + _dd_lg(1 - c + B[j])
    end
    for j in (nactive + 1):p
        logC = logC - _dd_lg(1 - c + A[j])
    end
    for j in (m + 1):q
        logC = logC - _dd_lg(c - B[j])
    end
    α = ntuple(j -> 1 - c + B[j], length(B))
    β = ntuple(idx -> (j = idx < h ? idx : idx + 1; 1 - c + A[j]), length(A) - 1)
    invZ = inv(Z)
    x = isodd(q - m - nactive) ? -invZ : invZ
    return _dd_family_series(logC, logZ, c - 1, α, β, x; maxterms)
end

# Replace the second entry of a parameter tuple with the tuple type inferable (a splat of
# `t[3:end]` is not, and an `ntuple` closure over a reassigned variable boxes it).
_dd_with_second(t::Tuple, v) = ntuple(j -> j == 2 ? v : t[j], length(t))

# Unpaired residues of the first active family at levels ell < L.
function _dd_lower_unpaired_residue(A, B, m, nactive, logZ::_CDD, ell::Integer)
    p = length(A)
    q = length(B)
    s = B[1] + ell
    lv = -_dd_lg(ell + 1)
    for j in 2:m; lv = lv + _dd_lg(B[j] - s); end
    for j in 1:nactive; lv = lv + _dd_lg(1 + s - A[j]); end
    for j in (m + 1):q; lv = lv - _dd_lg(1 + s - B[j]); end
    for j in (nactive + 1):p; lv = lv - _dd_lg(A[j] - s); end
    return _dd_exp(lv + s * logZ) * (isodd(ell) ? -1.0 : 1.0)
end

function _dd_upper_unpaired_residue(A, B, m, nactive, logZ::_CDD, ell::Integer)
    p = length(A)
    q = length(B)
    s = A[1] - 1 - ell
    lv = -_dd_lg(ell + 1)
    for j in 1:m; lv = lv + _dd_lg(B[j] - s); end
    for j in 2:nactive; lv = lv + _dd_lg(1 - A[j] + s); end
    for j in (m + 1):q; lv = lv - _dd_lg(1 - B[j] + s); end
    for j in (nactive + 1):p; lv = lv - _dd_lg(A[j] - s); end
    return _dd_exp(lv + s * logZ) * (isodd(ell) ? -1.0 : 1.0)
end

# Unpaired residues and ordinary families shared by both lanes.  Returns
# `(total, absum, ok)`.
function _dd_ordinary_terms(A, B, m, nactive, Z::_CDD, logZ::_CDD, L::Integer, lower::Bool; maxterms)
    total = zero(_CDD)
    absum = 0.0
    for ell in 0:(L - 1)
        v = lower ? _dd_lower_unpaired_residue(A, B, m, nactive, logZ, ell) :
                    _dd_upper_unpaired_residue(A, B, m, nactive, logZ, ell)
        total = total + v
        absum += abs(v)
    end
    for k in 3:(lower ? m : nactive)
        v, va, ok = lower ? _dd_lower_family_term(A, B, m, nactive, Z, logZ, k; maxterms) :
                            _dd_upper_family_term(A, B, m, nactive, Z, logZ, k; maxterms)
        ok || return total, absum, false
        total = total + v
        absum += va
    end
    return total, absum, true
end

"""
    _logratio_jet_exact_extended(a, b, m, nactive, z, L, branch; maxterms)

Double-double evaluation of the exact-pair grouped series for canonical
`(a, b, L)` (as returned by `_detect_exact_logratio_jet_pair`).  Returns
`(status, value::ComplexF64, kappa)`.
"""
function _logratio_jet_exact_extended(a::Tuple, b::Tuple, m::Integer, nactive::Integer, z,
                                      L::Integer, branch::Symbol; maxterms::Integer=10_000)
    lower = branch === :lower
    A = map(_CDD, a)
    B = map(_CDD, b)
    # Canonical exact pair: b2 = b1 + L (lower) or a2 = a1 - L (upper).
    lower ? (B = _dd_with_second(B, B[1] + L)) : (A = _dd_with_second(A, A[1] - L))
    Z = _CDD(z)
    logZ = _dd_log(Z)
    p = length(A)
    q = length(B)
    total, absum, ok = _dd_ordinary_terms(A, B, m, nactive, Z, logZ, L, lower; maxterms)
    ok || return (status=:terms_exceeded, value=ComplexF64(NaN, NaN), kappa=Inf)

    # Start state of the grouped pair at n = 0 (K = L), direct sums.
    K = L
    if lower
        s = B[1] + K
        logX = -_dd_lg(K + 1) - _dd_lg(1)
        D = logZ - _dd_psi(_CDD(K + 1)) - _dd_psi(_CDD(1))
        for j in 3:m
            logX = logX + _dd_lg(B[j] - s); D = D - _dd_psi(B[j] - s)
        end
        for j in 1:nactive
            logX = logX + _dd_lg(1 - A[j] + s); D = D + _dd_psi(1 - A[j] + s)
        end
        for j in (m + 1):q
            logX = logX - _dd_lg(1 + s - B[j]); D = D - _dd_psi(1 + s - B[j])
        end
        for j in (nactive + 1):p
            logX = logX - _dd_lg(A[j] - s); D = D + _dd_psi(A[j] - s)
        end
        P = _dd_exp(logX + s * logZ) * (isodd(L + 1) ? -1.0 : 1.0)
    else
        s = A[1] - 1 - K
        logX = -_dd_lg(K + 1) - _dd_lg(1)
        D = logZ + _dd_psi(_CDD(K + 1)) + _dd_psi(_CDD(1))
        for j in 1:m
            logX = logX + _dd_lg(B[j] - s); D = D - _dd_psi(B[j] - s)
        end
        for j in 3:nactive
            logX = logX + _dd_lg(1 - A[j] + s); D = D + _dd_psi(1 - A[j] + s)
        end
        for j in (m + 1):q
            logX = logX - _dd_lg(1 - B[j] + s); D = D - _dd_psi(1 - B[j] + s)
        end
        for j in (nactive + 1):p
            logX = logX - _dd_lg(A[j] - s); D = D + _dd_psi(A[j] - s)
        end
        P = _dd_exp(logX + s * logZ) * (isodd(L) ? -1.0 : 1.0)
    end

    tail_constants = _logratio_tail_constants(a, b, L, zero(real(float(first(b)))), branch)
    invZ = inv(Z)
    for nn in 0:(maxterms - 1)
        K = L + nn
        if nn > 0
            if lower
                s = B[1] + K
                sold = s - 1
                ratio = Z / (K * nn)
                for j in 3:m
                    ratio = ratio / (B[j] - s); D = D + inv(B[j] - s)
                end
                for j in 1:nactive
                    ratio = ratio * (1 - A[j] + sold); D = D + inv(1 - A[j] + sold)
                end
                for j in (m + 1):q
                    ratio = ratio / (1 + sold - B[j]); D = D - inv(1 + sold - B[j])
                end
                for j in (nactive + 1):p
                    ratio = ratio * (A[j] - s); D = D - inv(A[j] - s)
                end
                D = D - inv(_CDD(nn)) - inv(_CDD(K))
            else
                s = A[1] - 1 - K
                sold = s + 1
                ratio = invZ / (K * nn)
                for j in 1:m
                    ratio = ratio * (B[j] - sold); D = D - inv(B[j] - sold)
                end
                for j in 3:nactive
                    ratio = ratio / (1 - A[j] + s); D = D - inv(1 - A[j] + s)
                end
                for j in (m + 1):q
                    ratio = ratio * (1 - B[j] + s); D = D + inv(1 - B[j] + s)
                end
                for j in (nactive + 1):p
                    ratio = ratio / (A[j] - sold); D = D + inv(A[j] - sold)
                end
                D = D + inv(_CDD(K)) + inv(_CDD(nn))
            end
            P = P * ratio
        end
        term = P * D
        total = total + term
        absum += abs(term)
        isfinite(absum) || return (status=:nonfinite, value=ComplexF64(NaN, NaN), kappa=Inf)
        cert = _logratio_pair_tail_bound(ComplexF64(P), ComplexF64(D), a, b, z, L,
                                         zero(real(float(first(b)))), nn, branch;
                                         constants=tail_constants)
        if cert.tail_bound <= _LOGRATIO_JET_EXTENDED_TOL * abs(total)
            value = ComplexF64(total)
            return (status=:ok, value=value, kappa=absum / abs(value))
        end
    end
    return (status=:terms_exceeded, value=ComplexF64(total), kappa=Inf)
end

# Signed divided difference (log Γ(x + σε) - log Γ(x)) / ε in double-double. The argument is
# shifted to real part >= `_DD_SHIFT` with scalar `log1p` terms, so the jet needs no shift loop.
function _dd_loggamma_divdiff(x::_CDD, E::_CDD, sigma::Integer)
    S = max(0, ceil(Int, _DD_SHIFT - Float64(real(x))))
    val = _logratio_loggamma_divdiff_signed(x + S, E, sigma)
    for k in 0:(S - 1)
        val = val - log1p(sigma * E / (x + k)) / E
    end
    return val
end

"""
    _logratio_jet_near_extended(a, b, m, nactive, z, L, epsv, branch; maxterms)

Double-double evaluation of the finite-separation grouped series for a
canonical near pair (`b[2] = b[1] + L + epsv` lower, `a[2] = a[1] - L + epsv`
upper, as returned by `_detect_near_logratio_jet_pair`).  Returns
`(status, value::ComplexF64, kappa)`.
"""
function _logratio_jet_near_extended(a::Tuple, b::Tuple, m::Integer, nactive::Integer, z,
                                     L::Integer, epsv, branch::Symbol; maxterms::Integer=10_000)
    lower = branch === :lower
    A = map(_CDD, a)
    B = map(_CDD, b)
    E = _CDD(epsv)
    invE = inv(E)
    # The pair is carried as the exact triple (first parameter, L, ε).
    lower ? (B = _dd_with_second(B, B[1] + L + E)) : (A = _dd_with_second(A, A[1] - L + E))
    Z = _CDD(z)
    logZ = _dd_log(Z)
    p = length(A)
    q = length(B)
    total, absum, ok = _dd_ordinary_terms(A, B, m, nactive, Z, logZ, L, lower; maxterms)
    ok || return (status=:terms_exceeded, value=ComplexF64(NaN, NaN), kappa=Inf)

    # Start state at n = 0 (K = L): mirrors `_logratio_lower/upper_pair_state`.
    K = L
    onec = one(_CDD)
    if lower
        s = B[1] + K
        logX = _dd_lg(onec + E) - _dd_lg(K + 1)
        D = _dd_loggamma_divdiff(onec, E, -1) - _dd_loggamma_divdiff(onec, E, +1) + logZ
        for j in 3:m
            logX = logX + _dd_lg(B[j] - s); D = D + _dd_loggamma_divdiff(B[j] - s, E, -1)
        end
        for j in 1:nactive
            logX = logX + _dd_lg(1 - A[j] + s); D = D + _dd_loggamma_divdiff(1 - A[j] + s, E, +1)
        end
        for j in (m + 1):q
            logX = logX - _dd_lg(1 + s - B[j]); D = D - _dd_loggamma_divdiff(1 + s - B[j], E, +1)
        end
        for j in (nactive + 1):p
            logX = logX - _dd_lg(A[j] - s); D = D - _dd_loggamma_divdiff(A[j] - s, E, -1)
        end
        D = D - _logratio_logpoch_divdiff_signed(onec, K, E, +1)
        P = _dd_exp(logX + s * logZ) * (isodd(L + 1) ? -1.0 : 1.0)
    else
        s = A[1] - 1 - K
        logX = _dd_lg(onec - E) - _dd_lg(K + 1)
        D = _dd_loggamma_divdiff(onec, E, +1) - _dd_loggamma_divdiff(onec, E, -1) + logZ
        for j in 1:m
            logX = logX + _dd_lg(B[j] - s); D = D + _dd_loggamma_divdiff(B[j] - s, E, -1)
        end
        for j in 3:nactive
            logX = logX + _dd_lg(1 - A[j] + s); D = D + _dd_loggamma_divdiff(1 - A[j] + s, E, +1)
        end
        for j in (m + 1):q
            logX = logX - _dd_lg(1 - B[j] + s); D = D - _dd_loggamma_divdiff(1 - B[j] + s, E, +1)
        end
        for j in (nactive + 1):p
            logX = logX - _dd_lg(A[j] - s); D = D - _dd_loggamma_divdiff(A[j] - s, E, -1)
        end
        D = D - _logratio_logpoch_divdiff_signed(onec, K, E, -1)
        P = _dd_exp(logX + s * logZ) * (isodd(L) ? -1.0 : 1.0)
    end

    tail_constants = _logratio_tail_constants(a, b, L, epsv, branch)
    invZ = inv(Z)
    for nn in 0:(maxterms - 1)
        K = L + nn
        if nn > 0
            if lower
                s = B[1] + K
                sold = s - 1
                ratio = Z / (K * (nn - E))
                for j in 3:m
                    ratio = ratio / (B[j] - s); D = D - log1p(-E / (B[j] - s)) * invE
                end
                for j in 1:nactive
                    ratio = ratio * (1 - A[j] + sold); D = D + log1p(E / (1 - A[j] + sold)) * invE
                end
                for j in (m + 1):q
                    ratio = ratio / (1 + sold - B[j]); D = D - log1p(E / (1 + sold - B[j])) * invE
                end
                for j in (nactive + 1):p
                    ratio = ratio * (A[j] - s); D = D + log1p(-E / (A[j] - s)) * invE
                end
                D = D + log1p(-E / nn) * invE - log1p(E / K) * invE
            else
                s = A[1] - 1 - K
                sold = s + 1
                ratio = invZ / (K * (nn + E))
                for j in 1:m
                    ratio = ratio * (B[j] - sold); D = D + log1p(-E / (B[j] - sold)) * invE
                end
                for j in 3:nactive
                    ratio = ratio / (1 - A[j] + s); D = D - log1p(E / (1 - A[j] + s)) * invE
                end
                for j in (m + 1):q
                    ratio = ratio * (1 - B[j] + s); D = D + log1p(E / (1 - B[j] + s)) * invE
                end
                for j in (nactive + 1):p
                    ratio = ratio / (A[j] - sold); D = D - log1p(-E / (A[j] - sold)) * invE
                end
                D = D - log1p(-E / K) * invE + log1p(E / nn) * invE
            end
            P = P * ratio
        end
        term = P * D * _logratio_exprel(E * D)
        total = total + term
        absum += abs(term)
        isfinite(absum) || return (status=:nonfinite, value=ComplexF64(NaN, NaN), kappa=Inf)
        cert = _logratio_pair_tail_bound(ComplexF64(P), ComplexF64(D), a, b, z, L, epsv, nn,
                                         branch; constants=tail_constants)
        if cert.tail_bound <= _LOGRATIO_JET_EXTENDED_TOL * abs(total)
            value = ComplexF64(total)
            return (status=:ok, value=value, kappa=absum / abs(value))
        end
    end
    return (status=:terms_exceeded, value=ComplexF64(total), kappa=Inf)
end
