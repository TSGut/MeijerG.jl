# Private higher-order exact Laurent residue evaluator.

const _HIGHER_LAURENT_MAXELL = 40
const _HIGHER_LAURENT_MAXMULT = 4

"""
    _HigherGammaFactor(offset, slope)

Represent a Mellin-Barnes gamma factor `Gamma(offset + slope*s)`.
"""
struct _HigherGammaFactor{T}
    offset::T
    slope::T
end

"""
    _higher_mb_factors(a, b, m, n)

Return numerator and denominator gamma factors of the Mellin-Barnes kernel.
"""
function _higher_mb_factors(a::Tuple, b::Tuple, m::Integer, n::Integer,
                            ::Type{T}=_cluster_parameter_type(a, b)) where {T}
    num = _HigherGammaFactor{T}[]
    den = _HigherGammaFactor{T}[]
    # The parameters are converted before the offsets are formed, so that in a wider
    # working type `1 - a[j]` is exact and the factor's pole coincides with the pole
    # position formed from the parameter.
    for j in 1:m
        push!(num, _HigherGammaFactor{T}(T(b[j]), -one(T)))
    end
    for j in 1:n
        push!(num, _HigherGammaFactor{T}(one(T) - T(a[j]), one(T)))
    end
    for j in (m + 1):length(b)
        push!(den, _HigherGammaFactor{T}(one(T) - T(b[j]), one(T)))
    end
    for j in (n + 1):length(a)
        push!(den, _HigherGammaFactor{T}(T(a[j]), -one(T)))
    end
    return (num=num, den=den)
end

"""
    _higher_pole_order(arg)

Return the nonnegative integer pole order of `Gamma(arg)`, otherwise `nothing`.
"""
_higher_pole_order(arg) = _cluster_pole_order(arg)

"""
    _higher_complex_type(z)

Return the native complex coefficient type associated with `z`.
"""
_higher_complex_type(z) = typeof(complex(float(real(z)), float(imag(z))))

"""
    _higher_hjet(r, N, T)

Return the jet of `u*Gamma(-r+u)` using
`u*Gamma(-r+u)=Gamma(1+u)/prod_{k=1}^r (u-k)`.
"""
function _higher_hjet(r::Integer, N::Integer, ::Type{T}) where {T<:Number}
    u = _jet_var(T, N)
    g = _jet_exp(_higher_loggamma_jet(one(T), one(T), N))
    denom = _jet_const(one(T), N)
    for k in 1:r
        denom = denom * (u - T(k))
    end
    return g * _jet_inv(denom)
end

"""
    _higher_singular_factor_jet(slope, r, N)

Return the jet of `t*Gamma(-r+slope*t)`.
"""
function _higher_singular_factor_jet(slope::T, r::Integer, N::Integer) where {T<:Number}
    h = _higher_hjet(r, N, T)
    return (one(T) / slope) * _jet_compose_scale(h, slope)
end

"""
    _higher_laurent_cluster_coefficients(num, den, s0, N)

Return negative Laurent coefficients at `s0`, or `nothing` for unsupported
local structures such as denominator poles.
"""
function _higher_laurent_cluster_coefficients(num, den, s0::T, N::Integer) where {T<:Number}
    logR = _jet_const(zero(T), N)
    Gprod = _jet_const(one(T), N)
    nsing = 0
    for f in num
        arg0 = f.offset + f.slope * s0
        po = _higher_pole_order(arg0)
        if po === nothing
            logR = logR + _higher_loggamma_jet(arg0, f.slope, N)
        else
            Gprod = Gprod * _higher_singular_factor_jet(f.slope, po, N)
            nsing += 1
        end
    end
    for f in den
        arg0 = f.offset + f.slope * s0
        _higher_pole_order(arg0) === nothing || return nothing
        logR = logR - _higher_loggamma_jet(arg0, f.slope, N)
    end
    nsing >= 1 || return nothing
    nsing <= N + 1 || return nothing
    G = Gprod * _jet_exp(logR)
    coeffs = Dict{Int,T}()
    for q in 0:(nsing - 1)
        coeffs[-1 - q] = _jet_coeff(G, nsing - 1 - q)
    end
    return (coeffs=coeffs, multiplicity=nsing, G=G)
end

"""
    _higher_laurent_residue_folded(num, den, s0, z, orientation, r)

Return one folded log-polynomial residue contribution with `z^s0` included in
the jet exponent.  This avoids overflow in negligible tail terms.
"""
function _higher_laurent_residue_folded(num, den, s0::T, z::T,
                                        orientation::Integer, r::Integer) where {T<:Complex}
    N = r - 1
    logz = log(z)
    logR = _jet_const(s0 * logz, N)
    Gprod = _jet_const(one(T), N)
    nsing = 0
    for f in num
        arg0 = f.offset + f.slope * s0
        po = _higher_pole_order(arg0)
        if po === nothing
            logR = logR + _higher_loggamma_jet(arg0, f.slope, N)
        else
            Gprod = Gprod * _higher_singular_factor_jet(f.slope, po, N)
            nsing += 1
        end
    end
    for f in den
        arg0 = f.offset + f.slope * s0
        _higher_pole_order(arg0) === nothing || return nothing
        logR = logR - _higher_loggamma_jet(arg0, f.slope, N)
    end
    nsing == r || return nothing
    Gz = Gprod * _jet_exp(logR)
    val = zero(T)
    for q in 0:(r - 1)
        val += orientation * _jet_coeff(Gz, r - 1 - q) * logz^q / factorial(q)
    end
    return val
end

"""
    _higher_laurent_residue_logpoly(coeffs, s0, z, orientation)

Evaluate the direct log-polynomial residue from negative Laurent coefficients.
"""
function _higher_laurent_residue_logpoly(coeffs::Dict{Int,T}, s0, z,
                                         orientation::Integer) where {T}
    r = length(coeffs)
    val = zero(z * one(T) + 0im)
    logz = log(z)
    for q in 0:(r - 1)
        val += orientation * coeffs[-1 - q] * z^s0 * logz^q / factorial(q)
    end
    return val
end

"""
    _higher_laurent_branch_sum(a, b, m, n, z, branch; maxell, maxmult)

Sum exact higher-order Laurent residue clusters for one branch.
"""
function _higher_laurent_branch_sum(a::Tuple, b::Tuple, m::Integer, n::Integer, z,
                                    branch::Symbol;
                                    maxell::Integer=_HIGHER_LAURENT_MAXELL,
                                    maxmult::Integer=_HIGHER_LAURENT_MAXMULT)
    C = _higher_complex_type(z)
    factors = _higher_mb_factors(a, b, m, n, C)
    num = factors.num
    den = factors.den
    pole_clusters = _detect_active_confluence_clusters(a, b, m, n, branch;
                                                       maxell, near_radius=0.0,
                                                       include_singletons=true)
    if branch === :upper
        sort!(pole_clusters; by=c -> real(c.center), rev=true)
    end
    orientation = branch === :lower ? -1 : 1
    total = zero(C)
    absum = zero(real(C))
    maxmult_seen = 0
    smallcount = 0
    converged = false
    R = _cluster_real_float_type(z)
    termtol = convert(R, 10) * eps(R)
    for cluster in pole_clusters
        cluster.kind === :exact || continue
        s0 = convert(C, cluster.center)
        r = count(f -> _higher_pole_order(f.offset + f.slope * s0) !== nothing, num)
        r >= 1 || continue
        maxmult_seen = max(maxmult_seen, r)
        r <= maxmult || return (status=:unsupported_multiplicity, value=missing,
                                max_multiplicity=r)
        contrib = _higher_laurent_residue_folded(num, den, s0, convert(C, z),
                                                 orientation, r)
        contrib === nothing &&
            return (status=:unsupported_local_structure, value=missing,
                    max_multiplicity=maxmult_seen)
        if !_logratio_finite(contrib)
            smallcount >= 1 && (converged = true; break)
            return (status=:nonfinite, value=missing, max_multiplicity=maxmult_seen)
        end
        total += contrib
        absum += abs(contrib)
        if abs(contrib) <= termtol * max(abs(total), floatmin(R))
            smallcount += 1
            smallcount >= 2 && (converged = true; break)
        else
            smallcount = 0
        end
    end
    maxmult_seen >= 2 || return (status=:not_near_confluent, value=missing,
                                 max_multiplicity=maxmult_seen)
    # Running out of enumerated pole levels before the tail test passes means the sum
    # was truncated (e.g. maxell too small for the working precision).
    converged || return (status=:not_converged, value=total, max_multiplicity=maxmult_seen,
                         absolute_sum=absum)
    return (status=:ok, value=total, max_multiplicity=maxmult_seen, absolute_sum=absum)
end

"""
    _try_meijerg_slater_higher_laurent_exact(a, b, m, n, z; maxell, maxmult)

Try the standalone exact higher-order Laurent evaluator.  Returns a structured
result and never falls back to another method.
"""
function _try_meijerg_slater_higher_laurent_exact(a::Tuple, b::Tuple,
        m::Integer, n::Integer, z; maxell::Integer=_HIGHER_LAURENT_MAXELL,
        maxmult::Integer=_HIGHER_LAURENT_MAXMULT, working_precision::Symbol=:auto)
    p = length(a)
    q = length(b)
    p == q && return (status=:unsupported_finite_p_eq_q, value=missing,
                      max_multiplicity=0)
    branch = p < q ? :lower : :upper
    summary = _active_confluence_summary(a, b, m, n, z; branch,
                                         maxell, near_radius=0.0)
    summary.status === :ok || return (status=summary.status, value=missing,
                                      max_multiplicity=summary.max_multiplicity)
    only(summary.clusters).kind === :exact ||
        return (status=:not_near_confluent, value=missing,
                max_multiplicity=summary.max_multiplicity)
    # Enumerate more pole levels until the tail test passes (higher precision needs more)
    function sum_converged(aa, bb, zz)
        ell = maxell
        while true
            result = _higher_laurent_branch_sum(aa, bb, m, n, zz, branch; maxell=ell, maxmult)
            (result.status === :not_converged && ell < 8maxell) || return result
            ell *= 2
        end
    end
    result = sum_converged(a, b, z)
    result.status === :ok || return result
    # Cancellation between residue terms: for BigFloat, rerun at a precision raised
    # by log2(kappa). Other types keep the value (the caller decides).
    kappa = _logratio_summation_kappa(result.absolute_sum, result.value)
    _logratio_summation_ok(kappa) && return result
    hi = working_precision !== :input ?
         _logratio_escalated(sum_converged, (a, b, z), z, kappa) : nothing
    hi === nothing && return merge(result, (status=:summation_cancellation,))
    return hi
end
