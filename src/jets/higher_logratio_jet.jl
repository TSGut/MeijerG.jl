# Private higher-order exact and near log-ratio jet evaluator.

const _HIGHER_LOGRATIO_MAXMULT = 4
const _HIGHER_LOGRATIO_NEAR_RADIUS = 1.1e-2
const _HIGHER_LOGRATIO_MAXELL = 40
const _HIGHER_LOGRATIO_TOL = 1e-14
const _HIGHER_LOGRATIO_MAX_GROUPED_TERMS = 80
const _HIGHER_LOGRATIO_LOCAL_PAD = 16
const _HIGHER_LOGRATIO_LOCAL_PAD_MAX = 48
const _HIGHER_LOGRATIO_DIRECT_RECIPGAMMA_RADIUS = 1e-1
const _HIGHER_LOGRATIO_EXTENDED_MAXELL = 400

"""
    _higher_logratio_maxell(z)

Return a convergence-adaptive cap on the pole index enumerated by the active
cluster detector.  Significant clusters of the lower/upper Slater expansion decay
like `ρ^ℓ` with `ρ = min(|z|, 1/|z|)`, so we only need to enumerate enough levels
for `ρ^maxell` to fall well below machine precision. The summation's tail break
still governs termination.  Replaces a flat scan and is the dominant speedup for
isolated higher-order clusters, while `_HIGHER_LOGRATIO_MAXELL` remains the cap
near `|z| = 1` where convergence is slow.
"""
function _higher_logratio_maxell(z)
    ρ = min(abs(z), inv(abs(z)))
    ρ >= 0.95 && return _HIGHER_LOGRATIO_MAXELL
    return clamp(ceil(Int, -24 / log10(ρ)) + 2,
                 2 * _HIGHER_LOGRATIO_MAXMULT, _HIGHER_LOGRATIO_MAXELL)
end

"""
    _higher_logratio_maxell_extended(z, ν, branch; decades=32)

Pole-level cap for the double-double lane, valid on and off the convergent
regime: the residues of the branch series behave like `ρ^ℓ / (ℓ!)^ν` with
`ρ = |z|` (lower) or `1/|z|` (upper) and `ν = |q - p|`, and levels are enumerated
until that envelope has fallen `decades` decades below its peak.
"""
function _higher_logratio_maxell_extended(z, ν::Integer, branch::Symbol; decades::Real=32)
    logρ = branch === :lower ? log(abs(z)) : -log(abs(z))
    envelope(ℓ) = ℓ * logρ - ν * loggamma(ℓ + 1.0)
    best = envelope(0)
    ℓ = 0
    while ℓ < _HIGHER_LOGRATIO_EXTENDED_MAXELL
        ℓ += 1
        v = envelope(ℓ)
        best = max(best, v)
        best - v >= decades * log(10) && break
    end
    return clamp(ℓ + 4, 2 * _HIGHER_LOGRATIO_MAXMULT, _HIGHER_LOGRATIO_EXTENDED_MAXELL)
end

"""
    _higher_logratio_maxell(z, ν, branch)

Pole-level cap for the Float64 jet in the convergent regime: the envelope rule
of `_higher_logratio_maxell_extended` with 22 decades, capped at
`_HIGHER_LOGRATIO_MAXELL`.  A series that runs out of levels is re-run at the cap.
"""
_higher_logratio_maxell(z, ν::Integer, branch::Symbol) =
    min(_higher_logratio_maxell_extended(z, ν, branch; decades=22), _HIGHER_LOGRATIO_MAXELL)

# Largest integer offset between active families: the levels before every family has a
# pole in the cluster.
function _higher_logratio_span(active)
    span = 0
    for i in eachindex(active), j in (i + 1):length(active)
        d = real(active[j] - active[i])
        isfinite(d) && (span = max(span, abs(round(Int, d))))
    end
    return span
end

# In-place jet accumulation with the coefficient operations of `+(::Jet, ::Jet)`
# and `-(::Jet, ::Jet)` (`x + y`, `x + (-y)`), so results are unchanged.
function _jet_add!(acc, jet, sign::Integer)
    @inbounds for k in eachindex(acc.coeffs)
        acc.coeffs[k] = sign > 0 ? acc.coeffs[k] + jet.coeffs[k] : acc.coeffs[k] + (-jet.coeffs[k])
    end
    return acc
end
_jet_add_const!(acc, s) = (acc.coeffs[1] += s; acc)

const _HIGHER_LOGRATIO_ANCHOR_INTERVAL = 16
const _HIGHER_LOGRATIO_EXTENDED_ANCHOR_INTERVAL = 16

_higher_logratio_anchor_interval(::Type{C}) where {C} =
    C === _CDD ? _HIGHER_LOGRATIO_EXTENDED_ANCHOR_INTERVAL : _HIGHER_LOGRATIO_ANCHOR_INTERVAL

"""
    _HigherClusterState

State of one evaluated cluster from which the next cluster of the same shape
(same families, every `ell` one larger) is reached by one log jet per gamma
factor: `prefactor` is the exponential of the constant term and `Λ` the log
jet of `g(s0 + t)` with constant term zero, so the recurrence multiplies and
divides the prefactor instead of adding to a logarithm of large magnitude.
"""
struct _HigherClusterState{C}
    families::Vector{Int}
    ells::Vector{Int}
    s0::C
    deltas::Vector{C}
    prefactor::C
    Λ::Jet{C}
    direct::Union{Nothing,Jet{C}}
    steps::Int
end

# Orders `1:N` of `log(c + slope*t)` (the coefficients of `affine_log_jet`) added to `Λ` with
# the given sign. The order-0 term belongs to the prefactor.
function _higher_log_affine_rest!(Λ, c, slope, sign::Integer)
    ratio = -slope / c
    power = one(ratio)
    @inbounds for p in 1:(length(Λ.coeffs) - 1)
        power *= ratio
        term = -power / p
        Λ.coeffs[p + 1] = sign > 0 ? Λ.coeffs[p + 1] + term : Λ.coeffs[p + 1] - term
    end
    return Λ
end

# The decision rule of `reciprocal_gamma`, so the recurrence knows whether a denominator
# factor is carried as a direct reciprocal-gamma jet.
function _higher_direct_den(arg0)
    pole_order = max(0, round(Int, -real(arg0)))
    return abs(arg0 + pole_order) <= _HIGHER_LOGRATIO_DIRECT_RECIPGAMMA_RADIUS
end

function _higher_logratio_state_value(state::_HigherClusterState, orientation::Integer)
    dd = state.direct === nothing ? exp_divided_difference(state.Λ, state.deltas) :
         exp_product_divided_difference(state.direct, state.Λ, state.deltas)
    return orientation * (state.prefactor * dd)
end

# The log jets of the non-member gamma factors are Taylor series in the cluster variable
# whose radius of convergence is the distance `c` from the expansion point to the nearest
# gamma singularity, so their truncation at order `N` leaves a relative error of about
# `(ρ/c)^N` with `ρ` the cluster spread. This returns the pad that brings that error
# below the unit roundoff of the coefficient type over the `levels` the state can be
# advanced through. A non-member active family (same slope sign as the members) that
# comes within the detector's radius joins the cluster at that level, which starts a
# fresh state, so its scan stops there. Direct reciprocal-gamma factors are entire.
function _higher_logratio_singularity_pad(num, den, members, s0::C, deltas, shift::Integer,
                                          levels::Integer) where {C}
    ρ = Float64(maximum(abs, deltas; init=zero(real(C))))
    ρ > 0 || return 0
    member_factors = Set(member.factor_index for member in members)
    active_slope = sign(Float64(real(num[first(members).factor_index].slope)))
    cmin = Inf
    for (idx, f) in enumerate(num)
        idx in member_factors && continue
        merges = sign(Float64(real(f.slope))) == active_slope
        cmin = min(cmin, _higher_chain_pole_distance(f.offset + f.slope * s0, f.slope, shift,
                                                     levels, merges))
    end
    for f in den
        arg0 = f.offset + f.slope * s0
        _higher_direct_den(arg0) && continue
        cmin = min(cmin, _higher_chain_pole_distance(arg0, f.slope, shift, levels, false))
    end
    cmin > ρ || return _HIGHER_LOGRATIO_LOCAL_PAD_MAX
    isfinite(cmin) || return 0
    needed = ceil(Int, log(Float64(eps(real(C)))) / log(ρ / cmin))
    return clamp(needed, 0, _HIGHER_LOGRATIO_LOCAL_PAD_MAX)
end

# Smallest distance from `arg0 + j*slope*shift`, `0 <= j < levels`, to a pole of Γ (the
# nonpositive integers). A merging factor's scan stops at the first level where it is
# within the detector's radius.
function _higher_chain_pole_distance(arg0, slope, shift::Integer, levels::Integer,
                                     merges::Bool)
    x0 = Float64(real(arg0))
    y = Float64(imag(arg0))
    step = Float64(real(slope)) * shift
    d = Inf
    for j in 0:(levels - 1)
        x = x0 + j * step
        dj = hypot(x - min(round(x), 0.0), y)
        merges && dj < _HIGHER_LOGRATIO_NEAR_RADIUS && break
        d = min(d, dj)
    end
    return d
end

"""
    _higher_logratio_anchor(num, den, members, z, orientation)

Evaluate one exact or near active pole cluster from scratch as a stable
divided difference of `g(s)=P(s)z^s`, returning its contribution and the
`_HigherClusterState` for the recurrence.  The coefficient type is that of `z`
(`ComplexF64` on the Float64 path, `Complex{_DD}` in the double-double lane).
"""
function _higher_logratio_anchor(num, den, members, z::C, orientation::Integer;
                                 levels::Integer=_HIGHER_LOGRATIO_ANCHOR_INTERVAL) where {C}
    r = length(members)
    positions = C[member.position for member in members]
    s0 = sum(positions) / r
    deltas = positions .- s0
    member_factors = Set(member.factor_index for member in members)
    # Extra local terms keep near-pole gamma factors from under-resolving the Opitz
    # jet. The order is adaptive in the cluster spread and, through
    # `_higher_logratio_singularity_pad`, in the distance to the nearest non-member
    # singularity relative to the working precision.
    N = r - 1 + max(_higher_logratio_pad(deltas),
                    _higher_logratio_singularity_pad(num, den, members, s0, deltas,
                                                     -orientation, levels))
    logz = log(z)
    Λ = _jet_const(s0 * logz, N) +
        _Jet(C[k == 1 ? logz : zero(C) for k in 0:N])
    direct = nothing

    # Λ is accumulated in place (coefficientwise `x + y` / `x + (-y)`, the operations
    # of the jet `+` and `-`), avoiding a temporary jet per factor.
    for member in members
        f = num[member.factor_index]
        slope = f.slope
        δk = member.position - s0
        u0 = -slope * δk
        _jet_add_const!(Λ, -log(slope))
        _jet_add!(Λ, _higher_loggamma_jet(1 + u0, slope, N), 1)
        for i in 1:member.ell
            _jet_add!(Λ, _higher_loglinear_jet(u0 - i, slope, N), -1)
        end
    end
    for (idx, f) in enumerate(num)
        idx in member_factors && continue
        _jet_add!(Λ, _higher_loggamma_jet(f.offset + f.slope * s0, f.slope, N), 1)
    end
    for f in den
        arg0 = f.offset + f.slope * s0
        reciprocal = _higher_direct_recipgamma_jet(arg0, f.slope, N)
        if reciprocal === nothing
            _jet_add!(Λ, _higher_loggamma_jet(arg0, f.slope, N), -1)
        else
            direct = direct === nothing ? reciprocal : direct * reciprocal
        end
    end
    prefactor = exp(Λ.coeffs[1])
    Λ.coeffs[1] = zero(C)
    state = _HigherClusterState{C}(Int[member.family for member in members],
                                   Int[member.ell for member in members],
                                   s0, deltas, prefactor, Λ, direct, 0)
    return _higher_logratio_state_value(state, orientation), state
end

"""
    _higher_logratio_advance(state, num, den, members, z, orientation, shift)

Reach the cluster `members`, the translate of `state` by `shift` (`+1` lower,
`-1` upper), from the previous state: `z^s` gains a factor `z^shift`, each
member's pole product gains one linear factor, and every other gamma factor
moves by one through `Γ(x+1) = xΓ(x)`.  Returns `nothing` when `members` is
not that translate or when a denominator factor changes between the log and the
direct reciprocal form, in which case the caller re-anchors.
"""
function _higher_logratio_advance(state::_HigherClusterState{C}, num, den, members,
                                  z::C, orientation::Integer, shift::Integer) where {C}
    length(members) == length(state.families) || return nothing
    for (k, member) in enumerate(members)
        (member.family == state.families[k] && member.ell == state.ells[k] + 1) ||
            return nothing
    end
    s0_old = state.s0
    s0 = s0_old + shift
    prefactor = shift > 0 ? state.prefactor * z : state.prefactor / z
    Λ = state.Λ
    N = length(Λ.coeffs) - 1
    direct = state.direct
    for f in den
        arg_old = f.offset + f.slope * s0_old
        arg_new = arg_old + f.slope * shift
        up = real(f.slope) * shift > 0
        d_old = _higher_direct_den(arg_old)
        d_old == _higher_direct_den(arg_new) || return nothing
        if d_old
            # `1/(x + slope t)` near a pole of Γ cancels against the direct factor.
            up && return nothing
            direct = direct * affine_jet(arg_new, f.slope, N)
        elseif up
            prefactor /= arg_old
            _higher_log_affine_rest!(Λ, arg_old, f.slope, -1)
        else
            prefactor *= arg_new
            _higher_log_affine_rest!(Λ, arg_new, f.slope, 1)
        end
    end
    for (k, member) in enumerate(members)
        f = num[member.factor_index]
        c = -f.slope * state.deltas[k] - member.ell
        prefactor /= c
        _higher_log_affine_rest!(Λ, c, f.slope, -1)
    end
    for (idx, f) in enumerate(num)
        any(member.factor_index == idx for member in members) && continue
        arg_old = f.offset + f.slope * s0_old
        if real(f.slope) * shift > 0
            prefactor *= arg_old
            _higher_log_affine_rest!(Λ, arg_old, f.slope, 1)
        else
            arg_new = arg_old + f.slope * shift
            prefactor /= arg_new
            _higher_log_affine_rest!(Λ, arg_new, f.slope, -1)
        end
    end
    next = _HigherClusterState{C}(state.families, Int[member.ell for member in members],
                                  s0, state.deltas, prefactor, Λ, direct, state.steps + 1)
    return _higher_logratio_state_value(next, orientation), next
end

"""
    _higher_logratio_factor_index(member, branch, m)

Map an active branch-family index to a numerator gamma-factor index.
"""
_higher_logratio_factor_index(member, branch::Symbol, m::Integer) =
    branch === :lower ? member.family : m + member.family

"""
    _higher_logratio_members(cluster, branch, m, active, C)

Attach numerator factor indices to detector cluster members, with the pole
positions recomputed from the active parameters in the coefficient type `C`
(`b[j] + ell` lower, `a[j] - 1 - ell` upper).  In the double-double lane this
keeps the poles exact. In ComplexF64 it reproduces the detector's positions.
"""
function _higher_logratio_members(cluster, branch::Symbol, m::Integer, active, ::Type{C}) where {C}
    members = Vector{_HigherMember{C}}(undef, length(cluster.members))
    for (k, member) in enumerate(cluster.members)
        position = branch === :lower ? C(active[member.family]) + member.ell :
                                       C(active[member.family]) - 1 - member.ell
        members[k] = _HigherMember{C}(member.family, member.ell, position,
                                      _higher_logratio_factor_index(member, branch, m))
    end
    return members
end

"""
    _HigherMember{C}

A pole of a cluster as the evaluators see it: its family and level, its position in
the coefficient type and the index of its numerator gamma factor.
"""
struct _HigherMember{C}
    family::Int
    ell::Int
    position::C
    factor_index::Int
end

# Evaluate one cluster. Clusters of different family sets interleave along a branch (a
# triple and a singleton alternate, say), so one state is carried per family set,
# re-anchored after `interval` recurrence steps or whenever the next cluster of that set
# is not the translate of the previous one.
function _higher_logratio_cluster!(states, cluster, branch, m, active, ::Type{C}, num, den,
                                   zC, orientation, shift, interval) where {C}
    members = _higher_logratio_members(cluster, branch, m, active, C)
    k = findfirst(st -> length(st.families) == length(members) &&
                        all(st.families[j] == members[j].family for j in eachindex(members)),
                  states)
    advanced = (k === nothing || states[k].steps >= interval) ? nothing :
               _higher_logratio_advance(states[k], num, den, members, zC, orientation, shift)
    value, state = advanced === nothing ?
                   _higher_logratio_anchor(num, den, members, zC, orientation; levels=interval) :
                   advanced
    k === nothing ? push!(states, state) : (states[k] = state)
    return value
end

"""
    _higher_logratio_branch_sum(a, b, m, n, z, branch; maxell, maxmult, near_radius, termtol, C, maxterms)

Evaluate the higher-order exact/near log-ratio jet series for one branch in the
coefficient type `C` (`ComplexF64`, or `Complex{_DD}` for the double-double
lane).  Pole levels are enumerated beyond `maxell` by the largest integer offset
between active families, so that every cluster summed is complete: an incomplete
cluster (a colliding family whose pole at that level was not enumerated) would
leave a non-member gamma factor at a pole.
"""
function _higher_logratio_branch_sum(a::Tuple, b::Tuple, m::Integer, n::Integer, z,
        branch::Symbol; maxell::Integer=_higher_logratio_maxell(z),
        maxmult::Integer=_HIGHER_LOGRATIO_MAXMULT,
        near_radius=_HIGHER_LOGRATIO_NEAR_RADIUS,
        termtol=_HIGHER_LOGRATIO_TOL, C::Type=ComplexF64,
        maxterms::Integer=_HIGHER_LOGRATIO_MAX_GROUPED_TERMS)
    factors = _higher_mb_factors(a, b, m, n, C)
    num = factors.num
    den = factors.den
    active = C[(branch === :lower ? b[1:m] : a[1:n])...]
    maxell_detect = maxell + _higher_logratio_span(active)
    clusters = _detect_active_confluence_clusters(a, b, m, n, branch;
                                                  maxell=maxell_detect, near_radius,
                                                  include_singletons=true)
    if branch === :lower
        base = minimum(real, active)
        filter!(c -> round(Int, real(c.center) - base) <= maxell_detect, clusters)
    else
        base = maximum(real, active)
        filter!(c -> round(Int, base - 1 - real(c.center)) <= maxell_detect, clusters)
        sort!(clusters; by=c -> real(c.center), rev=true)
    end
    orientation = branch === :lower ? -1 : 1
    total = zero(C)
    absum = 0.0
    maxmult_seen = 0
    terms = 0
    smallcount = 0
    converged = false
    # Consecutive clusters of a family set are translates of one another, so they are
    # anchored once per `interval` levels and advanced by the recurrence.
    interval = _higher_logratio_anchor_interval(C)
    shift = branch === :lower ? 1 : -1
    states = _HigherClusterState{C}[]
    zC = C(z)
    for cluster in clusters
        r = cluster.multiplicity
        maxmult_seen = max(maxmult_seen, r)
        r <= maxmult || return (status=:unsupported_multiplicity, value=missing,
                                max_multiplicity=r, terms=terms)
        contrib = try
            _higher_logratio_cluster!(states, cluster, branch, m, active, C, num, den, zC,
                                      orientation, shift, interval)
        catch err
            return (status=:internal_error, value=missing,
                    max_multiplicity=maxmult_seen, terms=terms,
                    reason=sprint(showerror, err))
        end
        terms += 1
        if !_logratio_finite(contrib)
            smallcount >= 1 && (converged = true; break)
            return (status=:nonfinite, value=missing,
                    max_multiplicity=maxmult_seen, terms=terms)
        end
        total += contrib
        absum += abs(contrib)
        if abs(contrib) <= termtol * max(abs(total), floatmin(Float64))
            smallcount += 1
            smallcount >= 2 && (converged = true; break)
        else
            smallcount = 0
        end
        terms <= maxterms ||
            return (status=:terms_exceeded, value=total,
                    max_multiplicity=maxmult_seen, terms=terms)
    end
    # Running out of complete clusters before the tail test passes means the sum
    # was truncated. Report it instead of returning the partial sum.
    converged || return (status=:terms_exceeded, value=total,
                         max_multiplicity=maxmult_seen, terms=terms)
    return (status=:ok, value=total, max_multiplicity=maxmult_seen, terms=terms,
            absolute_sum=absum)
end

"""
    _higher_logratio_policy(summary)

Apply the first guarded support table for higher-order active clusters.
"""
function _higher_logratio_policy(summary)
    summary.status === :ok || return (selected=false, status=summary.status,
                                      reason=summary.status)
    cluster = only(summary.clusters)
    r = cluster.multiplicity
    r in (3, 4) || return (selected=false, status=:unsupported_multiplicity,
                           reason=:unsupported_multiplicity)
    cluster.denominator_pole &&
        return (selected=false, status=:unsupported_local_structure,
                reason=:unsupported_local_structure)
    return (selected=true, status=:ok, reason=:selected)
end

"""
    _try_meijerg_slater_higher_logratio_jet(a, b, m, n, z; branch=:auto)

Try the private higher-order log-ratio jet for isolated `r=3` or `r=4` active
clusters.  Unsupported cases return non-selected statuses so `meijerg` can
fall through to existing fallback methods.

The summation-conditioning policy follows the pair jets: with
`working_precision=:auto` a Float64 series off the convergent regime, or one
whose `κ_sum` exceeds `_LOGRATIO_JET_MAX_SUMMATION_KAPPA`, is evaluated in the
double-double lane (`_higher_logratio_branch_sum` with `C = Complex{_DD}`), and
Float32 inputs use the Float64 path. `:extended` always takes the double-double lane
for Float64 inputs.
"""
function _try_meijerg_slater_higher_logratio_jet(a::Tuple, b::Tuple,
        m::Integer, n::Integer, z; branch::Symbol=:auto,
        maxell::Union{Nothing,Integer}=nothing,
        return_diagnostics::Bool=false, working_precision::Symbol=:auto)
    extended = working_precision !== :input
    !_logratio_jet_enabled() &&
        return _logratio_result(:not_near_confluent, _logratio_nan_value(z),
                                :unknown, false, :disabled, (;))
    try
        p = length(a)
        q = length(b)
        p == q && return _logratio_result(:unsupported_finite_p_eq_q,
                                          _logratio_nan_value(z), :unknown,
                                          false, :unsupported_finite_p_eq_q, (;))
        # ComplexF64 with Float64-tuned tolerances and pads, so only Float32/Float64 inputs
        # are selected. TODO: BigFloat support (generic element type, precision-dependent
        # termtol/pad/radii, validation against mpmath).
        R = _logratio_real_type(z, a, b)
        R <: Union{Float32,Float64} ||
            return _logratio_result(:unsupported_precision, _logratio_nan_value(z),
                                    :unknown, false, :unsupported_precision, (;))
        br = branch === :auto ? (p < q ? :lower : :upper) : branch
        levels = maxell === nothing ? _higher_logratio_maxell(z, abs(q - p), br) : maxell
        span = _higher_logratio_span(br === :lower ? b[1:m] : a[1:n])
        summary = _active_confluence_summary(a, b, m, n, z; branch=br,
                                             maxell=levels + span,
                                             near_radius=_HIGHER_LOGRATIO_NEAR_RADIUS)
        policy = _higher_logratio_policy(summary)
        policy.selected || return _logratio_result(policy.status, _logratio_nan_value(z),
                                                   br, false, policy.reason,
                                                   return_diagnostics ? (summary=summary,) : (;))
        if extended && R === Float32
            r = _try_meijerg_slater_higher_logratio_jet(map(_logratio_widen, a),
                    map(_logratio_widen, b), m, n, _logratio_widen(z);
                    branch=br, return_diagnostics, working_precision=:auto)
            return merge(r, (value=convert(Complex{Float32}, r.value),))
        end
        lane = extended && R === Float64
        function extended_lane()
            ext = _higher_logratio_branch_sum(a, b, m, n, z, br;
                      maxell=_higher_logratio_maxell_extended(z, abs(q - p), br),
                      near_radius=_HIGHER_LOGRATIO_NEAR_RADIUS,
                      termtol=_LOGRATIO_JET_EXTENDED_TOL, C=_CDD,
                      maxterms=_HIGHER_LOGRATIO_EXTENDED_MAXELL + 8)
            value = ext.status === :ok ? ComplexF64(ext.value) : _logratio_nan_value(z)
            kappa = ext.status === :ok ? ext.absolute_sum / abs(value) : Inf
            diagnostics = return_diagnostics ?
                (status=ext.status, value=value, lane=:extended, summation_kappa=kappa,
                 terms=ext.terms, max_multiplicity=ext.max_multiplicity, summary=summary) : (;)
            ext.status === :ok ||
                return _logratio_result(ext.status, value, br, false, ext.status, diagnostics)
            kappa <= _LOGRATIO_JET_EXTENDED_MAX_KAPPA ||
                return _logratio_result(:summation_cancellation, value, br, false,
                                        :summation_cancellation, diagnostics)
            _logratio_finite(value) ||
                return _logratio_result(:nonfinite, value, br, false, :nonfinite, diagnostics)
            return _logratio_result(:ok, value, br, true, :selected, diagnostics)
        end
        lane && (working_precision === :extended || !_logratio_in_regime(br, z)) &&
            return extended_lane()
        diag = _higher_logratio_branch_sum(a, b, m, n, z, br;
                                           maxell=levels, near_radius=_HIGHER_LOGRATIO_NEAR_RADIUS)
        if diag.status === :terms_exceeded && levels < _HIGHER_LOGRATIO_MAXELL
            diag = _higher_logratio_branch_sum(a, b, m, n, z, br;
                                               maxell=_HIGHER_LOGRATIO_MAXELL,
                                               near_radius=_HIGHER_LOGRATIO_NEAR_RADIUS)
        end
        diagnostics = return_diagnostics ? merge(diag, (summary=summary,)) : (;)
        diag.status === :ok ||
            return _logratio_result(diag.status, _logratio_nan_value(z), br,
                                    false, diag.status, diagnostics)
        _logratio_finite(diag.value) ||
            return _logratio_result(:nonfinite, diag.value, br, false,
                                    :nonfinite, diagnostics)
        # Summation-conditioning certificate shared with the pair jet.
        kappa = _logratio_summation_kappa(diag.absolute_sum, diag.value)
        if !_logratio_summation_ok(kappa)
            lane && return extended_lane()
            return _logratio_result(:summation_cancellation, diag.value, br, false,
                                    :summation_cancellation,
                                    return_diagnostics ? merge(diagnostics, (summation_kappa=kappa,)) : (;))
        end
        return _logratio_result(:ok, diag.value, br, true, :selected,
                                return_diagnostics ? merge(diagnostics, (summation_kappa=kappa,)) : (;))
    catch e
        return _logratio_result(:internal_error, _logratio_nan_value(z),
                                branch === :auto ? :unknown : branch,
                                false, :internal_error, (error=sprint(showerror, e),))
    end
end
