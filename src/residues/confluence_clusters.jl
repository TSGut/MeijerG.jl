# Private active-pole confluence cluster detection. Geometry-only helpers for exact and near
# clusters of any multiplicity, without deciding which multiplicities are numerically supported.

const _CONFLUENCE_CLUSTER_MAXELL = 40
const _CONFLUENCE_CLUSTER_NEAR_RADIUS = 1.1e-2

"""
    _cluster_real_float_type(x)

Return the floating real type associated with `x`.
"""
_cluster_real_float_type(x) = typeof(float(real(x)))

"""
    _cluster_exact_tol(x)

Return the scale-aware tolerance used to classify exact pole coincidences.
"""
function _cluster_exact_tol(x)
    R = _cluster_real_float_type(x)
    R === BigFloat && return sqrt(eps(R)) * max(one(R), abs(x))
    return convert(R, 1e-14) * max(one(R), abs(x))
end

"""
    _cluster_position(branch, value, ell)

Return the Mellin-Barnes pole position for an active lower or upper family.
"""
function _cluster_position(branch::Symbol, value, ell::Integer)
    branch === :lower && return value + ell
    branch === :upper && return value - one(value) - ell
    throw(ArgumentError("unsupported active-pole branch $branch"))
end

"""
    _cluster_active_values(a, b, m, n, branch)

Return the active parameters for the selected branch.
"""
function _cluster_active_values(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                branch::Symbol)
    branch === :lower && return b[1:m]
    branch === :upper && return a[1:n]
    throw(ArgumentError("unsupported active-pole branch $branch"))
end

"""
    _cluster_denominator_factors(a, b, m, n)

Return denominator gamma factors of the Mellin-Barnes kernel as
`(offset, slope)` pairs representing `Gamma(offset + slope*s)`.
"""
function _cluster_denominator_factors(a::Tuple, b::Tuple, m::Integer, n::Integer)
    T = _cluster_parameter_type(a, b)
    den = Vector{NamedTuple{(:offset, :slope),Tuple{T,T}}}(undef, 0)
    for j in (m + 1):length(b)
        push!(den, (offset=one(b[j]) - b[j], slope=one(b[j])))
    end
    for j in (n + 1):length(a)
        push!(den, (offset=a[j], slope=-one(a[j])))
    end
    return den
end

"""
    _cluster_parameter_type(a, b)

Return the common type of the parameters, `Float64` when there are none.
"""
_cluster_parameter_type(a::Tuple, b::Tuple) =
    isempty(a) && isempty(b) ? Float64 : promote_type(map(typeof, (a..., b...))...)

"""
    _ClusterCandidate{P}

One active pole: its family index, its level `ell` and its position.
"""
const _ClusterCandidate{P} = NamedTuple{(:family, :ell, :position),Tuple{Int,Int,P}}

"""
    _ActiveCluster{P,R}

A group of active poles at a common position (`kind` is `:exact` or `:near`), as
returned by `_detect_active_confluence_clusters`.  The record is concretely typed
for every multiplicity, so the cluster vector and the evaluators that walk it
avoid dynamic dispatch.
"""
struct _ActiveCluster{P,R}
    kind::Symbol
    branch::Symbol
    multiplicity::Int
    center::P
    diameter::R
    members::Vector{_ClusterCandidate{P}}
    families::Vector{Int}
    ells::Vector{Int}
    shifts::Vector{Int}
    denominator_pole::Bool
end

"""
    _cluster_pole_order(arg)

Return the nonnegative integer pole order of `Gamma(arg)` when `arg` is an
integer pole, otherwise `nothing`.
"""
function _cluster_pole_order(arg)
    tol = _cluster_exact_tol(arg)
    abs(imag(arg)) <= tol || return nothing
    r = round(Int, -real(arg))
    (r >= 0 && abs(arg + r) <= tol) && return r
    return nothing
end

"""
    _cluster_has_denominator_pole(a, b, m, n, center)

Return whether a denominator gamma factor has a pole at `center`.
"""
function _cluster_has_denominator_pole(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                       center)
    return _cluster_has_denominator_pole(_cluster_denominator_factors(a, b, m, n), center)
end

_cluster_has_denominator_pole(den::AbstractVector, center) =
    any(f -> _cluster_pole_order(f.offset + f.slope * center) !== nothing, den)

"""
    _cluster_candidates(a, b, m, n, branch; maxell)

Return all active pole candidates up to `maxell` for the selected branch.
"""
function _cluster_candidates(a::Tuple, b::Tuple, m::Integer, n::Integer,
                             branch::Symbol; maxell::Integer)
    # Concretely typed values and candidates: the greedy clustering below is quadratic in
    # their number, and field access on an `Any` vector, or indexing a parameter tuple of
    # runtime length, dominated its cost.
    T = _cluster_parameter_type(a, b)
    active = T[_cluster_active_values(a, b, m, n, branch)...]
    P = typeof(_cluster_position(branch, zero(T), 0))
    candidates = Vector{_ClusterCandidate{P}}(undef, length(active) * (maxell + 1))
    k = 0
    for family in eachindex(active), ell in 0:maxell
        pos = _cluster_position(branch, active[family], ell)
        candidates[k += 1] = (family=family, ell=ell, position=pos)
    end
    return candidates
end

"""
    _cluster_kind(positions; near_radius)

Classify a group of pole positions as `:exact`, `:near`, or `:separated`.
"""
function _cluster_kind(positions; near_radius)
    length(positions) <= 1 && return :separated
    diameter = maximum(abs(positions[i] - positions[j])
                       for i in eachindex(positions) for j in eachindex(positions))
    center = sum(positions) / length(positions)
    diameter <= _cluster_exact_tol(center) && return :exact
    diameter <= near_radius && return :near
    return :separated
end

"""
    _make_active_cluster(a, b, m, n, branch, members; near_radius)

Build a structured active-cluster record from grouped candidate members.
"""
function _make_active_cluster(a::Tuple, b::Tuple, m::Integer, n::Integer,
                              branch::Symbol, members::Vector{_ClusterCandidate{P}};
                              near_radius,
                              den=_cluster_denominator_factors(a, b, m, n)) where {P}
    positions = P[member.position for member in members]
    center = sum(positions) / length(positions)
    diameter = maximum(abs(x - y) for x in positions for y in positions)
    kind = _cluster_kind(positions; near_radius)
    families = Int[member.family for member in members]
    ells = Int[member.ell for member in members]
    shifts = Int[round(Int, real(positions[k] - positions[1])) for k in eachindex(positions)]
    return _ActiveCluster(kind, branch, length(members), center, diameter, members,
                          families, ells, shifts, _cluster_has_denominator_pole(den, center))
end

# The singleton record of one candidate, built without the pairwise geometry.
function _make_singleton_cluster(branch::Symbol, member::_ClusterCandidate{P}, den) where {P}
    return _ActiveCluster(:exact, branch, 1, member.position, zero(abs(member.position)),
                          [member], Int[member.family], Int[member.ell], Int[0],
                          _cluster_has_denominator_pole(den, member.position))
end

# Greedy grouping of sorted candidates by position, each unused candidate seeding a group.
# A function barrier so the quadratic loop runs on the concretely typed candidate vector.
# Returns the member indices of each group in seed order.
function _cluster_groups(candidates::Vector{NT}, near_radius) where {NT}
    used = falses(length(candidates))
    groups = Vector{Int}[]
    @inbounds for i in eachindex(candidates)
        used[i] && continue
        seed = candidates[i]
        group = Int[i]
        used[i] = true
        tol = max(near_radius, _cluster_exact_tol(seed.position))
        for j in (i + 1):length(candidates)
            used[j] && continue
            if abs(candidates[j].position - seed.position) <= tol
                push!(group, j)
                used[j] = true
            end
        end
        push!(groups, group)
    end
    return groups
end

"""
    _detect_active_confluence_clusters(a, b, m, n, branch; maxell, near_radius, include_singletons)

Detect active-pole clusters of arbitrary multiplicity for one branch.  The
detector reports all exact or near clusters it sees. Evaluator support is a
separate policy decision.
"""
function _detect_active_confluence_clusters(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                            branch::Symbol;
                                            maxell::Integer=_CONFLUENCE_CLUSTER_MAXELL,
                                            near_radius=_CONFLUENCE_CLUSTER_NEAR_RADIUS,
                                            include_singletons::Bool=false)
    candidates = _cluster_candidates(a, b, m, n, branch; maxell)
    sort!(candidates; by=c -> (real(c.position), imag(c.position)))
    P = fieldtype(eltype(candidates), :position)
    clusters = _ActiveCluster{P,typeof(abs(zero(P)))}[]
    den = _cluster_denominator_factors(a, b, m, n)   # once, not per cluster
    for group in _cluster_groups(candidates, near_radius)
        if include_singletons && length(group) == 1
            push!(clusters, _make_singleton_cluster(branch, candidates[only(group)], den))
        else
            members = candidates[group]
            cluster = _make_active_cluster(a, b, m, n, branch, members; near_radius, den)
            cluster.kind === :separated || push!(clusters, cluster)
        end
    end
    sort!(clusters; by=c -> (real(c.center), imag(c.center), c.multiplicity))
    return clusters
end

"""
    _active_family_pair(branch, active, i, j)

Return the nearest integer shift and residual between two active families.
"""
function _active_family_pair(branch::Symbol, active, i::Integer, j::Integer)
    δ = active[j] - active[i]
    shift = round(Int, real(δ))
    residual = δ - shift
    return (i=i, j=j, shift=shift, residual=residual,
            abs_residual=abs(residual))
end

"""
    _family_cluster_center(branch, active, families, shifts)

Return one pole location where all families in an exact shifted cluster align.
For near clusters this is the corresponding reference location.
"""
function _family_cluster_center(branch::Symbol, active, families, shifts)
    if branch === :lower
        maxshift = maximum(shifts)
        ells = Tuple(maxshift - shifts[k] for k in eachindex(shifts))
        center = active[families[1]] + maxshift
        return (center=center, ells=ells)
    end
    minshift = minimum(shifts)
    ells = Tuple(shifts[k] - minshift for k in eachindex(shifts))
    center = active[families[1]] + minshift - one(active[families[1]])
    return (center=center, ells=ells)
end

"""
    _make_active_family_cluster(a, b, m, n, branch, active, families; near_radius)

Build a confluence cluster of active parameter families.  This is the dispatch
level object: one shifted active-family cluster may generate many residue
clusters in the Slater series.
"""
function _make_active_family_cluster(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                     branch::Symbol, active, families;
                                     near_radius)
    base = active[first(families)]
    shifts = Tuple(round(Int, real(active[f] - base)) for f in families)
    residuals = Tuple(active[families[k]] - base - shifts[k]
                      for k in eachindex(families))
    diameter = maximum(abs(residuals[i] - residuals[j])
                       for i in eachindex(residuals) for j in eachindex(residuals))
    exact_tol = maximum(_cluster_exact_tol(active[f]) for f in families)
    kind = diameter <= exact_tol ? :exact : diameter <= near_radius ? :near : :separated
    loc = _family_cluster_center(branch, active, families, shifts)
    members = Tuple((family=families[k], ell=loc.ells[k]) for k in eachindex(families))
    return (kind=kind, branch=branch, multiplicity=length(families),
            center=loc.center, diameter=diameter, families=Tuple(families),
            shifts=shifts, residuals=residuals, ells=loc.ells, members=members,
            denominator_pole=_cluster_has_denominator_pole(a, b, m, n, loc.center))
end

"""
    _detect_active_family_confluence_clusters(a, b, m, n, branch; near_radius)

Detect exact or near integer-separated active parameter-family clusters of
arbitrary multiplicity.  This is the classifier intended for dispatch.
"""
function _detect_active_family_confluence_clusters(a::Tuple, b::Tuple,
        m::Integer, n::Integer, branch::Symbol;
        near_radius=_CONFLUENCE_CLUSTER_NEAR_RADIUS)
    active = _cluster_active_values(a, b, m, n, branch)
    length(active) <= 1 && return Any[]
    parent = collect(eachindex(active))

    """Find the representative of a family-cluster component."""
    function findroot(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end

    """Merge two family-cluster components."""
    function unionroot(i, j)
        ri = findroot(i)
        rj = findroot(j)
        ri == rj && return nothing
        parent[rj] = ri
        return nothing
    end

    for i in eachindex(active), j in (i + 1):length(active)
        pair = _active_family_pair(branch, active, i, j)
        pair.abs_residual <= near_radius && unionroot(i, j)
    end

    groups = Dict{Int,Vector{Int}}()
    for i in eachindex(active)
        root = findroot(i)
        push!(get!(groups, root, Int[]), i)
    end

    clusters = Any[]
    for families in values(groups)
        length(families) <= 1 && continue
        sort!(families)
        cluster = _make_active_family_cluster(a, b, m, n, branch, active, families;
                                              near_radius)
        cluster.kind === :separated || push!(clusters, cluster)
    end
    sort!(clusters; by=c -> (real(c.center), imag(c.center), c.multiplicity))
    return clusters
end

"""
    _active_confluence_summary(a, b, m, n, z; branch=:auto, maxell, near_radius)

Return a geometry summary for active-pole confluence risk.  This is used by
tests and future dispatch planning. It does not perform numerical evaluation.
"""
function _active_confluence_summary(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
        branch::Symbol=:auto, maxell::Integer=_CONFLUENCE_CLUSTER_MAXELL,
        near_radius=_CONFLUENCE_CLUSTER_NEAR_RADIUS)
    p = length(a)
    q = length(b)
    p == q && return (status=:unsupported_finite_p_eq_q, branch=:unknown,
                      clusters=(), max_multiplicity=0)
    br = branch === :auto ? (p < q ? :lower : :upper) : branch
    br in (:lower, :upper) || return (status=:unsupported_branch, branch=br,
                                      clusters=(), max_multiplicity=0)
    clusters = _detect_active_family_confluence_clusters(a, b, m, n, br;
                                                         near_radius)
    isempty(clusters) && return (status=:not_near_confluent, branch=br,
                                 clusters=(), max_multiplicity=0)
    any(c -> c.denominator_pole, clusters) &&
        return (status=:unsupported_local_structure, branch=br,
                clusters=Tuple(clusters),
                max_multiplicity=maximum(c.multiplicity for c in clusters))
    length(clusters) > 1 &&
        return (status=:unsupported_multiple_clusters, branch=br,
                clusters=Tuple(clusters),
                max_multiplicity=maximum(c.multiplicity for c in clusters))
    cluster = only(clusters)
    return (status=:ok, branch=br, clusters=(cluster,),
            max_multiplicity=cluster.multiplicity)
end
