using QuadGK: quadgk

@doc raw"""
    meijerg_contour(a, b, m, n, z; kwargs...)
    meijerg_contour(a_left, a_right, b_left, b_right, z; kwargs...)

Evaluate the Meijer G-function by direct numerical integration of the Mellin-Barnes contour.

**Warning**: This is a low-level API intended for power users and validation workflows. For normal use,
call `meijerg()` instead, which uses reductions, guarded log-ratio jets, and Slater before contour
integration.

This method numerically integrates the defining contour integral using explicit DLMF
contour families (DLMF 16.17):
- Case (i): vertical Barnes contour `s = σ + i t` (robust when strip exists)
- Case (ii): loop contour around lower poles of `Γ(b_j - s)` (DLMF Case ii)
- Case (iii): loop contour around upper poles of `Γ(1 - a_j + s)` (DLMF Case iii)

By default (`contour_case = :auto`), case selection follows pole geometry:
- If a separating vertical strip exists and the vertical integral converges,
  `|arg z| < (m + n - (p + q)/2) π` (DLMF 16.17(i)): Case (i), falling back to the
  loop contour of the branch if the quadrature fails
- Otherwise: Case (ii) for lower mode or Case (iii) for upper mode

# Keyword arguments
- `rtol`: relative tolerance for quadrature and truncation refinement
- `atol`: absolute tolerance for quadrature and truncation refinement
- `contour_case`: one of `:auto`, `:i`, `:ii`, `:iii`
- `sigma`: real part of the Barnes contour. If omitted, a default value inside the admissible strip is chosen
- `initial_halfwidth`: initial truncation half-width in the `t` variable
- `max_truncation_steps`: maximum number of half-width doublings used to resolve the infinite tails

# Notes
- Real negative `z` must be passed explicitly as complex (`z + 0im`) to select the desired branch.
"""
function meijerg_contour(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z; kwargs...)
    meijerg_contour(Tuple(a), Tuple(b), m, n, z; kwargs...)
end

function meijerg_contour(a_left::ParameterInput, a_right::ParameterInput,
                         b_left::ParameterInput, b_right::ParameterInput, z; kwargs...)
    a_left_tuple = Tuple(a_left)
    a_right_tuple = Tuple(a_right)
    b_left_tuple = Tuple(b_left)
    b_right_tuple = Tuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg_contour(a, b, length(b_left_tuple), length(a_left_tuple), z; kwargs...)
end

function meijerg_contour(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
                         rtol=nothing,
                         atol=nothing,
                         contour_case::Symbol=:auto,
                         sigma=nothing,
                         initial_halfwidth=nothing,
                         max_truncation_steps::Integer=8)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg_contour is implemented for nonzero z only"))
    z isa Real && z < zero(z) && throw(DomainError(z, "Contour evaluation on the negative real axis requires an explicit complex argument; pass z + 0im."))
    max_truncation_steps >= 2 || throw(ArgumentError("max_truncation_steps must be at least 2"))

    _validate_pairing(a, b, m, n)
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, float(z))
    real_type = typeof(real(z_promoted))
    contour_sigma = sigma === nothing ? _default_contour_sigma(a_promoted, b_promoted, m, n, real_type) : convert(real_type, sigma)

    contour_rtol = _contour_tolerance(real_type, rtol, sqrt(eps(real_type)))
    contour_atol = _contour_tolerance(real_type, atol, zero(real_type))
    halfwidth = initial_halfwidth === nothing ? _default_contour_halfwidth(real_type, z_promoted) : convert(real_type, initial_halfwidth)

    selected_case = _resolve_contour_case(contour_case, a_promoted, b_promoted, m, n, z_promoted)
    loop(case) = case === :ii ?
        _contour_loop_lower(a_promoted, b_promoted, m, n, z_promoted; rtol=contour_rtol, atol=contour_atol, halfwidth=halfwidth, max_truncation_steps=max_truncation_steps) :
        _contour_loop_upper(a_promoted, b_promoted, m, n, z_promoted; rtol=contour_rtol, atol=contour_atol, halfwidth=halfwidth, max_truncation_steps=max_truncation_steps)

    value = if selected_case === :i && contour_case === :auto
        # The loop contour of the branch is the fallback when the vertical quadrature fails.
        try
            _contour_case_i(a_promoted, b_promoted, m, n, z_promoted, contour_sigma, contour_rtol, contour_atol, halfwidth, max_truncation_steps)
        catch err
            err isa Union{DomainError,ErrorException} || rethrow()
            loop(_expansion_mode(p, q, z_promoted) === :lower ? :ii : :iii)
        end
    elseif selected_case === :i
        _contour_case_i(a_promoted, b_promoted, m, n, z_promoted, contour_sigma, contour_rtol, contour_atol, halfwidth, max_truncation_steps)
    else
        loop(selected_case)
    end
    return _output_like(value, a, b, z)
end

"""
    _resolve_contour_case(contour_case::Symbol, a::Tuple, b::Tuple, m::Integer, n::Integer, z) -> Symbol

Resolve contour-case selection to one of `:i`, `:ii`, or `:iii`.

- `:auto`: chooses `:i` when a vertical strip exists, otherwise chooses loop case
  by expansion mode (`:ii` for lower, `:iii` for upper).
- `:i`, `:ii`, `:iii`: explicit case selection.
"""
function _resolve_contour_case(contour_case::Symbol, a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    contour_case in (:auto, :i, :ii, :iii) || throw(ArgumentError("contour_case must be one of :auto, :i, :ii, :iii"))

    if contour_case === :auto
        if _has_vertical_strip(a, b, m, n) && _vertical_contour_converges(a, b, m, n, z)
            return :i
        end
        mode = _expansion_mode(length(a), length(b), z)
        return mode === :lower ? :ii : :iii
    end

    return contour_case
end

"""
    _contour_case_i(a::Tuple, b::Tuple, m::Integer, n::Integer, z,
                    sigma, rtol, atol, halfwidth, max_truncation_steps::Integer)

Evaluate DLMF Case (i): vertical Barnes contour integration.
"""
function _contour_case_i(a::Tuple, b::Tuple, m::Integer, n::Integer, z,
                         sigma, rtol, atol, halfwidth, max_truncation_steps::Integer)
    _validate_vertical_contour(sigma, a, b, m, n)
    logz = log(z)
    integrand = t -> _contour_line_integrand(t, sigma, a, b, m, n, logz)
    return _quadgk_symmetric_truncation(integrand, halfwidth; rtol=rtol, atol=atol, max_truncation_steps=max_truncation_steps)
end

"""
    _contour_line_integrand(t, sigma, a::Tuple, b::Tuple, m::Integer, n::Integer, logz)

Return the Mellin-Barnes integrand evaluated on the vertical contour
`s = sigma + i*t`, including the `1/(2π)` normalization factor.
"""
function _contour_line_integrand(t, sigma, a::Tuple, b::Tuple, m::Integer, n::Integer, logz)
    real_type = typeof(sigma)
    s = Complex{real_type}(sigma, t)
    return _contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi))
end

"""
    _contour_kernel(s, a::Tuple, b::Tuple, m::Integer, n::Integer, logz)

Compute the core Mellin-Barnes kernel at complex point `s`:
the gamma-product ratio times `exp(s*logz)`.

The ratio is accumulated with explicit indexed loops in the same evaluation
order, sign and pole handling as the equivalent `_gamma_prod_ratio` call over
flattened generators. The kernel runs at every quadrature node, so the
generator machinery dominated its cost.
"""
function _contour_kernel(s, a::Tuple, b::Tuple, m::Integer, n::Integer, logz)
    T = typeof(s)
    log_ratio = zero(T)
    sign = one(T)
    num_inf = false
    denom_inf = false
    @inbounds for j in 1:m
        lv, sv, bad = _loggamma_factor_eval(b[j] - s, T)
        log_ratio += lv; sign *= sv; num_inf |= bad
    end
    @inbounds for j in 1:n
        lv, sv, bad = _loggamma_factor_eval(one(s) - a[j] + s, T)
        log_ratio += lv; sign *= sv; num_inf |= bad
    end
    @inbounds for j in (m + 1):length(b)
        lv, sv, bad = _loggamma_factor_eval(one(s) - b[j] + s, T)
        log_ratio -= lv; sign *= sv; denom_inf |= bad
    end
    @inbounds for j in (n + 1):length(a)
        lv, sv, bad = _loggamma_factor_eval(a[j] - s, T)
        log_ratio -= lv; sign *= sv; denom_inf |= bad
    end
    num_inf && denom_inf && _gamma_prod_pole_error()
    return sign * exp(log_ratio) * exp(s * logz)
end

"""
    _quadgk_symmetric_truncation(f, halfwidth; rtol, atol, max_truncation_steps::Integer)

Integrate `f` over `(-∞, ∞)` by symmetric finite-interval truncation.

Starts from `[-halfwidth, halfwidth]`, doubles the half-width each iteration,
and terminates when successive integral estimates satisfy the requested
`rtol`/`atol` criterion.
"""
function _quadgk_symmetric_truncation(f, halfwidth; rtol, atol, max_truncation_steps::Integer)
    previous_value = nothing
    current_halfwidth = halfwidth

    for _ in 1:max_truncation_steps
        value, _ = quadgk(f, -current_halfwidth, current_halfwidth; rtol=rtol / 4, atol=atol / 4, maxevals=10_000)
        if previous_value !== nothing
            scale = max(abs(value), abs(previous_value), one(abs(value)))
            if abs(value - previous_value) <= max(atol, rtol * scale)
                return value
            end
        end
        previous_value = value
        current_halfwidth *= 2
    end

    throw(ErrorException("Contour integral truncation did not converge within the configured refinement budget."))
end

"""
    _contour_loop_lower(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol, atol, halfwidth)

Evaluate DLMF Case (ii): loop contour around lower poles of `Γ(b_j - s)`.
"""
function _contour_loop_lower(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol, atol, halfwidth, max_truncation_steps::Integer)
    m == 0 && return zero(z)
    return _contour_pole_cluster_sum(a, b, m, n, z, :lower;
                                     rtol=rtol, atol=atol, halfwidth=halfwidth,
                                     max_truncation_steps=max_truncation_steps)
end


"""
    _contour_loop_upper(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol, atol, halfwidth)

Evaluate DLMF Case (iii): loop contour around upper poles of `Γ(1 - a_j + s)`.
"""
function _contour_loop_upper(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol, atol, halfwidth, max_truncation_steps::Integer)
    n == 0 && return zero(z)
    return _contour_pole_cluster_sum(a, b, m, n, z, :upper;
                                     rtol=rtol, atol=atol, halfwidth=halfwidth,
                                     max_truncation_steps=max_truncation_steps)
end

function _contour_pole_cluster_sum(a::Tuple, b::Tuple, m::Integer, n::Integer, z, family::Symbol;
                                   rtol, atol, halfwidth, max_truncation_steps::Integer)
    initial_terms = max(8, ceil(Int, Float64(halfwidth)))
    previous_value = nothing
    terms = initial_terms

    for _ in 1:max_truncation_steps
        value = _contour_pole_cluster_sum_fixed(a, b, m, n, z, family, terms; rtol=rtol, atol=atol)
        if previous_value !== nothing
            scale = max(abs(value), abs(previous_value), one(abs(value)))
            if abs(value - previous_value) <= max(atol, rtol * scale)
                return value
            end
        end
        previous_value = value
        terms *= 2
    end

    throw(ErrorException("Loop contour pole-cluster truncation did not converge within the configured refinement budget."))
end

function _contour_pole_cluster_sum_fixed(a::Tuple, b::Tuple, m::Integer, n::Integer, z, family::Symbol, maxell::Integer; rtol, atol)
    logz = log(z)
    clusters = _contour_active_pole_clusters(a, b, m, n, family, maxell)
    total = zero(_contour_accumulator_seed(a, b, z))
    for cluster in clusters
        radius = _contour_cluster_radius(cluster.center, cluster.members, a, b, m, n, maxell)
        total += _contour_small_circle_integral(cluster.center, radius, a, b, m, n, logz; rtol=rtol, atol=atol)
    end
    return family === :lower ? -total : total
end

function _contour_accumulator_seed(a::Tuple, b::Tuple, z)
    seed = z
    isempty(a) || (seed += first(a))
    isempty(b) || (seed += first(b))
    return seed
end

# Pole-cluster bookkeeping uses concretely-typed centers and member tuples.
# The greedy clustering pass is quadratic in the candidate count, which grows
# with truncation doubling, so `Any[]` boxing dominated its cost.
const _ContourMember = Tuple{Symbol, Int, Int}

function _contour_active_pole_clusters(a::Tuple, b::Tuple, m::Integer, n::Integer, family::Symbol, maxell::Integer)
    if family === :lower
        count = m * (maxell + 1)
        C = typeof(@inbounds b[1] + 0)
        centers = Vector{C}(undef, count)
        members = Vector{_ContourMember}(undef, count)
        idx = 0
        @inbounds for j in 1:m, ell in 0:maxell
            idx += 1
            centers[idx] = b[j] + ell
            members[idx] = (:lower, j, ell)
        end
        return _contour_cluster_candidates(centers, members)
    else
        count = n * (maxell + 1)
        C = typeof(@inbounds a[1] + 0)
        centers = Vector{C}(undef, count)
        members = Vector{_ContourMember}(undef, count)
        idx = 0
        @inbounds for j in 1:n, ell in 0:maxell
            idx += 1
            centers[idx] = a[j] - 1 - ell
            members[idx] = (:upper, j, ell)
        end
        return _contour_cluster_candidates(centers, members)
    end
end

function _contour_cluster_candidates(centers::Vector{C}, members::Vector{_ContourMember}) where {C}
    N = length(centers)
    used = falses(N)
    clusters = Vector{NamedTuple{(:center, :members), Tuple{C, Vector{_ContourMember}}}}()
    @inbounds for i in 1:N
        used[i] && continue
        center = centers[i]
        cluster_members = _ContourMember[members[i]]
        used[i] = true
        for j in (i + 1):N
            used[j] && continue
            if abs(centers[j] - center) <= _contour_cluster_tol(center, centers[j])
                push!(cluster_members, members[j])
                used[j] = true
            end
        end
        push!(clusters, (center=center, members=cluster_members))
    end
    return clusters
end

_contour_cluster_tol(x, y) = sqrt(eps(typeof(real(x + y)))) * max(abs(x), abs(y), one(abs(x)))

function _contour_cluster_radius(center, members::Vector{_ContourMember}, a::Tuple, b::Tuple, m::Integer, n::Integer, maxell::Integer)
    member_set = Set{_ContourMember}(members)
    nearest = Inf
    @inbounds for j in 1:m, ell in 0:(maxell + 2)
        (:lower, j, ell) in member_set && continue
        nearest = min(nearest, Float64(abs((b[j] + ell) - center)))
    end
    @inbounds for j in 1:n, ell in 0:(maxell + 2)
        (:upper, j, ell) in member_set && continue
        nearest = min(nearest, Float64(abs((a[j] - 1 - ell) - center)))
    end
    real_type = typeof(real(center))
    coincident_tol = Float64(_contour_cluster_tol(center, center))
    nearest_for_radius = nearest <= coincident_tol ? Inf : nearest
    geometric_radius = isfinite(nearest_for_radius) ? min(0.2, 0.35 * nearest_for_radius) : 0.2
    floor_radius = sqrt(eps(real_type))
    radius = max(geometric_radius, Float64(floor_radius))
    isfinite(nearest_for_radius) && (radius = min(radius, 0.35 * nearest_for_radius))
    return convert(real_type, radius)
end

function _contour_small_circle_integral(center, radius, a::Tuple, b::Tuple, m::Integer, n::Integer, logz; rtol, atol)
    real_type = typeof(real(center))
    r = convert(real_type, radius)
    # Near-but-not-yet-clustered poles force tiny residue circles.  The kernel
    # then scales like 1/r, so the default contour tolerance can leave a visible
    # pre-cluster error spike even though ordinary loop circles are well
    # resolved.  Tighten only these small-circle quadratures.
    small_radius = eps(real_type)^(one(real_type) / 4)
    tight_rtol = eps(real_type)^(3 * one(real_type) / 4)
    local_rtol = r < small_radius ? min(rtol, tight_rtol) : rtol
    integrand = θ -> begin
        phase = cis(θ)
        s = center + r * phase
        ds_dθ = im * r * phase
        _contour_kernel(s, a, b, m, n, logz) * ds_dθ / (2 * convert(real_type, pi) * im)
    end
    value, _ = quadgk(integrand, zero(real_type), 2 * convert(real_type, pi); rtol=local_rtol, atol=atol, maxevals=10000)
    return value
end

"""
    _default_contour_sigma(a::Tuple, b::Tuple, m::Integer, n::Integer, ::Type{R}) where {R}

Choose a default real part `sigma` for the vertical contour.

If both active pole bounds are known, returns their midpoint. Otherwise shifts
one unit away from the known side, or `0` when no bounds are active.
"""
function _default_contour_sigma(a::Tuple, b::Tuple, m::Integer, n::Integer, ::Type{R}) where {R}
    left_bound = _left_pole_bound(a, n)
    right_bound = _right_pole_bound(b, m)

    if left_bound !== nothing && right_bound !== nothing
        if left_bound < right_bound
            return convert(R, (left_bound + right_bound) / 2)
        end
        return convert(R, (left_bound + right_bound) / 2)
    elseif left_bound !== nothing
        return convert(R, left_bound + 1)
    elseif right_bound !== nothing
        return convert(R, right_bound - 1)
    else
        return zero(R)
    end
end

"""
    _left_pole_bound(a::Tuple, n::Integer)

Return the rightmost real-location bound of active upper poles,
`max(real(a[j]) - 1, j=1:n)`, or `nothing` when `n == 0`.
"""
_left_pole_bound(a::Tuple, n::Integer) = n == 0 ? nothing : maximum(real(a[j]) - 1 for j in 1:n)

"""
    _right_pole_bound(b::Tuple, m::Integer)

Return the leftmost real-location bound of active lower poles,
`min(real(b[j]), j=1:m)`, or `nothing` when `m == 0`.
"""
_right_pole_bound(b::Tuple, m::Integer) = m == 0 ? nothing : minimum(real(b[j]) for j in 1:m)

"""
    _has_vertical_strip(a::Tuple, b::Tuple, m::Integer, n::Integer)

Return `true` when a strict separating vertical strip exists between
active upper and lower pole families, or when one family is inactive.
"""
_has_vertical_strip(a::Tuple, b::Tuple, m::Integer, n::Integer) = begin
    left_bound = _left_pole_bound(a, n)
    right_bound = _right_pole_bound(b, m)
    left_bound === nothing || right_bound === nothing || left_bound < right_bound
end

"""
    _vertical_contour_converges(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Return `true` when the Mellin-Barnes integral along a vertical line converges,
`|arg z| < (m + n - (p + q)/2) π` (DLMF 16.17(i)). Otherwise the integrand grows
along the line and only the loop contours apply.
"""
function _vertical_contour_converges(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    δ = (2 * (m + n) - (length(a) + length(b))) / 2
    return δ > 0 && abs(angle(z)) < δ * pi
end

"""
    _validate_vertical_contour(sigma, a::Tuple, b::Tuple, m::Integer, n::Integer)

Validate that the proposed contour line `Re(s)=sigma` lies strictly between
active upper and lower pole bounds. Throws `DomainError` when invalid.
"""
function _validate_vertical_contour(sigma, a::Tuple, b::Tuple, m::Integer, n::Integer)
    left_bound = _left_pole_bound(a, n)
    right_bound = _right_pole_bound(b, m)
    left_bound !== nothing && sigma <= left_bound && throw(DomainError(sigma, "Contour sigma must lie strictly to the right of the active upper pole family."))
    right_bound !== nothing && sigma >= right_bound && throw(DomainError(sigma, "Contour sigma must lie strictly to the left of the active lower pole family."))
    return nothing
end

"""
    _default_contour_halfwidth(::Type{R}, z) where {R}

Return the initial symmetric truncation half-width for contour quadrature,
scaled mildly with `abs(z)`.
"""
function _default_contour_halfwidth(::Type{R}, z) where {R}
    magnitude = max(one(R), convert(R, abs(z)))
    return max(convert(R, 8), convert(R, 4) + log1p(magnitude))
end

"""
    _contour_tolerance(::Type{R}, value, default) where {R}

Convert an optional tolerance keyword to type `R`, using `default` when
the caller provides `nothing`.
"""
function _contour_tolerance(::Type{R}, value, default) where {R}
    value === nothing && return convert(R, default)
    return convert(R, value)
end
