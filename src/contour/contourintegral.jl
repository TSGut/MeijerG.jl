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

By default (`contour_case = :auto`), case selection follows the branch and pole geometry:
- In the convergent regime of the branch (`|z| <= 1` lower, `|z| >= 1` upper): Case (ii)
  for lower mode or Case (iii) for upper mode, falling back to Case (i) if the loop fails
  and the vertical integral converges
- Otherwise, if a separating vertical strip exists and the vertical integral converges,
  `|arg z| < (m + n - (p + q)/2) π` (DLMF 16.17(i)): Case (i), falling back to the
  loop contour of the branch if the quadrature fails
- Otherwise: Case (ii) for lower mode or Case (iii) for upper mode

# Keyword arguments
- `rtol`: relative tolerance for quadrature and truncation refinement
- `atol`: absolute tolerance for quadrature and truncation refinement
- `contour_case`: one of `:auto`, `:i`, `:ii`, `:iii`
- `sigma`: real part of the Barnes contour. If omitted, a default value inside the admissible strip is chosen
- `initial_halfwidth`: initial truncation half-width in the `t` variable
- `max_truncation_steps`: maximum number of half-width doublings used to resolve the infinite tails.
  For the loop cases it sets the pole budget `max(8, initial_halfwidth) * 2^(max_truncation_steps - 1)`
  levels per family.

# Notes
- Real negative `z` must be passed explicitly as complex (`z + 0im`) to select the desired branch.
- Loop contours group nearby poles when the surrounding gap permits. Each circle is integrated
  with the trapezoidal rule, which converges geometrically for the periodic integrand, and pole
  levels are added until the sum is stable to a few ulps. Acceptance compares the accumulated
  circle error estimates, including the rounding scale of each circle, with
  `max(atol, rtol * max(abs(value), 1))`. These are numerical checks, not rigorous error bounds.
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

    vertical() = _contour_case_i(a_promoted, b_promoted, m, n, z_promoted, contour_sigma, contour_rtol, contour_atol, halfwidth, max_truncation_steps)

    value = if selected_case === :i && contour_case === :auto
        # The loop contour of the branch is the fallback when the vertical quadrature fails.
        try
            vertical()
        catch err
            err isa Union{DomainError,ErrorException} || rethrow()
            loop(_expansion_mode(p, q, z_promoted) === :lower ? :ii : :iii)
        end
    elseif contour_case === :auto && _has_vertical_strip(a_promoted, b_promoted, m, n) &&
           _vertical_contour_converges(a_promoted, b_promoted, m, n, z_promoted)
        # A loop preferred in the convergent regime falls back to the vertical quadrature.
        try
            loop(selected_case)
        catch err
            err isa Union{DomainError,ErrorException} || rethrow()
            vertical()
        end
    elseif selected_case === :i
        vertical()
    else
        loop(selected_case)
    end
    return _output_like(value, a, b, z)
end

"""
    _resolve_contour_case(contour_case::Symbol, a::Tuple, b::Tuple, m::Integer, n::Integer, z) -> Symbol

Resolve contour-case selection to one of `:i`, `:ii`, or `:iii`.

- `:auto`: chooses the loop case of the expansion mode (`:ii` for lower, `:iii` for
  upper) when `z` is in the convergent regime of that branch
  (`_contour_loop_in_regime`). Otherwise it chooses `:i` when a vertical strip exists
  and the vertical integral converges, and the loop case of the mode if not.
- `:i`, `:ii`, `:iii`: explicit case selection.
"""
function _resolve_contour_case(contour_case::Symbol, a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    contour_case in (:auto, :i, :ii, :iii) || throw(ArgumentError("contour_case must be one of :auto, :i, :ii, :iii"))

    if contour_case === :auto
        mode = _expansion_mode(length(a), length(b), z)
        loop_case = mode === :lower ? :ii : :iii
        _contour_loop_in_regime(mode, z) && return loop_case
        if _has_vertical_strip(a, b, m, n) && _vertical_contour_converges(a, b, m, n, z)
            return :i
        end
        return loop_case
    end

    return contour_case
end

# Convergent regime of a branch (lower small |z|, upper large |z|). There the residue sum of the
# loop contour barely cancels, and the trapezoidal loops reach full accuracy at about the cost of
# the vertical quadrature at its default tolerance.
_contour_loop_in_regime(mode::Symbol, z) = mode === :lower ? abs(z) <= 1 : abs(z) >= 1

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

"""
    _contour_pole_cluster_sum(a, b, m, n, z, family; rtol, atol, halfwidth, max_truncation_steps)

Sum the loop integrals around the active pole clusters of `family` (`:lower` or
`:upper`), adding one pole level at a time.  Circles already integrated are
reused, so each circle is evaluated once.  The sum is truncated after two
successive levels change it by at most a few ulps, counted only from
`_contour_minimum_levels` on, and is accepted when the accumulated circle error
estimates, the rounding scale of each circle included, are within
`max(atol, rtol * max(abs(value), 1))`.  The pole budget is that of the former
doubling schedule, `max(8, halfwidth) * 2^(max_truncation_steps - 1)` levels.
"""
function _contour_pole_cluster_sum(a::Tuple, b::Tuple, m::Integer, n::Integer, z, family::Symbol;
                                   rtol, atol, halfwidth, max_truncation_steps::Integer)
    logz = log(z)
    total = zero(complex(_contour_accumulator_seed(a, b, z)))
    R = typeof(real(total))
    max_level = max(8, ceil(Int, Float64(halfwidth))) << (max_truncation_steps - 1)
    first_level = min(_contour_minimum_levels(a, b, m, n, family), max_level)
    C = typeof(first(family === :lower ? b : a) + 0)
    circles = Dict{Tuple{C, R}, Tuple{typeof(total), R}}()
    previous = total
    settled = 0
    for maxell in first_level:max_level
        total = zero(total)
        error_sum = zero(R)
        for cluster in _contour_active_pole_clusters(a, b, m, n, family, maxell)
            radius = convert(R, _contour_cluster_radius(cluster.center, cluster.members,
                                                        a, b, m, n, maxell))
            value, error = get!(circles, (cluster.center, radius)) do
                _contour_circle_trapezoid(cluster.center, radius, a, b, m, n, logz)
            end
            total += value
            error_sum += error
        end
        if maxell > first_level
            step = abs(total - previous)
            settled = step <= max(atol, 4 * eps(R) * abs(total)) ? settled + 1 : 0
            if settled >= 2
                value = family === :lower ? -total : total
                target = max(atol, rtol * max(abs(value), one(R)))
                isfinite(error_sum) && error_sum + step <= target && return value
                throw(ErrorException("Loop contour quadrature or rounding check failed."))
            end
        end
        previous = total
    end

    throw(ErrorException("Loop contour truncation or quadrature checks did not converge within the configured refinement budget."))
end

# Residues of an active family can be suppressed, exactly or nearly, on its leading levels by the
# zeros of a passive reciprocal gamma factor on the same side: 1/Γ(1 - b_k + s) for the lower
# poles b_j + ℓ up to ℓ = Re(b_k - b_j) - 1, and 1/Γ(a_k - s) for the upper poles a_j - 1 - ℓ up
# to ℓ = Re(a_j - a_k) - 1.  The truncation test starts after these levels, so that small leading
# increments cannot end the sum.
function _contour_minimum_levels(a::Tuple, b::Tuple, m::Integer, n::Integer, family::Symbol)
    first_level = 2
    if family === :lower
        for j in 1:m, k in (m + 1):length(b)
            first_level = max(first_level, floor(Int, real(b[k] - b[j])) + 1)
        end
    else
        for j in 1:n, k in (n + 1):length(a)
            first_level = max(first_level, floor(Int, real(a[j] - a[k])) + 1)
        end
    end
    return first_level
end

"""
    _contour_circle_trapezoid(center, radius, a, b, m, n, logz)

Return `(value, error)` for `(1/2πi)∮ K(s) ds` over `|s - center| = radius` by the
trapezoidal rule.  The integrand is periodic and analytic in an annulus around the
circle, so the rule converges geometrically, with ratio
`max(reach / radius, radius / nearest)` for the cluster bounds of
`_contour_cluster_radius`, at most 1/2.  The number of points is doubled, reusing
the previous points, until two estimates agree to a few ulps of the largest
sample. `error` is that difference plus one ulp of the largest sample.
"""
function _contour_circle_trapezoid(center, radius, a::Tuple, b::Tuple, m::Integer, n::Integer, logz)
    R = typeof(radius)
    sample(k, points) = (phase = cispi(convert(R, 2k) / points);
        _contour_kernel(center + radius * phase, a, b, m, n, logz) * (radius * phase))
    points = 8
    max_points = 16 * nextpow(2, precision(R))
    total = sample(0, points)
    scale = abs(total)
    for k in 1:(points - 1)
        v = sample(k, points)
        total += v
        scale = max(scale, abs(v))
    end
    estimate = total / points
    while points < max_points
        for k in 1:2:(2points - 1)
            v = sample(k, 2points)
            total += v
            scale = max(scale, abs(v))
        end
        points *= 2
        refined = total / points
        difference = abs(refined - estimate)
        isfinite(difference) || break
        difference <= 16 * eps(R) * scale && return refined, difference + eps(R) * scale
        estimate = refined
    end
    return estimate, R(Inf)
end

function _contour_accumulator_seed(a::Tuple, b::Tuple, z)
    seed = z
    isempty(a) || (seed += first(a))
    isempty(b) || (seed += first(b))
    return seed
end

# Pole-cluster bookkeeping uses concretely-typed centers and member tuples.
# Keep centers and member lists concrete as the truncation grows.
const _ContourMember = Tuple{Symbol, Int, Int}

function _contour_active_pole_clusters(a::Tuple, b::Tuple, m::Integer, n::Integer, family::Symbol, maxell::Integer)
    parameters, count, direction = family === :lower ? (b, m, 1) : (a, n, -1)
    C = typeof(first(parameters) + 0)
    centers = C[]
    members = _ContourMember[]
    for j in 1:count, ell in 0:maxell
        base = family === :lower ? parameters[j] : parameters[j] - 1
        push!(centers, base + direction * ell)
        push!(members, (family, j, ell))
    end
    required_count = length(centers)
    included = Set(members)
    # Do not cut a multiple/nearby pole group at the artificial tail boundary.
    # Add the local halo, but do not start new circles from halo-only poles.
    for i in 1:required_count, j in 1:count
        base = family === :lower ? parameters[j] : parameters[j] - 1
        ell = max(0, round(Int, direction * real(centers[i] - base)))
        member = (family, j, ell)
        pole = base + direction * ell
        if !(member in included) && abs(pole - centers[i]) <= 0.2
            push!(centers, pole)
            push!(members, member)
            push!(included, member)
        end
    end
    return _contour_cluster_candidates(centers, members; required_count=required_count,
        admissible=(center, group) -> begin
            reach, nearest = _contour_cluster_bounds(center, group, a, b, m, n)
            reach < nearest / 4
        end)
end

function _contour_cluster_candidates(centers::Vector{C}, members::Vector{_ContourMember};
                                     admissible=(center, group) -> true,
                                     required_count=length(centers)) where {C}
    N = length(centers)
    used = falses(N)
    clusters = Vector{NamedTuple{(:center, :members), Tuple{C, Vector{_ContourMember}}}}()
    # The pole lattices have unit spacing. Keep groups local, but join a tight
    # group well before floating-point coincidence if its surrounding gap allows it.
    @inbounds for i in 1:required_count
        used[i] && continue
        center = centers[i]
        order = sortperm(abs.(centers .- center))
        candidates = Int[]
        for j in order
            (used[j] || abs(centers[j] - center) > 0.2) && break
            push!(candidates, j)
        end
        accepted = false
        for count in length(candidates):-1:1
            selected = candidates[1:count]
            reach = maximum(abs(centers[j] - center) for j in selected)
            nearest = count < N ? abs(centers[order[count + 1]] - center) : Inf
            reach < nearest / 4 || continue
            group = members[selected]
            admissible(center, group) || continue
            used[selected] .= true
            push!(clusters, (center=center, members=group))
            accepted = true
            break
        end
        accepted || throw(ErrorException("Cannot isolate a loop pole cluster with adequate clearance."))
    end
    return clusters
end

# Distance to every assigned pole and the nearest unassigned pole, including
# poles beyond the finite truncation. Each family is an infinite unit lattice;
# search outward from its nearest nonnegative index until a nonmember is found.
function _contour_cluster_bounds(center, members::Vector{_ContourMember}, a::Tuple, b::Tuple, m::Integer, n::Integer)
    member_set = Set(members)
    reach = zero(real(center))
    for (family, j, ell) in members
        pole = family === :lower ? b[j] + ell : a[j] - 1 - ell
        reach = max(reach, abs(pole - center))
    end
    nearest = oftype(reach, Inf)
    for (family, parameters, count, direction) in ((:lower, b, m, 1), (:upper, a, n, -1))
        for j in 1:count
            base = family === :lower ? parameters[j] : parameters[j] - 1
            k = max(0, round(Int, direction * real(center - base)))
            for step in (-1, 1)
                ell = k
                while ell >= 0 && (family, j, ell) in member_set
                    ell += step
                end
                ell < 0 && continue
                nearest = min(nearest, abs(base + direction * ell - center))
            end
        end
    end
    return reach, nearest
end

function _contour_cluster_radius(center, members::Vector{_ContourMember}, a::Tuple, b::Tuple, m::Integer, n::Integer, maxell::Integer)
    reach, nearest = _contour_cluster_bounds(center, members, a, b, m, n)
    reach < nearest / 4 || throw(ErrorException("Loop pole cluster has insufficient clearance."))
    radius = max(2 * reach, min(oftype(reach, 0.2), nearest * oftype(reach, 0.35)))
    reach < radius < nearest || throw(ErrorException("Invalid loop pole cluster radius."))
    return radius
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
