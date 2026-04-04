using QuadGK: quadgk

@doc raw"""
    meijerg_contour(a, b, m, n, z; kwargs...)
    meijerg_contour(a_left, a_right, b_left, b_right, z; kwargs...)

Evaluate the Meijer G-function by direct numerical integration of the Mellin-Barnes contour.

**Warning**: This is a low-level API intended for power users and validation workflows. For normal use,
call `meijerg()` instead, which uses adaptive perturbation for confluent cases and deploys contour 
integration only as an emergency fallback.

This method numerically integrates the defining contour integral using explicit DLMF
contour families (DLMF 16.17):
- Case (i): vertical Barnes contour `s = σ + i t` (robust when strip exists)
- Case (ii): loop contour around lower poles of `Γ(b_j - s)` (DLMF Case ii)
- Case (iii): loop contour around upper poles of `Γ(1 - a_j + s)` (DLMF Case iii)

By default (`contour_case = :auto`), case selection follows pole geometry:
- If a separating vertical strip exists: Case (i)
- Otherwise: Case (ii) for lower mode or Case (iii) for upper mode

# Keyword arguments
- `rtol`: relative tolerance for quadrature and truncation refinement
- `atol`: absolute tolerance for quadrature and truncation refinement
- `contour_case`: one of `:auto`, `:i`, `:ii`, `:iii`
- `sigma`: real part of the Barnes contour; if omitted, a default value inside the admissible strip is chosen
- `initial_halfwidth`: initial truncation half-width in the `t` variable
- `max_truncation_steps`: maximum number of half-width doublings used to resolve the infinite tails
- `allow_regularization_fallback`: (Optional, default: false) if `true`, falls back to perturbative regularization if contour integration fails

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
                         max_truncation_steps::Integer=8,
                         allow_regularization_fallback::Bool=false)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg_contour is implemented for nonzero z only"))
    z isa Real && z < zero(z) && throw(DomainError(z, "Contour evaluation on the negative real axis requires an explicit complex argument; pass z + 0im."))
    max_truncation_steps >= 2 || throw(ArgumentError("max_truncation_steps must be at least 2"))

    _validate_pairing(a, b, m, n)
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    real_type = typeof(real(z_promoted))
    contour_sigma = sigma === nothing ? _default_contour_sigma(a_promoted, b_promoted, m, n, real_type) : convert(real_type, sigma)

    contour_rtol = _contour_tolerance(real_type, rtol, sqrt(eps(real_type)))
    contour_atol = _contour_tolerance(real_type, atol, zero(real_type))
    halfwidth = initial_halfwidth === nothing ? _default_contour_halfwidth(real_type, z_promoted) : convert(real_type, initial_halfwidth)

    selected_case = _resolve_contour_case(contour_case, a_promoted, b_promoted, m, n, z_promoted)

    try
        if selected_case === :i
            return _contour_case_i(a_promoted, b_promoted, m, n, z_promoted, contour_sigma, contour_rtol, contour_atol, halfwidth, max_truncation_steps)
        elseif selected_case === :ii
            return _contour_loop_lower(a_promoted, b_promoted, m, n, z_promoted; rtol=contour_rtol, atol=contour_atol, halfwidth=halfwidth)
        else
            return _contour_loop_upper(a_promoted, b_promoted, m, n, z_promoted; rtol=contour_rtol, atol=contour_atol, halfwidth=halfwidth)
        end
    catch
        allow_regularization_fallback || rethrow()
        return _meijerg_perturb(a_promoted, b_promoted, m, n, z_promoted)
    end
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
        if _has_vertical_strip(a, b, m, n)
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
"""
function _contour_kernel(s, a::Tuple, b::Tuple, m::Integer, n::Integer, logz)
    kernel = _gamma_prod_ratio(
        Iterators.flatten(((b[j] - s for j in 1:m),
                           (one(s) - a[j] + s for j in 1:n))),
        Iterators.flatten(((one(s) - b[j] + s for j in m+1:length(b)),
                           (a[j] - s for j in n+1:length(a)))),
        typeof(s),
    )
    return kernel * exp(s * logz)
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
function _contour_loop_lower(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol, atol, halfwidth)
    m == 0 && return zero(z)

    real_type = typeof(real(z))
    logz = log(z)

    pole_reals = [real(b[j]) for j in 1:m]
    pole_max = maximum(pole_reals)
    pole_min = minimum(pole_reals)
    pole_center = (pole_min + pole_max) / 2
    pole_radius = (pole_max - pole_min) / 2 + convert(real_type, 0.6)

    σ_right = pole_center + pole_radius + one(real_type)
    σ_left = pole_center - pole_radius - one(real_type)
    L = max(convert(real_type, halfwidth), convert(real_type, 8), convert(real_type, 4) + log1p(convert(real_type, abs(z))))

    integrand_right = t -> begin
        s = Complex{real_type}(σ_right, t)
        _contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi))
    end

    integrand_left = t -> begin
        s = Complex{real_type}(σ_left, t)
        _contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi))
    end

    integrand_bottom = u -> begin
        angle = -convert(real_type, pi) * u
        s = Complex{real_type}(pole_center + pole_radius * cos(angle), -L + pole_radius * sin(angle))
        ds_du = pole_radius * (-convert(real_type, pi)) * Complex{real_type}(-sin(angle), cos(angle))
        _contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi)) * ds_du
    end

    integrand_top = u -> begin
        angle = convert(real_type, pi) * u
        s = Complex{real_type}(pole_center + pole_radius * cos(angle), L + pole_radius * sin(angle))
        ds_du = pole_radius * convert(real_type, pi) * Complex{real_type}(-sin(angle), cos(angle))
        -_contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi)) * ds_du
    end

    val_right, _ = quadgk(integrand_right, -L, L; rtol=rtol, atol=atol, maxevals=10000)
    val_left, _ = quadgk(integrand_left, -L, L; rtol=rtol, atol=atol, maxevals=10000)
    val_bottom, _ = quadgk(integrand_bottom, zero(real_type), one(real_type); rtol=rtol, atol=atol, maxevals=10000)
    val_top, _ = quadgk(integrand_top, zero(real_type), one(real_type); rtol=rtol, atol=atol, maxevals=10000)

    return val_right + val_left + val_bottom + val_top
end

"""
    _contour_loop_upper(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol, atol, halfwidth)

Evaluate DLMF Case (iii): loop contour around upper poles of `Γ(1 - a_j + s)`.
"""
function _contour_loop_upper(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol, atol, halfwidth)
    n == 0 && return zero(z)

    real_type = typeof(real(z))
    logz = log(z)

    pole_reals = [real(one(real_type) - a[j]) for j in 1:n]
    pole_max = maximum(pole_reals)
    pole_min = minimum(pole_reals)
    pole_center = (pole_min + pole_max) / 2
    pole_radius = (pole_max - pole_min) / 2 + convert(real_type, 0.6)

    σ_left = pole_center - pole_radius - one(real_type)
    σ_right = pole_center + pole_radius + one(real_type)
    L = max(convert(real_type, halfwidth), convert(real_type, 8), convert(real_type, 4) + log1p(convert(real_type, abs(z))))

    integrand_left = t -> begin
        s = Complex{real_type}(σ_left, t)
        _contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi))
    end

    integrand_right = t -> begin
        s = Complex{real_type}(σ_right, t)
        _contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi))
    end

    integrand_bottom = u -> begin
        angle = convert(real_type, pi) * u
        s = Complex{real_type}(pole_center + pole_radius * cos(angle), -L + pole_radius * sin(angle))
        ds_du = pole_radius * convert(real_type, pi) * Complex{real_type}(-sin(angle), cos(angle))
        _contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi)) * ds_du
    end

    integrand_top = u -> begin
        angle = -convert(real_type, pi) * u
        s = Complex{real_type}(pole_center + pole_radius * cos(angle), L + pole_radius * sin(angle))
        ds_du = pole_radius * (-convert(real_type, pi)) * Complex{real_type}(-sin(angle), cos(angle))
        -_contour_kernel(s, a, b, m, n, logz) / (2 * convert(real_type, pi)) * ds_du
    end

    val_left, _ = quadgk(integrand_left, -L, L; rtol=rtol, atol=atol, maxevals=10000)
    val_right, _ = quadgk(integrand_right, -L, L; rtol=rtol, atol=atol, maxevals=10000)
    val_bottom, _ = quadgk(integrand_bottom, zero(real_type), one(real_type); rtol=rtol, atol=atol, maxevals=10000)
    val_top, _ = quadgk(integrand_top, zero(real_type), one(real_type); rtol=rtol, atol=atol, maxevals=10000)

    return val_left + val_right + val_bottom + val_top
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
