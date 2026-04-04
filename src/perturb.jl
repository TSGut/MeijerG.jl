"""
    _meijerg_perturb(a::Tuple, b::Tuple, m::Integer, n::Integer, z; rtol=nothing, atol=nothing, max_refinements=8)

Compute the Meijer G-function using adaptive parameter perturbation for confluent-pole cases.

**Role in the algorithm**: This is the primary robust fallback when Slater fails or when 
confluence risk is high. It avoids expensive contour integration and scales well to arbitrary precision.

# Method

1. **Perturbation-based regularization**: When active parameters have integer differences,
   the residue expansion has logarithmic singularities. Perturbation removes these by slightly
   separating parameters.
2. **Adaptive convergence loop**: Starts with perturbation size ε, then halves it repeatedly.
   Each level computes 4-node Lagrange interpolation, checking convergence against `rtol`/`atol`.
3. **Termination**: Stops when successive refinements agree to requested tolerance, or after
   `max_refinements` levels (default 8). Converged result is returned without further exceptions.
4. **Type handling**: Computation happens in BigFloat precision for stability; result is converted
   to match input promotion type.

# Keyword arguments
- `rtol`: Relative tolerance for convergence (default: `sqrt(eps(real(z)))` for auto-precision)
- `atol`: Absolute tolerance for convergence (default: 0)
- `max_refinements`: Maximum perturbation halvings (default: 8 levels)
"""
function _meijerg_perturb(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
                          rtol=nothing,
                          atol=nothing,
                          max_refinements::Integer=8)
    _validate_pairing(a, b, m, n)
    
    # We perform interpolation in BigFloat precision
    a_big = big.(a)
    b_big = big.(b)
    z_big = big(z)

    p = length(a)
    q = length(b)
    mode = _expansion_mode(p, q, z)

    target_rtol = rtol === nothing ? sqrt(eps(real(z_big))) : big(rtol)
    target_atol = atol === nothing ? zero(target_rtol) : big(atol)

    eps0 = _confluent_epsilon(a_big, b_big, z_big)
    previous_value = nothing
    result_big = nothing
    converged = false

    for level in 0:max_refinements
        eps_level = eps0 / (big(2)^level)
        value = _meijerg_perturb_fixed_step(a_big, b_big, m, n, z_big, eps_level, mode)

        if previous_value !== nothing
            scale = max(abs(value), abs(previous_value), one(abs(value)))
            if abs(value - previous_value) <= max(target_atol, target_rtol * scale)
                result_big = value
                converged = true
                break
            end
        end

        previous_value = value
        result_big = value
    end

    if result_big === nothing
        throw(ErrorException("Adaptive perturbation did not produce a result."))
    end
    converged || throw(ErrorException("Adaptive perturbation did not converge within max_refinements; increase max_refinements or relax rtol/atol."))

    promoted = promote(z, a..., b...)
    target_type = eltype(promoted)
    return convert(target_type, result_big)
end

function _meijerg_perturb_fixed_step(a::Tuple, b::Tuple, m::Integer, n::Integer, z, eps0, mode::Symbol)
    steps = (-eps0, -eps0 / 2, eps0 / 2, eps0)
    values = map(steps) do eps
        a_perturbed, b_perturbed = _perturb_active_parameters(a, b, m, n, eps, mode)
        meijerg_slater(a_perturbed, b_perturbed, m, n, z)
    end
    return _lagrange_interpolate_zero(steps, values)
end

"""
    _perturb_active_parameters(a::Tuple, b::Tuple, m::Integer, n::Integer, eps, mode::Symbol)

Shift active Meijer G-parameters by small increments proportional to `eps` to remove confluent poles.

For `:lower` mode, displaces `b[1:m]` symmetrically about their center index.
For `:upper` mode, displaces `a[1:n]` symmetrically about their center index.
Returns the modified `(a, b)` pair.
"""
function _perturb_active_parameters(a::Tuple, b::Tuple, m::Integer, n::Integer, eps, mode::Symbol)
    if mode === :lower
        center = (m + 1) / 2
        b_active = ntuple(j -> b[j] + (j - center) * eps, m)
        b_rest = b[m+1:length(b)]
        return a, (b_active..., b_rest...)
    else
        center = (n + 1) / 2
        a_active = ntuple(j -> a[j] + (j - center) * eps, n)
        a_rest = a[n+1:length(a)]
        return (a_active..., a_rest...), b
    end
end

"""
    _confluent_epsilon(a::Tuple, b::Tuple, z) -> Real

Compute a perturbation step size ε scaled to the magnitude of the input parameters.

For confluent limits we cancel singular terms across symmetric perturbations; using
`sqrt(eps)` is often too small at high precision and amplifies cancellation noise.
We therefore use a larger scale, `eps^(1/4)`, which is empirically more stable while
remaining in the asymptotic regime.
"""
function _confluent_epsilon(a::Tuple, b::Tuple, z)
    # Compute a scale-aware perturbation size without downcasting precision.
    magnitudes = map(x -> abs(complex(big(x))), (a..., b..., z))
    scale = maximum((one(first(magnitudes)), magnitudes...))
    return sqrt(sqrt(eps(real(scale)))) * scale
end

"""
    _lagrange_interpolate_zero(xs::NTuple{N}, ys) -> Number

Interpolate the value at x = 0 using two-sided Lagrange polynomial interpolation through `N` points.

Given node abscissas `xs` and function values `ys`, constructs the unique interpolating
polynomial and evaluates it at 0.
"""
function _lagrange_interpolate_zero(xs::NTuple{N}, ys) where {N}
    total = zero(first(ys))
    for i in 1:N
        weight = one(xs[i])
        for j in 1:N
            i == j && continue
            weight *= -xs[j] / (xs[i] - xs[j])
        end
        total += ys[i] * weight
    end
    return total
end