@doc raw"""
    meijerg(a, b, m, n, z)
    meijerg(a_left, a_right, b_left, b_right, z)

Compute the Meijer G-function with automatic reduction to simpler functions and robust handling
of confluent poles.

# Overview

This is the **primary public API** and recommended entry point. It implements a hierarchical 
polyalgorithm strategy:

1. **Special-case reduction**: Recognizes patterns (exponential, trigonometric, Bessel, 
   incomplete gamma, etc.) and reduces to elementary or special functions
2. **Parameter validation**: Checks for valid parameter pairings (skipped for special reductions)
3. **Confluence assessment**: Quantifies proximity to confluent poles (integer-separated parameters)
4. **Method selection** (in order):
    - If confluence risk is **low** (< 0.5): Try Slater residue expansion (fast, accurate)
    - Otherwise or if Slater fails: **Try adaptive perturbation regularization** (robust for confluent poles)
      - Adaptively converges to user tolerance without numerical integration
      - Preserves full precision without BigFloat slowdowns
    - If perturbation fails: Emergency contour methods (rare in practice)

# Arguments
- `a`, `b`: complete upper and lower parameter collections (tuple/vector form)
- `m`, `n`: Meijer G indices with `0 <= m <= length(b)` and `0 <= n <= length(a)`
- `z`: evaluation point (any Julia number type: Float64, BigFloat, Complex, etc.)

# Alternative signature
The split-parameter form `meijerg(a_left, a_right, b_left, b_right, z)` is equivalent to
`meijerg((a_left..., a_right...), (b_left..., b_right...), length(b_left), length(a_left), z)`.
This convention follows the standard notation where `n = |a_left|` and `m = |b_left|`.

# Returns
The computed G-function value, with return type matching the promoted input types.
Supports arbitrary precision via `BigFloat` and complex arguments.

# Examples
```julia
# Exponential: G_{0,1}^{1,0}(z | - ; 0) = exp(-z)
meijerg((), (), (0,), (), -0.5)  # ≈ exp(0.5)

# Sine: √π · G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2√z)
sqrt(pi) * meijerg((), (), (0.5,), (0,), 0.1^2/4)

# Generic input: uses robust confluence-aware polyalgorithm
meijerg((0.3, 0.8), (0.2, 1.1), 1, 1, 0.7)
```

# See Also
- `meijerg_slater`: Pure Slater residue expansion (power users, no reductions or perturbation; unsafe for confluent poles)
- `meijerg_contour`: Direct Mellin-Barnes integration (power users only; use explicit `contour_case=` to select method)
"""
function meijerg(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z)
    meijerg(Tuple(a), Tuple(b), m, n, z)
end

function meijerg(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg is implemented for nonzero z only"))
    
    # Step 1: Try special-case reductions (exponential, Bessel, etc.)
    reduced = _reduce_special_case(a, b, m, n, z)
    if reduced !== nothing
        return reduced
    end
    
    # Validate parameter pairing after special reductions; required for standard algorithm.
    _validate_pairing(a, b, m, n)
    
    # Step 2: Promote inputs to a common floating-point type
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)

    # Step 3: Assess confluence risk and select evaluation method
    mode = _expansion_mode(p, q, z)
    confluence_risk = _confluence_risk(a, b, m, n, mode)
    
    # Try Slater if confluence risk is low
    if confluence_risk < 0.5
        try
            return meijerg_slater(a_promoted, b_promoted, m, n, z_promoted)
        catch e
            # Slater failed; continue to methods below
        end
    end
    
    # Try adaptive perturbation regularization before contour methods.
    try
        return _meijerg_perturb(a_promoted, b_promoted, m, n, z_promoted)
    catch e
        # Adaptive perturbation failed; continue to contour methods.
    end

    # Try contour methods in geometry-aware order:
    # - if a vertical strip exists: Case i first, then loop case
    # - otherwise: loop case first, then Case i
    has_strip = _has_vertical_strip(a_promoted, b_promoted, m, n)
    loop_case = mode === :lower ? :ii : :iii

    if has_strip
        try
            return meijerg_contour(a_promoted, b_promoted, m, n, z_promoted; contour_case=:i)
        catch e
            # Case i failed; try loop contour.
        end

        try
            return meijerg_contour(a_promoted, b_promoted, m, n, z_promoted; contour_case=loop_case)
        catch e
            # Loop contour failed after Case i.
        end
    else
        try
            return meijerg_contour(a_promoted, b_promoted, m, n, z_promoted; contour_case=loop_case)
        catch e
            # Loop contour failed; try vertical line.
        end

        try
            return meijerg_contour(a_promoted, b_promoted, m, n, z_promoted; contour_case=:i)
        catch e
            # Case i failed after loop contour.
        end
    end
    
    throw(ErrorException("meijerg: all evaluation methods failed (Slater, adaptive perturbation, loop contour, vertical contour)."))
end

function meijerg(a_left::ParameterInput, a_right::ParameterInput,
                 b_left::ParameterInput, b_right::ParameterInput, z)
    a_left_tuple = Tuple(a_left)
    a_right_tuple = Tuple(a_right)
    b_left_tuple = Tuple(b_left)
    b_right_tuple = Tuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg(a, b, length(b_left_tuple), length(a_left_tuple), z)
end

"""
    _confluence_risk(a::Tuple, b::Tuple, m::Integer, n::Integer, mode::Symbol) -> Real

Heuristic to quantify how close the parameters are to creating confluent poles.

Returns a value in [0, 1] where:
- 0 = no confluence risk, safe to use Slater
- 1 = severe confluence risk, must use robust contour methods

The heuristic measures the minimum distance from actual parameters to nearby integers,
scaled to account for the number of at-risk parameters.

# Algorithm

For lower expansion mode (`mode === :lower`):
  - Check all pairs of active lower parameters b[j], b[k] for j < k
  - For each pair, compute min_distance = min(|frac(b[j] - b[k])|, 1 - |frac(b[j] - b[k])|)
  - where frac(x) = x - floor(x) is the fractional part
  - Average all pairwise distances and invert to get risk score

For upper expansion mode (`mode === :upper`):
  - Same logic applied to active upper parameters a[1:n]

# Risk Scaling

- If minimum distance is > 0.2: risk ≈ 0 (safe)
- If minimum distance is 0.05-0.2: risk gradually increases
- If minimum distance is < 0.05: risk ≈ 1 (dangerous)
"""
function _confluence_risk(a::Tuple, b::Tuple, m::Integer, n::Integer, mode::Symbol)::Real
    parameters = mode === :lower ? b[1:m] : a[1:n]
    num_active = length(parameters)
    
    # Trivial cases: no confluence risk possible
    if num_active <= 1
        return 0.0
    end
    
    # Collect pairwise distances to nearest integer
    min_distances = Float64[]
    @inbounds for j in 1:num_active
        for k in j+1:num_active
            diff = real(parameters[j] - parameters[k])
            frac = abs(diff - round(diff))
            frac = min(frac, 1.0 - frac)  # distance to nearest integer
            push!(min_distances, frac)
        end
    end
    
    if isempty(min_distances)
        return 0.0
    end
    
    # Risk is inverse of minimum distance, scaled to [0, 1]
    # Threshold: distance > 0.2 is safe, distance < 0.05 is dangerous
    min_dist = minimum(min_distances)
    
    if min_dist > 0.2
        return 0.0
    elseif min_dist < 0.05
        return 1.0
    else
        # Linear interpolation between 0.05 (risk=1) and 0.2 (risk=0)
        return (0.2 - min_dist) / 0.15
    end
end

