@doc raw"""
    meijerg(a, b, m, n, z)
    meijerg(a_left, a_right, b_left, b_right, z)

Compute the Meijer G-function with automatic reduction to simpler functions and guarded
handling of confluent active-pole clusters.

# Overview

This is the **primary public API** and recommended entry point. It implements a hierarchical 
polyalgorithm strategy:

1. **Special-case reduction**: Recognizes patterns (exponential, trigonometric, Bessel, 
   incomplete gamma, etc.) and reduces to elementary or special functions
2. **Parameter validation**: Checks for valid parameter pairings (skipped for special reductions)
3. **Log-ratio jet selection**: Uses guarded exact and near-confluent active-pole
   log-ratio jets when their support policy selects
4. **Raw Slater expansion**: Uses the separated-pole residue expansion otherwise
5. **Contour fallback**: Falls back to Mellin-Barnes contour integration if Slater fails
   or if a residue sum cannot be certified

Every residue sum (jet series and separated Slater expansion) records its summation
conditioning `κ_sum = Σ|T_n| / |Σ T_n|`, the digits lost to cancellation between
its terms. This certificate decides between the input-precision result, an internal
doubled-precision re-evaluation, and contour integration.  A sum is re-evaluated
when it lies off the series' convergent regime (`|z| > 1` for the lower branch,
`|z| < 1` for the upper branch) or when `κ_sum` exceeds `1e2`: Float32 inputs on
the Float64 path, Float64 inputs in an internal double-double kernel (about five
times the cost of the Float64 kernel for the exact-pair series, more for the
separated Slater expansion and the finite-separation and higher-order jets), and
BigFloat inputs with more bits.  Contour integration remains the fallback beyond
`κ_sum ≈ 1e14`.  The returned type is unchanged by this.

# Arguments
- `a`, `b`: complete upper and lower parameter collections (tuple/vector form)
- `m`, `n`: Meijer G indices with `0 <= m <= length(b)` and `0 <= n <= length(a)`
- `z`: evaluation point (any Julia number type: Float64, BigFloat, Complex, etc.)

# Alternative signature
The split-parameter form `meijerg(a_left, a_right, b_left, b_right, z)` is equivalent to
`meijerg((a_left..., a_right...), (b_left..., b_right...), length(b_left), length(a_left), z)`.
This convention follows the standard notation where `n = |a_left|` and `m = |b_left|`.

# Returns
The computed G-function value in the promoted floating-point type of `a`, `b` and `z`:
real (e.g. `Float64`, `BigFloat`) when all inputs are real, complex otherwise, whichever
evaluation path is taken. Supports arbitrary
precision via `BigFloat` and complex arguments.

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
- `meijerg_slater`: Pure Slater residue expansion (power users, no reductions, unsafe for confluent poles)
- `meijerg_laurent`: Standalone exact-confluence Laurent/log-polynomial residue evaluator
- `meijerg_contour`: Direct Mellin-Barnes integration (power users only, use explicit `contour_case=` to select method)
"""
function meijerg(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z)
    meijerg(Tuple(a), Tuple(b), m, n, z)
end

meijerg(a::Tuple, b::Tuple, m::Integer, n::Integer, z) =
    _output_like(_meijerg(a, b, m, n, z), a, b, z)

function _meijerg(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
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
    
    # Validate parameter pairing after special reductions, required for the standard algorithm.
    _validate_pairing(a, b, m, n)
    
    # Step 2: Promote inputs to a common floating-point type
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)

    # Step 3: Select an evaluation method
    mode = _expansion_mode(p, q, z)

    # Try the log-ratio jets for exact and near-confluent active-pole clusters. The Laurent
    # evaluator stays out of this chain (`meijerg_laurent`).
    logratio_jet = _try_meijerg_slater_logratio_jet(a_promoted, b_promoted, m, n, z_promoted)
    if logratio_jet.status === :ok && logratio_jet.selected && isfinite(logratio_jet.value)
        return logratio_jet.value
    end
    # A detected but declined confluent geometry must not fall back to the separated
    # expansion.
    confluent_declined = logratio_jet.status in _LOGRATIO_DECLINED_STATUSES

    # The higher-order (r = 3, 4) jet is Float32/Float64-only. For other precisions use
    # the type-generic exact Laurent residue evaluator on exact clusters of that size.
    if !(typeof(float(real(z_promoted))) <: Union{Float32,Float64})
        laurent = _try_meijerg_slater_higher_laurent_exact(a_promoted, b_promoted, m, n, z_promoted)
        if laurent.status === :ok && laurent.max_multiplicity >= 3 && isfinite(laurent.value)
            return laurent.value
        end
    end

    # Otherwise use the ordinary Slater expansion with its summation certificate. Unsupported
    # confluent geometries (e.g. p == q) fall through when it throws or cannot certify its sum.
    if !confluent_declined
        slater = _try_meijerg_slater_certified(a_promoted, b_promoted, m, n, z_promoted)
        slater.status === :ok && isfinite(slater.value) && return slater.value
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
            # Case i failed, try loop contour.
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
            # Loop contour failed, try vertical line.
        end

        try
            return meijerg_contour(a_promoted, b_promoted, m, n, z_promoted; contour_case=:i)
        catch e
            # Case i failed after loop contour.
        end
    end
    
    z isa Real && z < zero(z) && throw(DomainError(z, "meijerg on the negative real axis requires an explicit complex argument. Pass z + 0im."))
    throw(ErrorException("meijerg: all evaluation methods failed (log-ratio jet, Slater, loop contour, vertical contour)."))
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
