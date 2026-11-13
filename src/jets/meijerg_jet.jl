@doc raw"""
    meijerg_jet(a, b, m, n, z; working_precision=:auto)
    meijerg_jet(a_left, a_right, b_left, b_right, z; working_precision=:auto)

Evaluate a Meijer G-function with an exact or near-confluent active-pole cluster by
the log-ratio jet residue series.

This is a standalone advanced evaluator, analogous to `meijerg_slater` and
`meijerg_laurent`: it performs no special-function reductions and never falls back to
another algorithm. It selects the exact-pair, finite-separation or higher-order
(three or four poles) jet from the pole geometry. `meijerg` calls the same jets with
`working_precision=:auto` and continues to the Slater expansion or contour integration
when they decline.

# Working precision
Every jet series records its summation conditioning `κ_sum = Σ|T_n| / |Σ T_n|`, the
digits lost to cancellation between its terms. The keyword decides how that is used.

- `:auto` (default, the behaviour of `meijerg`): the series is re-evaluated in doubled
  working precision when it lies off its convergent regime (`|z| > 1` for the lower
  branch, `|z| < 1` for the upper branch) or when `κ_sum` exceeds `1e2`. Float32
  inputs always run on the Float64 path.
- `:input`: the series is summed in the input precision only. Off the convergent
  regime the result can be inaccurate by the factor `κ_sum`, and nothing reports it.
- `:extended`: the series is always summed in doubled working precision. Float32 inputs
  run on the Float64 path, Float64 inputs in an internal double-double kernel, and
  BigFloat inputs at twice the current `precision(BigFloat)` or more when `κ_sum`
  requires it.

The returned type is the promoted type of the inputs in every mode. Off the convergent
regime a sum whose cancellation exceeds even the doubled precision is still returned.

# Errors
Inputs without a cluster the jets support throw an `ArgumentError`. This includes
`p == q`, separated poles, clusters of more than four poles and higher-order clusters
in BigFloat (see `meijerg_laurent` for exact ones). A supported cluster whose series
fails to converge throws an `ErrorException`.

# See Also
- `meijerg`: Main public API with reductions, certified jets and contour fallback
- `meijerg_slater`: Raw separated Slater expansion
- `meijerg_laurent`: Exact-confluence Laurent residue evaluator
"""
function meijerg_jet(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z;
                     working_precision::Symbol=:auto)
    return meijerg_jet(Tuple(a), Tuple(b), m, n, z; working_precision)
end

function meijerg_jet(a_left::ParameterInput, a_right::ParameterInput,
                     b_left::ParameterInput, b_right::ParameterInput, z;
                     working_precision::Symbol=:auto)
    a = (Tuple(a_left)..., Tuple(a_right)...)
    b = (Tuple(b_left)..., Tuple(b_right)...)
    return meijerg_jet(a, b, length(b_left), length(a_left), z; working_precision)
end

function meijerg_jet(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
                     working_precision::Symbol=:auto)
    _check_working_precision(working_precision)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg_jet is implemented for nonzero z only"))
    z isa Real && z < zero(z) &&
        throw(DomainError(z, "meijerg_jet on the negative real axis requires an explicit complex argument. Pass z + 0im."))
    _validate_pairing(a, b, m, n)

    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    result = _try_meijerg_slater_logratio_jet(a_promoted, b_promoted, m, n, z_promoted;
                                              working_precision)
    accepted = (result.status === :ok && result.selected) ||
               result.status === :summation_cancellation
    accepted && _logratio_finite(result.value) &&
        return _output_like(result.value, a, b, z)
    result.status in _LOGRATIO_DECLINED_STATUSES &&
        throw(ErrorException("meijerg_jet: the jet series for this cluster failed ($(result.status))"))
    throw(ArgumentError("meijerg_jet supports only exact or near-confluent active-pole clusters of at most four poles ($(result.status))"))
end
