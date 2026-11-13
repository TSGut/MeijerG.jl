# Certified Slater dispatch: the separated expansion with the summation-conditioning certificate
# and its double-double lane. Off the convergent regime Float64 summation loses log10(κ_sum)
# digits, so the lane re-evaluates every family as a double-double Maclaurin series.

"""
    _slater_extended(a, b, m, n, z, branch; maxterms)

Double-double Slater expansion for Float64 inputs.  Returns
`(status, value::ComplexF64, kappa)`. A family prefactor at a gamma pole
(unresolved confluence) makes the result non-finite, which is reported as
`:nonfinite`.
"""
function _slater_extended(a::Tuple, b::Tuple, m::Integer, n::Integer, z, branch::Symbol;
                          maxterms::Integer=10_000)
    A = map(_CDD, a)
    B = map(_CDD, b)
    Z = _CDD(z)
    logZ = _dd_log(Z)
    total = zero(_CDD)
    absum = 0.0
    for k in 1:(branch === :lower ? m : n)
        v, va, ok = branch === :lower ?
            _dd_lower_family_term(A, B, m, n, Z, logZ, k; maxterms) :
            _dd_upper_family_term(A, B, m, n, Z, logZ, k; maxterms)
        ok || return (status=:terms_exceeded, value=ComplexF64(NaN, NaN), kappa=Inf)
        total = total + v
        absum += va
    end
    value = ComplexF64(total)
    (_logratio_finite(value) && isfinite(absum)) ||
        return (status=:nonfinite, value=value, kappa=Inf)
    kappa = iszero(value) ? (iszero(absum) ? 1.0 : Inf) : absum / abs(value)
    return (status=:ok, value=value, kappa=kappa)
end

"""
    _try_meijerg_slater_certified(a, b, m, n, z; working_precision=:auto)

Dispatch step of `meijerg` for the separated Slater expansion.  Returns
`(status, value, kappa, lane)`. The value is accepted (`status === :ok`) when
its summation conditioning is certified:

  * in the input precision when `κ_sum <= _LOGRATIO_JET_MAX_SUMMATION_KAPPA`
  * with `working_precision=:auto`, Float64 inputs off the convergent regime, or
    with a larger `κ_sum`, in the double-double lane (`_slater_extended`) while
    its own `κ_sum <= _LOGRATIO_JET_EXTENDED_MAX_KAPPA`. Float32 inputs use the
    Float64 path and BigFloat inputs rerun with raised precision.

Otherwise the status is `:summation_cancellation` (or the failure status) and
`meijerg` continues to contour integration.
"""
function _try_meijerg_slater_certified(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
                                       working_precision::Symbol=:auto)
    extended = working_precision !== :input
    p = length(a)
    q = length(b)
    R = typeof(float(real(z)))
    branch = _expansion_mode(p, q, z)
    if extended && R === Float32
        r = _try_meijerg_slater_certified(map(_logratio_widen, a), map(_logratio_widen, b),
                                          m, n, _logratio_widen(z); working_precision)
        return merge(r, (value=_output_like(r.value, a, b, z),))
    end
    lane = extended && R === Float64
    function extended_lane()
        # An unresolved confluence is singular in any precision. Report `:nonfinite` (or `:failed`
        # on an exception) so that `meijerg` continues to contour integration.
        ext = try
            _slater_extended(a, b, m, n, z, branch)
        catch
            nothing
        end
        ext === nothing &&
            return (status=:failed, value=_logratio_nan_value(z), kappa=Inf, lane=:extended)
        status = ext.status === :ok ?
            (ext.kappa <= _LOGRATIO_JET_EXTENDED_MAX_KAPPA ? :ok : :summation_cancellation) :
            ext.status
        return (status=status, value=ext.value, kappa=ext.kappa, lane=:extended)
    end
    lane && p != q && !_logratio_in_regime(branch, z) && return extended_lane()
    native = try
        _meijerg_slater_certified(a, b, m, n, z)
    catch
        nothing
    end
    native === nothing &&
        return (status=:failed, value=_logratio_nan_value(z), kappa=Inf, lane=:native)
    isfinite(native.value) ||
        return (status=:nonfinite, value=native.value, kappa=Inf, lane=:native)
    _logratio_summation_ok(native.kappa) &&
        return (status=:ok, value=native.value, kappa=native.kappa, lane=:native)
    lane && return extended_lane()
    if extended && R === BigFloat
        hi = _logratio_escalated((a, b, z), z, native.kappa) do aa, bb, zz
            (value=_meijerg_slater_certified(aa, bb, m, n, zz).value,)
        end
        hi === nothing ||
            return (status=:ok, value=hi.value, kappa=native.kappa, lane=:escalated)
    end
    return (status=:summation_cancellation, value=native.value, kappa=native.kappa,
            lane=:native)
end
