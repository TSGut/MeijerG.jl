"""
    _erf(w)

Error function `erf(w)`, evaluated through the regularized lower incomplete
gamma function via `erf(w) = P(1/2, w²)`.

Routing through `Gamma.gamma_inc` keeps `SpecialFunctions` out of the runtime
dependencies and supports real, complex and `BigFloat` arguments uniformly.
"""
function _erf(w)
    wf = float(w)
    R = typeof(real(wf))
    isfinite(wf) || return _erf_nonfinite(wf)
    abs(wf) < sqrt(floatmin(R)) && return 2 * wf / sqrt(convert(R, pi))
    real(wf) < zero(R) && return -_erf(-wf)
    return first(gamma_inc(one(wf) / 2, wf * wf))
end

"""
    _erfc(w)

Complementary error function `erfc(w) = 1 - erf(w)`, evaluated through the
regularized upper incomplete gamma function via `erfc(w) = Q(1/2, w²)`.
"""
function _erfc(w)
    wf = float(w)
    R = typeof(real(wf))
    isfinite(wf) || return one(wf) - _erf_nonfinite(wf)
    real(wf) < zero(R) && return 2 * one(wf) - _erfc(-wf)
    abs(wf) < sqrt(floatmin(R)) && return one(wf) - 2 * wf / sqrt(convert(R, pi))
    return last(gamma_inc(one(wf) / 2, wf * wf))
end

"""
    _erf_nonfinite(w)

Limiting values of `erf` at the non-finite arguments that `Gamma.gamma_inc` has
no entry point for.
"""
_erf_nonfinite(w::Real) = isnan(w) ? w : oftype(w, sign(w))

function _erf_nonfinite(w)
    re, im = reim(w)
    isinf(re) && isfinite(im) && return complex(oftype(re, sign(re)), zero(im))
    isfinite(re) && isinf(im) && return complex(zero(re), im)
    return complex(oftype(re, NaN), oftype(im, NaN))
end
