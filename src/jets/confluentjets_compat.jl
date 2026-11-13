# Private compatibility names for the extracted ConfluentJets implementation, kept so the
# move does not reorder any floating-point arithmetic.

const _Jet = Jet
const _jet_order = jet_order
const _jet_coeff = coefficient
const _jet_const = constant_jet
const _jet_var = variable_jet
const _jet_inv = inv
const _jet_exp = exp
const _jet_log = log
const _jet_compose_scale = compose_scale

_higher_loggamma_jet(arg0::T, slope::T, N::Integer) where {T<:Number} =
    loggamma(affine_jet(arg0, slope, N))

_higher_loglinear_jet(c::T, slope::T, N::Integer) where {T<:Number} =
    affine_log_jet(c, slope, N)

function _higher_direct_recipgamma_jet(
        arg0::T, slope::T, N::Integer) where {T<:Number}
    return reciprocal_gamma(
        affine_jet(arg0, slope, N);
        radius=_HIGHER_LOGRATIO_DIRECT_RECIPGAMMA_RADIUS)
end

_higher_hsym(p::Integer, vars::AbstractVector) =
    complete_homogeneous(p, vars)
_higher_taylor_matrix(values::AbstractVector, deltas::AbstractVector) =
    taylor_matrix(values, deltas)
_higher_divided_difference_exp(values::AbstractVector,
                               deltas::AbstractVector) =
    exp_divided_difference(Jet(values), deltas)
_higher_logratio_pad(deltas) =
    default_opitz_pad(deltas; maximum_pad=_HIGHER_LOGRATIO_LOCAL_PAD)

const _logratio_finite = finite_number
const _logratio_real_float_type = real_float_type
const _logratio_int_as = integer_as
const _logratio_exprel = exprel
const _logratio_log1p = safe_log1p
const _logratio_tail_h = scaled_exprel
const _logratio_logpoch = logpochhammer
const _logratio_loggamma_divdiff_limit = loggamma_divdiff_limit
const _logratio_logpoch_divdiff_limit = logpochhammer_divdiff_limit

function _logratio_loggamma_divdiff_signed(
        a, epsv, sigma::Integer; tol=nothing)
    return loggamma_divdiff(a, epsv, sigma; tolerance=tol)
end

_logratio_loggamma_divdiff_pair(a, epsv) = loggamma_divdiff_pair(a, epsv)

_logratio_logpoch_divdiff_signed(a, n::Integer, epsv, sigma::Integer) =
    logpochhammer_divdiff(a, n, epsv, sigma)
