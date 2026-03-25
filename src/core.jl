const ParameterInput = Union{Tuple, AbstractVector}

_totuple(xs::Tuple) = xs
_totuple(xs::AbstractVector) = Tuple(xs)

@doc raw"""
    meijerg(a, b, m, n, z)
    meijerg(a_left, a_right, b_left, b_right, z)

Compute the Meijer G-function

```math
G_{p,q}^{m,n}\!\left(z\;\middle|\;\begin{matrix} a_1,\ldots,a_p \\ b_1,\ldots,b_q \end{matrix}\right)
```

using Slater's residue expansions (DLMF §16.17). This is the **core evaluator** that always
uses the true mathematical definition, never applying explicit reductions to simpler functions.

# Overview

The Meijer G-function is defined via a Mellin–Barnes contour integral whose residues yield
finite sums of generalized hypergeometric functions (pFq). This implementation:
- Uses the lower expansion (Slater, DLMF 16.17.2) for p < q or p = q, |z| ≤ 1
- Uses the upper expansion (Slater, DLMF 16.17.3) for p > q or p = q, |z| > 1
- Handles confluent poles via parameter perturbation limits
- Supports arbitrary-precision arithmetic via BigFloat and complex numbers

# Arguments
- `a`, `b`: complete upper and lower parameter collections (tuple or vector)
  - `p = length(a)` are the total upper parameters
  - `q = length(b)` are the total lower parameters
- `m` (0 ≤ m ≤ q): number of active lower parameters (contributing residues)
- `n` (0 ≤ n ≤ p): number of active upper parameters
- `z`: evaluation point (Float64, BigFloat, Complex variants, nonzero)

# Split-parameter form
`meijerg(a_left, a_right, b_left, b_right, z)` is equivalent to
```
meijerg((a_left..., a_right...), (b_left..., b_right...),
        length(b_left), length(a_left), z)
```

This follows the standard indexing convention where active top parameters come first,
active bottom parameters come first, etc.

# See Also
- `meijerg_reduce`: Wrapper that attempts reduction to special functions, then falls back to this
"""
function meijerg(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z)
    meijerg(_totuple(a), _totuple(b), m, n, z)
end

function meijerg(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg is implemented for nonzero z only"))

    return _meijerg(a, b, m, n, z)
end

function meijerg(a_left::ParameterInput, a_right::ParameterInput,
                 b_left::ParameterInput, b_right::ParameterInput, z)
    a_left_tuple = _totuple(a_left)
    a_right_tuple = _totuple(a_right)
    b_left_tuple = _totuple(b_left)
    b_right_tuple = _totuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg(a, b, length(b_left_tuple), length(a_left_tuple), z)
end

function _meijerg(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    mode = _expansion_mode(p, q, z)
    _validate_pairing(a, b, m, n)
    _has_confluent_poles(a, b, m, n, mode) && return _confluent_expansion(a, b, m, n, z, mode)
    return _meijerg_simple(a, b, m, n, z, mode)
end

function _meijerg_simple(a::Tuple, b::Tuple, m::Integer, n::Integer, z, mode::Symbol)
    if mode === :lower
        return _lower_expansion(a, b, m, n, z)
    else
        return _upper_expansion(a, b, m, n, z)
    end
end

_expansion_mode(p::Integer, q::Integer, z) = p < q ? :lower : p > q ? :upper : abs(z) <= one(abs(z)) ? :lower : :upper

function _promote_inputs(a::Tuple, b::Tuple, z)
    a_float = map(float, a)
    b_float = map(float, b)
    promoted = promote(float(z), a_float..., b_float...)
    p = length(a)
    q = length(b)
    z_promoted = promoted[1]
    a_promoted = ntuple(i -> promoted[i + 1], p)
    b_promoted = ntuple(i -> promoted[p + i + 1], q)
    return a_promoted, b_promoted, z_promoted
end

function _lower_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    m == 0 && return zero(z)
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(z_promoted, p - m - n)
    total = zero(z_promoted)
    for k in 1:m
        @inbounds begin
            b_k = b_promoted[k]
            α = ntuple(j -> one(b_k) + b_k - a_promoted[j], p)
            β = Tuple(one(b_k) + b_k - b_promoted[j] for j in 1:q if j != k)
            term = _gamma_product((b_promoted[j] - b_k for j in 1:m if j != k))
            term *= _gamma_product((one(b_k) + b_k - a_promoted[j] for j in 1:n))
            term /= _gamma_product((one(b_k) + b_k - b_promoted[j] for j in m+1:q))
            term /= _gamma_product((a_promoted[j] - b_k for j in n+1:p))
            term *= _pow_meijerg(z_promoted, b_k)
            term *= pFq(α, β, argument)
            total += term
        end
    end
    return something(total, zero(z_promoted))
end

function _upper_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    n == 0 && return zero(z)
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(inv(z_promoted), q - m - n)
    total = zero(z_promoted)
    for h in 1:n
        @inbounds begin
            a_h = a_promoted[h]
            α = ntuple(j -> one(a_h) - a_h + b_promoted[j], q)
            β = Tuple(one(a_h) - a_h + a_promoted[j] for j in 1:p if j != h)
            term = _gamma_product((a_h - a_promoted[j] for j in 1:n if j != h))
            term *= _gamma_product((one(a_h) - a_h + b_promoted[j] for j in 1:m))
            term /= _gamma_product((one(a_h) - a_h + a_promoted[j] for j in n+1:p))
            term /= _gamma_product((a_h - b_promoted[j] for j in m+1:q))
            term *= _pow_meijerg(z_promoted, a_h - one(a_h))
            term *= pFq(α, β, argument)
            total += term
        end
    end
    return something(total, zero(z_promoted))
end

_signed_argument(z, exponent::Integer) = isodd(exponent) ? -z : z

function _pow_meijerg(z::Real, exponent)
    if _isintegerlike(exponent)
        return z^Int(round(_realpart(exponent)))
    elseif z >= zero(z)
        return z^exponent
    else
        return complex(z)^exponent
    end
end
_pow_meijerg(z, exponent) = z^exponent

function _gamma_product(values)
    result = nothing
    for value in values
        gamma_value = gamma(value)
        result = result === nothing ? gamma_value : result * gamma_value
    end
    return something(result, 1)
end

function _validate_pairing(a::Tuple, b::Tuple, m::Integer, n::Integer)
    @inbounds for j in 1:n, k in 1:m
        difference = a[j] - b[k]
        _ispositiveintegerlike(difference) && throw(DomainError((a[j], b[k]), "a[j] - b[k] must not be a positive integer for j <= n and k <= m"))
    end
end

_realpart(x::Complex) = real(x)
_realpart(x) = x

_isintegerlike(x::Integer) = true
_isintegerlike(x::Rational) = denominator(x) == 1
_isintegerlike(x::Real) = isfinite(x) && isinteger(x)
_isintegerlike(x::Complex) = iszero(imag(x)) && _isintegerlike(real(x))
_isintegerlike(::Any) = false

function _ispositiveintegerlike(x)
    return _isintegerlike(x) && _realpart(x) > zero(_realpart(x))
end

function _has_confluent_poles(a::Tuple, b::Tuple, m::Integer, n::Integer, mode::Symbol)
    if mode === :lower
        @inbounds for j in 1:m, k in j+1:m
            _isintegerlike(b[j] - b[k]) && return true
        end
    else
        @inbounds for j in 1:n, k in j+1:n
            _isintegerlike(a[j] - a[k]) && return true
        end
    end
    return false
end

function _confluent_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z, mode::Symbol)
    if _needs_high_precision_confluent(a, b, z)
        return setprecision(256) do
            a_big = map(big, a)
            b_big = map(big, b)
            z_big = big(z)
            value_big = _confluent_expansion_impl(a_big, b_big, m, n, z_big, mode)
            _convert_confluent_result(value_big, a, b, z)
        end
    end
    return _confluent_expansion_impl(a, b, m, n, z, mode)
end

function _confluent_expansion_impl(a::Tuple, b::Tuple, m::Integer, n::Integer, z, mode::Symbol)
    eps0 = _confluent_epsilon(a, b, z)
    steps = (eps0, eps0 / 2, eps0 / 4, eps0 / 8)
    values = map(steps) do eps
        a_perturbed, b_perturbed = _perturb_active_parameters(a, b, m, n, eps, mode)
        _meijerg_simple(a_perturbed, b_perturbed, m, n, z, mode)
    end
    return _lagrange_extrapolate_zero(steps, values)
end

function _perturb_active_parameters(a::Tuple, b::Tuple, m::Integer, n::Integer, eps, mode::Symbol)
    if mode === :lower
        b_values = collect(b)
        center = (m + 1) / 2
        for j in 1:m
            b_values[j] = b_values[j] + (j - center) * eps
        end
        return a, Tuple(b_values)
    else
        a_values = collect(a)
        center = (n + 1) / 2
        for j in 1:n
            a_values[j] = a_values[j] + (j - center) * eps
        end
        return Tuple(a_values), b
    end
end

function _confluent_epsilon(a::Tuple, b::Tuple, z)
    magnitudes = map(x -> abs(complex(float(x))), (a..., b..., z))
    scale = maximum((one(first(magnitudes)), magnitudes...))
    return sqrt(eps(float(real(scale)))) * scale
end

function _lagrange_extrapolate_zero(xs::NTuple{N}, ys) where {N}
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

function _needs_high_precision_confluent(a::Tuple, b::Tuple, z)
    _contains_bigfloat(a) && return false
    _contains_bigfloat(b) && return false
    _contains_bigfloat((z,)) && return false
    return true
end

_contains_bigfloat(values::Tuple) = any(_isbigfloatlike, values)
_isbigfloatlike(x::BigFloat) = true
_isbigfloatlike(x::Complex{BigFloat}) = true
_isbigfloatlike(::Any) = false

function _convert_confluent_result(value, a::Tuple, b::Tuple, z)
    a_float = map(float, a)
    b_float = map(float, b)
    target_type = typeof(promote(float(z), a_float..., b_float...)[1])
    if target_type <: Real && value isa Complex
        return convert(target_type, real(value))
    end
    return convert(target_type, value)
end
