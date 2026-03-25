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

using Slater's residue expansions. This implementation supports arbitrary Julia
number types, including `BigFloat` and complex numbers built from them, by
forwarding the hypergeometric pieces to `HypergeometricFunctions.jl` and the
Gamma factors to `SpecialFunctions.jl`.

Arguments:
- `a`, `b`: complete upper and lower parameter collections
- `m`, `n`: Meijer G indices with `0 <= m <= length(b)` and `0 <= n <= length(a)`
- `z`: evaluation point

The split-parameter form follows the common convention
`meijerg(a_left, a_right, b_left, b_right, z)`, where
`n == length(a_left)` and `m == length(b_left)`.

Current limitation: only simple-pole cases are implemented. If repeated poles
would contribute logarithmic terms, a `DomainError` is thrown.
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

    a_reduced, b_reduced, m_reduced, n_reduced = _reduce_orders(a, b, m, n)
    return _meijerg(a_reduced, b_reduced, m_reduced, n_reduced, z)
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
    if mode === :lower
        _validate_simple_lower(a, b, m, n)
        return _lower_expansion(a, b, m, n, z)
    else
        _validate_simple_upper(a, b, m, n)
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
    # Parallel evaluate each residue-term independently, then reduce serially.
    # Allocate a concrete-typed buffer for terms to avoid type-instability
    # and excessive allocations. Determine an appropriate element type by
    # promoting the runtime types of the inputs (z and the parameters).
    types = (typeof(z_promoted),)
    if p > 0
        types = (types..., map(typeof, a_promoted)...)
    end
    if q > 0
        types = (types..., map(typeof, b_promoted)...)
    end
    term_eltype = length(types) > 1 ? promote_type(types...) : types[1]
    terms = Vector{term_eltype}(undef, m)
    Threads.@threads for k in 1:m
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
            terms[k] = term
        end
    end
    total = zero(z_promoted)
    for k in 1:m
        total += terms[k]
    end
    return something(total, zero(z_promoted))
end

function _upper_expansion(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    n == 0 && return zero(z)
    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    p = length(a_promoted)
    q = length(b_promoted)
    argument = _signed_argument(inv(z_promoted), q - m - n)
    # Allocate a concrete-typed buffer for terms to avoid type-instability
    # and excessive allocations. Determine an appropriate element type by
    # promoting the runtime types of the inputs (z and the parameters).
    types = (typeof(z_promoted),)
    if p > 0
        types = (types..., map(typeof, a_promoted)...)
    end
    if q > 0
        types = (types..., map(typeof, b_promoted)...)
    end
    term_eltype = length(types) > 1 ? promote_type(types...) : types[1]
    terms = Vector{term_eltype}(undef, n)
    # Parallel evaluate each residue-term independently, then reduce serially.
    Threads.@threads for h in 1:n
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
            terms[h] = term
        end
    end
    total = zero(z_promoted)
    for h in 1:n
        total += terms[h]
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

function _reduce_orders(a::Tuple, b::Tuple, m::Integer, n::Integer)
    a_left = collect(a[1:n])
    a_right = collect(a[n+1:end])
    b_left = collect(b[1:m])
    b_right = collect(b[m+1:end])

    changed = true
    while changed
        changed = false

        for i in eachindex(a_left)
            j = findfirst(x -> isequal(x, a_left[i]), b_right)
            if j !== nothing
                deleteat!(a_left, i)
                deleteat!(b_right, j)
                n -= 1
                changed = true
                break
            end
        end
        changed && continue

        for i in eachindex(a_right)
            j = findfirst(x -> isequal(x, a_right[i]), b_left)
            if j !== nothing
                deleteat!(a_right, i)
                deleteat!(b_left, j)
                m -= 1
                changed = true
                break
            end
        end
    end

    return (Tuple(vcat(a_left, a_right)), Tuple(vcat(b_left, b_right)), m, n)
end

function _validate_simple_lower(a::Tuple, b::Tuple, m::Integer, n::Integer)
    _validate_pairing(a, b, m, n)
    @inbounds for j in 1:m, k in j+1:m
        _isintegerlike(b[j] - b[k]) && throw(DomainError((b[j], b[k]), "simple-pole lower expansion requires distinct bottom parameters among the first m entries"))
    end
end

function _validate_simple_upper(a::Tuple, b::Tuple, m::Integer, n::Integer)
    _validate_pairing(a, b, m, n)
    @inbounds for j in 1:n, k in j+1:n
        _isintegerlike(a[j] - a[k]) && throw(DomainError((a[j], a[k]), "simple-pole upper expansion requires distinct top parameters among the first n entries"))
    end
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
