"""
    _meijerg_perturb(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Compute the Meijer G-function using parameter perturbation for confluent-pole cases.

This private method handles the special case where the active parameters have integer differences,
creating logarithmic singularities in the residue expansion. Instead of direct summation,
we perturb parameters slightly to remove the poles, then extrapolate back to zero perturbation
using Lagrange interpolation.

# Method
1. Select four perturbation magnitudes: ε, ε/2, ε/4, ε/8 (where ε scales with problem magnitude)
2. For each perturbation, evaluate `meijerg_slater` with slightly separated parameters
3. Apply Lagrange extrapolation on perturbation magnitudes to estimate limit as perturbation → 0
4. Convert result back to original type (handles high-precision BigFloat → Float64 roundtrip)

# References
- Slater, L. J. (1966). Generalized Hypergeometric Functions. Cambridge University Press. §6.9
- Springer Handbook of Special Functions (2015), §16.2.1
"""
function _meijerg_perturb(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    if _needs_high_precision_confluent(a, b, z)
        return setprecision(256) do
            a_big = map(big, a)
            b_big = map(big, b)
            z_big = big(z)
            value_big = _meijerg_perturb_impl(a_big, b_big, m, n, z_big)
            _convert_confluent_result(value_big, a, b, z)
        end
    end
    return _meijerg_perturb_impl(a, b, m, n, z)
end

"""
    _meijerg_perturb_impl(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Core perturbation-and-extrapolation loop for confluent-pole cases.

Evaluates `meijerg_slater` at four perturbation magnitudes (ε, ε/2, ε/4, ε/8) and
applies Lagrange extrapolation to estimate the limit as perturbation → 0.
"""
function _meijerg_perturb_impl(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    mode = _expansion_mode(p, q, z)
    _validate_pairing(a, b, m, n)
    
    eps0 = _confluent_epsilon(a, b, z)
    steps = (eps0, eps0 / 2, eps0 / 4, eps0 / 8)
    values = map(steps) do eps
        a_perturbed, b_perturbed = _perturb_active_parameters(a, b, m, n, eps, mode)
        _meijerg_simple(a_perturbed, b_perturbed, m, n, z, mode)
    end
    return _lagrange_extrapolate_zero(steps, values)
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

"""
    _confluent_epsilon(a::Tuple, b::Tuple, z) -> Real

Compute a perturbation step size ε scaled to the magnitude of the input parameters.

The step size is proportional to the square root of machine epsilon at the dominant scale,
so perturbations are neither too large nor too small relative to the parameter values.
"""
function _confluent_epsilon(a::Tuple, b::Tuple, z)
    magnitudes = map(x -> abs(complex(float(x))), (a..., b..., z))
    scale = maximum((one(first(magnitudes)), magnitudes...))
    return sqrt(eps(float(real(scale)))) * scale
end

"""
    _lagrange_extrapolate_zero(xs::NTuple{N}, ys) -> Number

Extrapolate the limit at x = 0 using Lagrange polynomial interpolation through `N` points.

Given node abscissas `xs` and function values `ys`, constructs the unique interpolating
polynomial and evaluates it at 0.
"""
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

"""
    _needs_high_precision_confluent(a::Tuple, b::Tuple, z) -> Bool

Return `true` if the perturbation evaluation should be promoted to `BigFloat` precision.

Returns `false` if any parameter already contains a `BigFloat` value,
indicating that precision has already been elevated by the caller.
"""
function _needs_high_precision_confluent(a::Tuple, b::Tuple, z)
    _contains_bigfloat(a) && return false
    _contains_bigfloat(b) && return false
    _contains_bigfloat((z,)) && return false
    return true
end

"""
    _contains_bigfloat(values::Tuple) -> Bool

Return `true` if any element of `values` is a `BigFloat` or `Complex{BigFloat}`.
"""
_contains_bigfloat(values::Tuple) = any(_isbigfloatlike, values)

"""
    _isbigfloatlike(x) -> Bool

Return `true` if `x` is a `BigFloat` or `Complex{BigFloat}`, `false` otherwise.
"""
_isbigfloatlike(x::BigFloat) = true
_isbigfloatlike(x::Complex{BigFloat}) = true
_isbigfloatlike(::Any) = false

"""
    _convert_confluent_result(value, a::Tuple, b::Tuple, z) -> Number

Convert a `BigFloat` or `Complex{BigFloat}` result back to the type of the original inputs.

Discards the imaginary part if the promoted type of the original inputs is real.
"""
function _convert_confluent_result(value, a::Tuple, b::Tuple, z)
    a_float = map(float, a)
    b_float = map(float, b)
    target_type = typeof(promote(float(z), a_float..., b_float...)[1])
    if target_type <: Real && value isa Complex
        return convert(target_type, real(value))
    end
    return convert(target_type, value)
end
