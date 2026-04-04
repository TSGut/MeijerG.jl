# You will need to have the PyCall package installed
# and use something like conda calls to install mpmath.
using PyCall

# Uncomment the following line to install mpmath via Conda if you haven't already.
# pyimport_conda("mpmath", "mpmath")

mp = pyimport("mpmath")
builtins = pyimport("builtins")
mp.mp.dps = 100 # Set mpmath precision to 100 decimal places (edit as needed)

"""
    pyseq(xs)

Convert a Julia iterable to a Python tuple for mpmath calls.
"""
pyseq(xs) = builtins.tuple(collect(xs))

"""
    pystr(x)

Convert a Python object to its Python `str(...)` representation.
"""
pystr(x) = pycall(builtins.str, String, x)

"""
    _to_mpmath_number(x)

Convert a Julia numeric value to an mpmath numeric object, preserving precision
for BigFloat/Complex{BigFloat} inputs by routing through strings.
"""
function _to_mpmath_number(x)
    if x isa Complex
        return mp.mpc(string(real(x)), string(imag(x)))
    end
    return mp.mpf(string(x))
end

"""
    _to_mpmath_tuple(xs::Tuple)

Convert a Julia tuple of numbers into a Python tuple of mpmath numbers.
"""
_to_mpmath_tuple(xs::Tuple) = pyseq(map(_to_mpmath_number, xs))

"""
    _from_mpmath_number(x, T::Type)

Convert an mpmath result `x` back to Julia type `T`.
"""
function _from_mpmath_number(x, ::Type{T}) where {T}
    re = parse(BigFloat, pystr(mp.re(x)))
    im = parse(BigFloat, pystr(mp.im(x)))
    if T <: Complex
        return convert(T, Complex{BigFloat}(re, im))
    end
    return convert(T, re)
end

"""
    meijerg_mpmath(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Evaluate `mpmath.meijerg` using Julia's full-parameter Meijer G convention
`(a, b, m, n, z)`.
"""
function meijerg_mpmath(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    target_type = eltype(promote(z, a..., b...))
    a_groups = (_to_mpmath_tuple(a[1:n]), _to_mpmath_tuple(a[n+1:end]))
    b_groups = (_to_mpmath_tuple(b[1:m]), _to_mpmath_tuple(b[m+1:end]))
    result = mp.meijerg(a_groups, b_groups, _to_mpmath_number(z))
    return _from_mpmath_number(result, target_type)
end

"""
    meijerg_mpmath(a_left::Tuple, a_right::Tuple, b_left::Tuple, b_right::Tuple, z)

Evaluate `mpmath.meijerg` using the split-parameter convention
`(a_left, a_right, b_left, b_right, z)`, matching the alternate API style
used by `meijerg`.
"""
function meijerg_mpmath(a_left::Tuple, a_right::Tuple, b_left::Tuple, b_right::Tuple, z)
    a = (a_left..., a_right...)
    b = (b_left..., b_right...)
    return meijerg_mpmath(a, b, length(b_left), length(a_left), z)
end

# Test that this works
using Test, MeijerG
@test meijerg((0.2,), (0.1,), 1, 0, 0.5) ≈ meijerg_mpmath((0.2,), (0.1,), 1, 0, 0.5)
@test meijerg((0.2,), (0.1,), 1, 0, 0.5+im*0.1) ≈ meijerg_mpmath((0.2,), (0.1,), 1, 0, 0.5+im*0.1)
@test meijerg((big"0.2",), (big"0.1",), 1, 0, big"0.5") ≈ meijerg_mpmath((big"0.2",), (big"0.1",), 1, 0, big"0.5")
@test meijerg((big"0.2",), (big"0.1",), 1, 0, big"0.5"+im*big"0.1") ≈ meijerg_mpmath((big"0.2",), (big"0.1",), 1, 0, big"0.5"+im*big"0.1")
@test meijerg((0.2,), (), (0.1,), (), 0.5) ≈ meijerg_mpmath((0.2,), (), (0.1,), (), 0.5)
@test meijerg((0.2,), (), (0.1,), (), 0.5 + 0.1im) ≈ meijerg_mpmath((0.2,), (), (0.1,), (), 0.5 + 0.1im)