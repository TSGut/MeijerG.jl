# Scalar forward recurrences adapted from ClassicalOrthogonalPolynomials.jl
# v0.15.16 and its RecurrenceRelationships-based evaluation path.
#
# Copyright (c) 2021 JuliaApproximation.
# Licensed under the MIT License.

"""
    _laguerrel(n, α, z)

Evaluate the generalized Laguerre polynomial `L_n^α(z)`
"""
function _laguerrel(n::Integer, α, z)
    n >= 0 || throw(DomainError(n, "Laguerre degree must be nonnegative"))
    T = promote_type(typeof(α), typeof(z))
    αT = convert(T, α)
    zT = convert(T, z)
    p0 = one(T)
    n == 0 && return p0

    p1 = convert(T, muladd(-one(T), zT, αT + one(T)) * p0)
    @inbounds for k in 2:n
        kT = convert(T, k)
        A = -one(T) / kT
        B = (αT + convert(T, 2k - 1)) / kT
        C = (αT + convert(T, k - 1)) / kT
        p1, p0 = muladd(muladd(A, zT, B), p1, -C * p0), p1
    end
    return p1
end

"""
    _jacobip(n, α, β, z)

Evaluate the Jacobi polynomial `P_n^(α,β)(z)`
"""
function _jacobip(n::Integer, α, β, z)
    n >= 0 || throw(DomainError(n, "Jacobi degree must be nonnegative"))
    S = promote_type(typeof(α), typeof(β))
    T = promote_type(S, typeof(z))
    αS = convert(S, α)
    βS = convert(S, β)
    zT = convert(T, z)
    p0 = one(T)
    n == 0 && return p0

    αβ = αS + βS
    A1 = (αβ + convert(S, 2)) / convert(S, 2)
    B1 = (αS - βS) / convert(S, 2)
    p1 = convert(T, muladd(A1, zT, B1) * p0)
    @inbounds for k in 2:n
        j = k - 1
        twoj = convert(S, 2j)
        kS = convert(S, k)
        jS = convert(S, j)
        A = (twoj + αβ + one(S)) * (twoj + αβ + convert(S, 2)) /
            (convert(S, 2) * kS * (jS + αβ + one(S)))
        B = (αS^2 - βS^2) * (twoj + αβ + one(S)) /
            (convert(S, 2) * kS * (jS + αβ + one(S)) *
             (twoj + αβ))
        C = (jS + αS) * (jS + βS) *
            (twoj + αβ + convert(S, 2)) /
            (kS * (jS + αβ + one(S)) * (twoj + αβ))
        p1, p0 = muladd(muladd(A, zT, B), p1, -C * p0), p1
    end
    return p1
end
