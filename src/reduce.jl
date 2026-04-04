"""
    _has_confluent_poles(a::Tuple, b::Tuple, m::Integer, n::Integer, mode::Symbol)

Detect if the Meijer G-function input exhibits confluent poles.

Confluent poles arise when active parameters have integer differences, creating logarithmic
singularities in the residue expansion. Detection depends on expansion mode:
- **Lower expansion** (`mode === :lower`): check integer differences in active lower parameters `b[1:m]`
    (this includes `p < q` and also `p = q` with `|z| <= 1`)
- **Upper expansion** (`mode === :upper`): check integer differences in active upper parameters `a[1:n]`
    (this includes `p > q` and also `p = q` with `|z| > 1`)

# Returns
`true` if confluent poles detected, `false` otherwise.
"""
function _has_confluent_poles(a::Tuple, b::Tuple, m::Integer, n::Integer, mode::Symbol)
    if mode === :lower
        @inbounds for j in 1:m, k in j+1:m
            _isintegerlike(b[j] - b[k]) && return true
        end
    else # mode === :upper
        @inbounds for j in 1:n, k in j+1:n
            _isintegerlike(a[j] - a[k]) && return true
        end
    end
    return false
end

"""
    _reduce_special_case(a::Tuple, b::Tuple, m::Integer, n::Integer, z)

Attempt to recognize special cases and reduce to simpler functions. 

For square-root-based explicit reductions, real inputs are only reduced on
strictly positive real `z`; pass `z + 0im` to evaluate on the complex branch.

Returns `nothing` if no rule matches (caller will fall back to full `meijerg` evaluation).

This is the core dispatcher for the reduction system. It checks for reduction rules in order:

0. **Order cancellation**: Parameters appearing in both `a_left` and `b_right` (or vice versa)
   cancel out, reducing problem size recursively.
1. **Exponential**: `G_{0,1}^{1,0}(z | - ; 0) = exp(-z)`
2. **Sine**: `G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2√z) / √π`
3. **Cosine**: `G_{0,2}^{1,0}(z | - ; 0, 1/2) = cos(2√z) / √π`
4. **Bessel I continuation at order 0**:
    `G_{0,2}^{1,0}(z | - ; 0, 0) = I_0(2√(-z))` for negative real `z`
5. **Bessel J**: `G_{0,2}^{1,0}(z | - ; ν/2, -ν/2) = J_ν(2√z)`
6. **Bessel K (symmetric and shifted form)**:
    `G_{0,2}^{2,0}(z | - ; b₁, b₂) = z^{-c} 2K_ν(2√z)` with
    `c = (b₁+b₂)/2`, `ν = |b₁-b₂|`
7. **Erf half-order reduction**:
    `G_{1,2}^{1,1}(z | 1 ; 1/2, 0) = √π erf(√z)`
8. **Erfc half-order reduction**:
    `G_{1,2}^{2,1}(z | 1 ; 1/2, 0) = √π erfc(√z)`
9. **Lower incomplete gamma**: `G_{1,2}^{1,1}(z | 1 ; a, 0) = γ(a,z)`
10. **Upper incomplete gamma**: `G_{1,2}^{2,1}(z | 1 ; a, 0) = Γ(a,z)`
11. **Laguerre polynomial**:
     `G_{1,2}^{1,0}(z | n+1 ; 0, -α) = L_n^(α)(z)/Γ(n+α+1)` for integer `n ≥ 0`
12. **Hermite polynomial channels** via Laguerre specializations:
     even branch `b=(0, 1/2)` and odd branch `b=(0, -1/2)`
13. **Jacobi polynomial**:
     `P_n^(α,β)(x) = Γ(n+α+1)Γ(-n-α-β) G_{2,2}^{1,0}((x-1)/2 | n+1, -n-α-β ; 0, -α)`
14. **Logarithm (1+z)**: `G_{2,2}^{1,2}(z | 1, 1 ; 1, 0) = log(1+z) / z` (direct evaluation via `log1p`)
"""
function _reduce_special_case(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    # Cancellation/order reduction: eliminate exact parameter matches across
    # left/right partitions, then continue reducing recursively.
    reduced_orders = _reduce_orders(a, b, m, n)
    if reduced_orders !== nothing
        a_reduced, b_reduced, m_reduced, n_reduced = reduced_orders
        return meijerg(a_reduced, b_reduced, m_reduced, n_reduced, z)
    end

    half_zero_tuple = length(b) == 2 && begin
        half_b1 = one(b[1]) / (one(b[1]) + one(b[1]))
        b[1] == half_b1 && iszero(b[2])
    end
    zero_half_tuple = length(b) == 2 && begin
        half_b2 = one(b[2]) / (one(b[2]) + one(b[2]))
        iszero(b[1]) && b[2] == half_b2
    end
    zero_neg_half_tuple = length(b) == 2 && begin
        half_b2 = one(b[2]) / (one(b[2]) + one(b[2]))
        iszero(b[1]) && b[2] == -half_b2
    end

    # G_{0,1}^{1,0}(z | - ; 0) = exp(-z)
    if a === () && m == 1 && n == 0 && isequal(b, (0,))
        z_promoted = float(z)
        return exp(-z_promoted)
    end

    # G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2sqrt(z))/sqrt(pi)
    # Restrict to positive real z to keep real-branch semantics explicit.
    if a === () && m == 1 && n == 0 && half_zero_tuple && z isa Real && z <= zero(z)
        throw(DomainError(z, "Sine reduction requires positive real z; pass z + 0im for complex-branch evaluation."))
    end
    if _can_use_sqrt_reduction(z) && a === () && m == 1 && n == 0 && half_zero_tuple
        z_promoted = float(z)
        return sin(2 * sqrt(z_promoted)) * inv(sqrt(convert(typeof(z_promoted), pi)))
    end

    # G_{0,2}^{1,0}(z | - ; 0, 1/2) = cos(2sqrt(z))/sqrt(pi)
    # Restrict to positive real z to keep real-branch semantics explicit.
    if a === () && m == 1 && n == 0 && zero_half_tuple && z isa Real && z <= zero(z)
        throw(DomainError(z, "Cosine reduction requires positive real z; pass z + 0im for complex-branch evaluation."))
    end
    if _can_use_sqrt_reduction(z) && a === () && m == 1 && n == 0 && zero_half_tuple
        z_promoted = float(z)
        return cos(2 * sqrt(z_promoted)) * inv(sqrt(convert(typeof(z_promoted), pi)))
    end

    # G_{0,2}^{1,0}(z | - ; 0, 0) = J_0(2sqrt(z)) = I_0(2sqrt(-z))
    # For negative-real arguments this explicit I_0 form avoids avoidable branch-noise
    # and returns the expected real-valued result for the principal continuation.
    if a === () && m == 1 && n == 0 && length(b) == 2 && iszero(b[1]) && iszero(b[2])
        z_promoted = float(z)
        if z_promoted isa Real && z_promoted < 0
            return besseli(zero(z_promoted), 2 * sqrt(-z_promoted))
        end
    end

    # G_{0,2}^{1,0}(z | - ; ν/2, -ν/2) = J_ν(2sqrt(z))
    # Restrict to positive real z to avoid branch ambiguities for sqrt on real negatives.
    if a === () && m == 1 && n == 0 && length(b) == 2 && isequal(b[2], -b[1]) && z isa Real && z <= zero(z)
        throw(DomainError(z, "Bessel-J reduction requires positive real z; pass z + 0im for complex-branch evaluation"))
    end
    if _can_use_sqrt_reduction(z) && a === () && m == 1 && n == 0 && length(b) == 2 && isequal(b[2], -b[1])
        z_promoted = float(z)
        ν = float(b[1] - b[2])
        try
            return besselj(ν, 2 * sqrt(z_promoted))
        catch err
            # SpecialFunctions currently lacks the BigFloat non-integer-order
            # input path required here, so fall back to Slater for BigFloat
            # (real or complex) until direct Bessel-J support is available.
            if err isa MethodError && (z isa BigFloat || z isa Complex{BigFloat})
                return nothing
            end
            rethrow()
        end
    end

    # G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2sqrt(z))
    # G_{0,2}^{2,0}(z | - ; b1, b2) = z^{-c} * 2K_ν(2√z)  where c = (b1+b2)/2, ν = |b1-b2|
    # Derived from DLMF 16.19.1 parameter-shift formula + DLMF 10.32.10 Bessel K identity.
    # Handles any b1, b2 including the confluent-pole case b1-b2 ∈ ℤ.
    if a === () && m == 2 && n == 0 && length(b) == 2 && z isa Real && z <= zero(z)
        throw(DomainError(z, "Bessel-K reduction requires positive real z; pass z + 0im for complex-branch evaluation"))
    end
    if _can_use_sqrt_reduction(z) && a === () && m == 2 && n == 0 && length(b) == 2
        z_promoted = float(z)
        b1f, b2f = float(b[1]), float(b[2])
        c  = (b1f + b2f) / (one(b1f) + one(b1f))
        ν  = abs(float(b1f - b2f))
        return z_promoted^(-c) * 2 * besselk(ν, 2 * sqrt(z_promoted))
    end

    # G_{1,2}^{1,1}(z | 1 ; a, 0) = γ(a,z)
    # For a=1/2 this is sqrt(pi)*erf(sqrt(z)); keep an explicit branch for clarity
    # and better numerical behavior at small |z|.
    if n == 1 && m == 1 && isequal(a, (1,)) && half_zero_tuple && z isa Real && z <= zero(z)
        throw(DomainError(z, "Erf reduction requires positive real z; pass z + 0im for complex-branch evaluation"))
    end
    if _can_use_sqrt_reduction(z) && n == 1 && m == 1 && isequal(a, (1,)) && half_zero_tuple
        z_promoted = float(z)
        return sqrt(convert(typeof(z_promoted), pi)) * erf(sqrt(z_promoted))
    end

    # G_{1,2}^{2,1}(z | 1 ; 1/2, 0) = sqrt(pi)*erfc(sqrt(z))
    if n == 1 && m == 2 && isequal(a, (1,)) && half_zero_tuple && z isa Real && z <= zero(z)
        throw(DomainError(z, "Erfc reduction requires positive real z; pass z + 0im for complex-branch evaluation"))
    end
    if _can_use_sqrt_reduction(z) && n == 1 && m == 2 && isequal(a, (1,)) && half_zero_tuple
        z_promoted = float(z)
        return sqrt(convert(typeof(z_promoted), pi)) * erfc(sqrt(z_promoted))
    end

    # Orthogonal polynomial reductions (integer n >= 0):
    # G_{1,2}^{1,0}(z | n+1 ; 0, -α) = L_n^(α)(z) / Γ(n+α+1)
    # and Hermite-even/odd specializations via α = -1/2, +1/2.
    if n == 0 && m == 1 && length(a) == 1 && length(b) == 2 && iszero(b[1])
        n_candidate = a[1] - one(a[1])
        if _isintegerlike(n_candidate) && real(n_candidate) >= 0
            n_int = Int(round(real(n_candidate)))
            z_promoted = float(z)
            α = -float(b[2])

            # General Laguerre reduction.
            if !zero_half_tuple && !zero_neg_half_tuple
                lag_arg = float(n_int) + α + 1
                inv_gamma = if lag_arg isa Real && lag_arg < zero(lag_arg)
                    logabsγ, signγ = logabsgamma(lag_arg)
                    inv(signγ * exp(logabsγ))
                else
                    exp(-loggamma(lag_arg))
                end
                return laguerrel(n_int, α, z_promoted) * inv_gamma
            end

            # Hermite-even channel: b = (0, 1/2), equivalent to α = -1/2.
            if zero_half_tuple
                coeff = 2.0^(2n_int) * factorial(big(n_int))
                denom = factorial(big(2n_int)) * sqrt(convert(typeof(z_promoted), pi))
                half = one(z_promoted) / (one(z_promoted) + one(z_promoted))
                return (coeff / denom) * laguerrel(n_int, -half, z_promoted)
            end

            # Hermite-odd channel (x factored outside): b = (0, -1/2), α = +1/2.
            if zero_neg_half_tuple
                coeff = 2.0^(2n_int + 1) * factorial(big(n_int))
                denom = factorial(big(2n_int + 1)) * sqrt(convert(typeof(z_promoted), pi))
                half = one(z_promoted) / (one(z_promoted) + one(z_promoted))
                return (coeff / denom) * laguerrel(n_int, half, z_promoted)
            end
        end
    end

    # Jacobi polynomial reduction (integer n >= 0):
    # P_n^(α,β)(x) = Γ(n+α+1)Γ(-n-α-β) * G_{2,2}^{1,0}((x-1)/2 | n+1, -n-α-β ; 0, -α)
    # with x = 2z + 1 for the Meijer-G argument z.
    if n == 0 && m == 1 && length(a) == 2 && length(b) == 2 && iszero(b[1])
        n_candidate = a[1] - one(a[1])
        if _isintegerlike(n_candidate) && real(n_candidate) >= 0
            n_int = Int(round(real(n_candidate)))
            z_promoted = float(z)

            α = -float(b[2])
            β = -float(a[2]) - float(n_int) - α
            two = one(z_promoted) + one(z_promoted)
            x_jacobi = two * z_promoted + one(z_promoted)

            jac_arg1 = float(n_int) + α + 1
            jac_arg2 = -float(n_int) - α - β
            inv_g1 = if jac_arg1 isa Real && jac_arg1 < zero(jac_arg1)
                logabsγ1, signγ1 = logabsgamma(jac_arg1)
                inv(signγ1 * exp(logabsγ1))
            else
                exp(-loggamma(jac_arg1))
            end
            inv_g2 = if jac_arg2 isa Real && jac_arg2 < zero(jac_arg2)
                logabsγ2, signγ2 = logabsgamma(jac_arg2)
                inv(signγ2 * exp(logabsγ2))
            else
                exp(-loggamma(jac_arg2))
            end
            inv_pref = inv_g1 * inv_g2
            return jacobip(n_int, α, β, x_jacobi) * inv_pref
        end
    end

    # G_{1,2}^{1,1}(z | 1 ; a, 0) = γ(a,z)
    if n == 1 && m == 1 && isequal(a, (1,)) && length(b) == 2 && isequal(b[2], zero(b[2]))
        z_promoted = float(z)
        α = float(b[1])
        full_gamma = if α isa Real && α < zero(α)
            logabsγ, signγ = logabsgamma(α)
            signγ * exp(logabsγ)
        else
            exp(loggamma(α))
        end
        return full_gamma - gamma(α, z_promoted)
    end

    # G_{1,2}^{2,1}(z | 1 ; a, 0) = Γ(a,z)
    if n == 1 && m == 2 && isequal(a, (1,)) && length(b) == 2 && isequal(b[2], zero(b[2]))
        z_promoted = float(z)
        α = float(b[1])
        return gamma(α, z_promoted)
    end

    # Logarithmic confluent poles: G_{2,2}^{1,2}(z | 1, 1 ; 1, 0)
    # Maps to log(1+z) / z via parameter perturbation limit (Gradshteyn & Ryzhik 9.353, DLMF 15.8.2)
    # Direct computation avoids expensive 4× evaluation + two-sided Lagrange interpolation
    if length(a) == 2 && length(b) == 2 && m == 1 && n == 2 && isequal(a, (1, 1)) && isequal(b, (1, 0))
        z_promoted = float(z)
        return log1p(z_promoted) / z_promoted
    end

    return nothing
end

"""
    _reduce_orders(a::Tuple, b::Tuple, m::Integer, n::Integer)

Eliminate parameter cancellations to reduce Meijer G-function order.

When a parameter appears in both `a_left = a[1:n]` and `b_right = b[m+1:end]`,
the two cancel out. Similarly for parameters in `a_right = a[n+1:end]` and `b_left = b[1:m]`.

This function repeatedly removes such pairs until no more cancellations remain.

# Returns
- `nothing` if no reductions occur (i.e., input is already reduced)
- `(a_reduced, b_reduced, m_reduced, n_reduced)` tuple with cancelled parameters removed

# Example
```
a = (1.0, 0.25)
b = (0.25, 1.0)
m = 1, n = 1

# a_left = (1.0), a_right = (0.25), b_left = (0.25), b_right = (1.0)
# Match: a_left[1] == b_right[1] (both 1.0)
# Match: a_right[1] == b_left[1] (both 0.25)
# Result: a_reduced=(), b_reduced=(), m_reduced=0, n_reduced=0
```
"""
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

    reduced = (Tuple(vcat(a_left, a_right)), Tuple(vcat(b_left, b_right)), m, n)
    return reduced == (a, b, m, n) ? nothing : reduced
end

"""
    _can_use_sqrt_reduction(z) -> Bool

Return whether sqrt-based real-form reductions are allowed for `z`.

This guard enforces Julia's real-domain convention for explicit reductions
that use `sqrt(z)`: real reductions are only used for strictly positive real `z`,
while complex inputs are always allowed and follow principal-branch continuation.

# Returns
`true` when a sqrt-based reduction path is permitted, otherwise `false`.
"""
_can_use_sqrt_reduction(z) = !(z isa Real) || (z > zero(z))
