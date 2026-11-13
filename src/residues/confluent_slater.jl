# Exact-confluence Laurent Slater support for lower (`p < q`) and upper (`p > q`) active-pole
# collisions. No near-confluence, multiplicity above two or balanced `p == q` cases.
# Not called by `meijerg`, exposed as the standalone `meijerg_laurent` evaluator.

const _CONFLUENT_SLATER_MAXELL = 32

_confluent_real_float_type(x) = typeof(float(real(x)))

function _confluent_exact_tol(x, y)
    R = _confluent_real_float_type(x + y)
    return sqrt(eps(R)) * max(abs(x), abs(y), one(abs(x)))
end

function _confluent_near_tol(x, y)
    R = _confluent_real_float_type(x + y)
    return eps(R)^(one(R) / 4) * max(abs(x), abs(y), one(abs(x)))
end

"""
    _detect_exact_confluent_clusters(a, b, m, n, branch; maxell)

Build active-pole clusters for the chosen Slater side.

For `branch === :lower`, the active exponents are `b[j] + ell`, `j <= m`.
For `branch === :upper`, they are `1 - a[j] + ell`, `j <= n`.

The returned `clusters` include singleton clusters because the confluent
residue sum still needs the ordinary simple-pole residues. The `has_near`
flag is set when a pair is close but not exact. Production code treats this
as a fallback case rather than applying the exact double-pole limit.
"""
function _detect_exact_confluent_clusters(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                          branch::Symbol; maxell::Integer=_CONFLUENT_SLATER_MAXELL)
    candidates = Any[]
    if branch === :lower
        for j in 1:m, ell in 0:maxell
            push!(candidates, (center=b[j] + ell, member=(j, ell)))
        end
    elseif branch === :upper
        for j in 1:n, ell in 0:maxell
            push!(candidates, (center=one(a[j]) - a[j] + ell, member=(j, ell)))
        end
    else
        throw(ArgumentError("branch must be :lower or :upper"))
    end

    used = falses(length(candidates))
    clusters = Any[]
    has_near = false
    for i in eachindex(candidates)
        used[i] && continue
        center = candidates[i].center
        members = Any[candidates[i].member]
        used[i] = true
        for j in (i + 1):length(candidates)
            used[j] && continue
            dist = abs(candidates[j].center - center)
            if dist <= _confluent_exact_tol(candidates[j].center, center)
                push!(members, candidates[j].member)
                used[j] = true
            elseif dist <= _confluent_near_tol(candidates[j].center, center)
                has_near = true
            end
        end
        push!(clusters, (center=center, members=members, multiplicity=length(members)))
    end
    return (clusters=clusters, has_near=has_near, candidates=candidates)
end

function _confluent_exact_clusters_locally_isolated(detected)
    for cluster in detected.clusters
        cluster.multiplicity == 2 || continue
        for candidate in detected.candidates
            any(==(candidate.member), cluster.members) && continue
            dist = abs(candidate.center - cluster.center)
            if dist <= _confluent_exact_tol(candidate.center, cluster.center)
                return false
            end
        end
    end
    return true
end

function _confluent_exact_precheck(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                   branch::Symbol; maxell::Integer=_CONFLUENT_SLATER_MAXELL)
    active = branch === :lower ? b[1:m] : branch === :upper ? a[1:n] : ()
    length(active) <= 1 && return (has_exact=false, has_near=false)
    has_near = false
    for i in eachindex(active), j in (i + 1):length(active)
        δ = branch === :lower ? active[j] - active[i] : active[i] - active[j]
        L = round(Int, real(δ))
        abs(L) <= maxell || continue
        residual = δ - L
        if abs(residual) <= _confluent_exact_tol(δ, L)
            return (has_exact=true, has_near=has_near)
        elseif abs(residual) <= _confluent_near_tol(δ, L)
            has_near = true
        end
    end
    return (has_exact=false, has_near=has_near)
end

"""
    _confluent_mb_factors(a, b, m, n)

Return the numerator and denominator gamma factors in the Mellin-Barnes
kernel, without the `z^s` factor. Each factor is represented as
`Gamma(offset + slope*s)`.
"""
function _confluent_mb_factors(a::Tuple, b::Tuple, m::Integer, n::Integer)
    numerator = Any[]
    denominator = Any[]
    for j in 1:m
        push!(numerator, (offset=b[j], slope=-one(b[j])))
    end
    for j in 1:n
        push!(numerator, (offset=one(a[j]) - a[j], slope=one(a[j])))
    end
    for j in (m + 1):length(b)
        push!(denominator, (offset=one(b[j]) - b[j], slope=one(b[j])))
    end
    for j in (n + 1):length(a)
        push!(denominator, (offset=a[j], slope=-one(a[j])))
    end
    return (numerator=numerator, denominator=denominator)
end

function _confluent_pole_order(arg)
    R = _confluent_real_float_type(arg)
    abs(imag(arg)) <= sqrt(eps(R)) || return nothing
    r = round(Int, -real(arg))
    r >= 0 && abs(arg + r) <= sqrt(eps(R)) && return r
    return nothing
end

function _confluent_regular_product_and_logder(regular_num, regular_den, s0)
    R0 = one(s0 + 0im)
    L = zero(s0 + 0im)
    for f in regular_num
        arg = f.offset + f.slope * s0
        R0 *= gamma(arg)
        L += f.slope * digamma(arg)
    end
    for f in regular_den
        arg = f.offset + f.slope * s0
        R0 /= gamma(arg)
        L -= f.slope * digamma(arg)
    end
    return R0, L
end

function _confluent_singular_parts(sing)
    slope = sing.factor.slope
    r = sing.pole_order
    R = _confluent_real_float_type(slope)
    C = (-one(R))^r / gamma(convert(R, r + 1))
    H = sum(inv(convert(R, k)) for k in 1:r; init=zero(R))
    return (principal=C / slope, constant=C * (H - convert(R, Base.MathConstants.eulergamma)))
end

"""
    _confluent_analytic_laurent(factors, s0)

Compute the negative-power Laurent coefficients of the gamma-product kernel at
`s0`. Returns a dictionary with keys `-1` and, for double poles, `-2`.
Returns `nothing` if the local pole is unsupported by the internal exact
confluence path, such as denominator poles or multiplicity greater than two.
"""
function _confluent_analytic_laurent(factors, s0)
    singular_num = Any[]
    regular_num = Any[]
    singular_den = Any[]
    regular_den = Any[]

    for f in factors.numerator
        arg = f.offset + f.slope * s0
        pole_order = _confluent_pole_order(arg)
        if pole_order === nothing
            push!(regular_num, f)
        else
            push!(singular_num, (factor=f, pole_order=pole_order))
        end
    end
    for f in factors.denominator
        arg = f.offset + f.slope * s0
        pole_order = _confluent_pole_order(arg)
        if pole_order === nothing
            push!(regular_den, f)
        else
            push!(singular_den, (factor=f, pole_order=pole_order))
        end
    end

    isempty(singular_den) || return nothing
    1 <= length(singular_num) <= 2 || return nothing

    R0, L = _confluent_regular_product_and_logder(regular_num, regular_den, s0)
    coeffs = Dict{Int, Any}()
    if length(singular_num) == 1
        sp = _confluent_singular_parts(singular_num[1])
        coeffs[-1] = sp.principal * R0
        return coeffs
    end

    sp1 = _confluent_singular_parts(singular_num[1])
    sp2 = _confluent_singular_parts(singular_num[2])
    coeffs[-2] = sp1.principal * sp2.principal * R0
    coeffs[-1] = R0 * (sp1.principal * sp2.constant +
                       sp2.principal * sp1.constant +
                       sp1.principal * sp2.principal * L)
    return coeffs
end

"""
    _confluent_residue_logpoly(coeffs, s0, z, branch)

Convert Laurent coefficients of the gamma-product kernel into the local
log-polynomial residue contribution at `z`.
"""
function _confluent_residue_logpoly(coeffs, s0, z, branch::Symbol)
    orientation = branch === :lower ? -1 : 1
    value = zero(z + 0im)
    maxlog = max(0, -minimum(keys(coeffs)) - 1)
    for r in 0:maxlog
        ck = get(coeffs, -1 - r, zero(value))
        value += orientation * ck * z^s0 * log(z)^r / factorial(r)
    end
    return value
end

"""
    _confluent_slater_sum(a, b, m, n, z, branch, clusters)

Evaluate the Slater residue sum using analytic Laurent residues for singleton
and exact double active-pole clusters. Returns `nothing` if any cluster is
outside the supported internal path.
"""
function _confluent_slater_sum(a::Tuple, b::Tuple, m::Integer, n::Integer, z,
                               branch::Symbol, clusters)
    factors = _confluent_mb_factors(a, b, m, n)
    total = zero(z + 0im)
    R = _confluent_real_float_type(z)
    termtol = convert(R, 10) * eps(R)
    smallcount = 0
    # The legacy double detector stores upper centers as `-s`. Sorting by this
    # center keeps both lower and upper Slater streams in their decaying order.
    ordered = sort(collect(clusters); by=c -> (real(c.center), imag(c.center)))
    for cluster in ordered
        cluster.multiplicity <= 2 || return nothing
        s0 = branch === :lower ? cluster.center : -cluster.center
        coeffs = _confluent_analytic_laurent(factors, s0)
        coeffs === nothing && return nothing
        contrib = _confluent_residue_logpoly(coeffs, s0, z, branch)
        if !_logratio_finite(contrib)
            smallcount >= 1 && break
            return nothing
        end
        total += contrib
        if abs(contrib) <= termtol * max(abs(total), floatmin(R))
            smallcount += 1
            smallcount >= 2 && break
        else
            smallcount = 0
        end
    end
    return total
end

"""
    _try_meijerg_slater_confluent_exact(a, b, m, n, z)

Try the internal exact active-pole double-confluence Slater path.

Returns the computed value when all conservative guards pass. Returns `nothing`
for nonconfluent inputs, near-confluent inputs, multiplicity greater than two,
finite balanced `p == q` cases, or any unsupported local Laurent structure.
"""
function _try_meijerg_slater_confluent_exact(a::Tuple, b::Tuple, m::Integer, n::Integer, z;
                                             maxell::Integer=_CONFLUENT_SLATER_MAXELL)
    p = length(a)
    q = length(b)
    p == q && return nothing

    branch = p < q ? :lower : :upper
    precheck = _confluent_exact_precheck(a, b, m, n, branch; maxell)
    precheck.has_exact || return nothing
    detected = _detect_exact_confluent_clusters(a, b, m, n, branch; maxell)
    any(c -> c.multiplicity > 2, detected.clusters) && return nothing
    any(c -> c.multiplicity == 2, detected.clusters) || return nothing
    # Only poles indistinguishable at the exact-cluster scale invalidate the
    # local double-pole Laurent formula. Broader near checks belong to dispatch.
    _confluent_exact_clusters_locally_isolated(detected) || return nothing

    return _confluent_slater_sum(a, b, m, n, z, branch, detected.clusters)
end

"""
    _laurent_exact_multiplicity(a, b, m, n, branch; maxell)

Return the largest exact active-pole cluster multiplicity seen by the
standalone Laurent evaluator.
"""
function _laurent_exact_multiplicity(a::Tuple, b::Tuple, m::Integer, n::Integer,
                                     branch::Symbol; maxell::Integer=_CONFLUENT_SLATER_MAXELL)
    detected = _detect_exact_confluent_clusters(a, b, m, n, branch; maxell)
    maxmult = maximum(c.multiplicity for c in detected.clusters; init=0)
    return (max_multiplicity=maxmult, detected=detected)
end

@doc raw"""
    meijerg_laurent(a, b, m, n, z)
    meijerg_laurent(a_left, a_right, b_left, b_right, z)

Evaluate exact active-pole confluence by the Laurent/log-polynomial
residue formula.

This is a standalone advanced evaluator, analogous in spirit to
`meijerg_slater`: it performs no special-function reductions and no contour
fallback. It is intended for exact lower or upper active-pole
collisions where the separated Slater terms are singular but the combined
Mellin-Barnes residue is finite.

# Supported cases
- lower branch (`p < q`) active-lower exact collisions
- upper branch (`p > q`) active-upper exact collisions
- active-pole cluster multiplicity at most four.

Unsupported inputs throw an `ArgumentError` rather than falling back to another
algorithm. The primary `meijerg` dispatch uses log-ratio jets for exact and
near confluence. This Laurent evaluator is kept available as an explicit
residue-form diagnostic.

# See Also
- `meijerg`: Main public API with reductions and guarded confluence handling
- `meijerg_slater`: Raw separated Slater expansion
"""
function meijerg_laurent(a::ParameterInput, b::ParameterInput, m::Integer, n::Integer, z)
    meijerg_laurent(Tuple(a), Tuple(b), m, n, z)
end

function meijerg_laurent(a_left::ParameterInput, a_right::ParameterInput,
                         b_left::ParameterInput, b_right::ParameterInput, z)
    a_left_tuple = Tuple(a_left)
    a_right_tuple = Tuple(a_right)
    b_left_tuple = Tuple(b_left)
    b_right_tuple = Tuple(b_right)
    a = (a_left_tuple..., a_right_tuple...)
    b = (b_left_tuple..., b_right_tuple...)
    return meijerg_laurent(a, b, length(b_left_tuple), length(a_left_tuple), z)
end

function meijerg_laurent(a::Tuple, b::Tuple, m::Integer, n::Integer, z)
    p = length(a)
    q = length(b)
    0 <= m <= q || throw(ArgumentError("m must satisfy 0 <= m <= length(b)"))
    0 <= n <= p || throw(ArgumentError("n must satisfy 0 <= n <= length(a)"))
    iszero(z) && throw(DomainError(z, "meijerg_laurent is implemented for nonzero z only"))
    _validate_pairing(a, b, m, n)

    a_promoted, b_promoted, z_promoted = _promote_inputs(a, b, z)
    branch = p < q ? :lower : :upper
    mult = _laurent_exact_multiplicity(a_promoted, b_promoted, m, n, branch)

    if mult.max_multiplicity == 2
        value = _try_meijerg_slater_confluent_exact(a_promoted, b_promoted, m, n, z_promoted)
        value !== nothing && _logratio_finite(value) && return _output_like(value, a, b, z)
    elseif 3 <= mult.max_multiplicity <= _HIGHER_LAURENT_MAXMULT
        result = _try_meijerg_slater_higher_laurent_exact(a_promoted, b_promoted,
                                                          m, n, z_promoted)
        # Standalone evaluator: an ill-conditioned (cancelling) sum is still returned.
        result.status in (:ok, :summation_cancellation) && _logratio_finite(result.value) &&
            return _output_like(result.value, a, b, z)
    end

    throw(ArgumentError("meijerg_laurent supports only exact active-pole confluence of multiplicity at most four"))
end
