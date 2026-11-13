using Test
using MeijerG

finite_value(v) = v isa Number && isfinite(real(v)) && isfinite(imag(v))

function lower_case(; L=1, epsilon=1e-6, z=3e-4 + 0im, m=2)
    a = (0.42 + 0.23im, 1.16 - 0.18im)
    b1 = -0.31 + 0.11im
    btail = (1.37 + 0.22im, 1.91 - 0.17im)
    b = m == 2 ? (b1, b1 + L + epsilon, btail...) :
        (0.82 - 0.13im, b1, b1 + L + epsilon, btail...)
    return (a=a, b=b, m=m, n=1, z=z)
end

function upper_case(; L=1, epsilon=1e-6, z=50.0 + 0im, n=2)
    a1 = -0.22 + 0.18im
    atail = (1.41 - 0.22im, 2.03 + 0.31im)
    a = n == 2 ? (a1, a1 - L + epsilon, atail...) :
        (0.94 + 0.12im, a1, a1 - L + epsilon, atail...)
    b = (0.26 + 0.31im, 1.12 - 0.16im)
    return (a=a, b=b, m=1, n=n, z=z)
end

@testset "Private log-ratio jet near-confluence" begin
    @testset "detection and policy gates" begin
        c = lower_case()
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected
        @test r.branch === :lower

        c = upper_case()
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected
        @test r.branch === :upper

        c = lower_case(epsilon=0.0)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :not_near_confluent
        @test !r.selected

        c = lower_case(epsilon=1e-14)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected

        c = lower_case(epsilon=1e-3)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected

        c = lower_case(epsilon=3e-3)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :not_near_confluent
        @test r.reason === :epsilon_too_large
        @test !r.selected

        c = upper_case(epsilon=3e-3)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :not_near_confluent
        @test r.reason === :epsilon_too_large
        @test !r.selected

        c = lower_case(L=4)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected

        a = (0.2 + 0.1im, 0.7 - 0.2im)
        b = (-0.3 + 0.2im, -0.3 + 0.2im + 1e-6)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(a, b, 2, 1, 3e-4 + 0im)
        @test r.status === :unsupported_finite_p_eq_q

        b1 = -0.31 + 0.11im
        a3 = (0.42 + 0.23im, 1.16 - 0.18im)
        b3 = (b1, b1 + 1e-6, b1 + 2e-6, 1.37 + 0.22im)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(a3, b3, 3, 1, 3e-4 + 0im)
        @test r.status === :unsupported_higher_multiplicity

        # A near pair outside the canonical positions (1, 2), or in the opposite
        # order, is canonicalized by the detector like an exact pair.
        c = lower_case(m=3)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected
        bperm = (c.b[2], c.b[3], c.b[1], c.b[4:end]...)
        rperm = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, bperm, c.m, c.n, c.z)
        @test rperm.status === :ok
        @test rperm.value ≈ r.value rtol=1e-14
        bswap = (c.b[1], c.b[3], c.b[2], c.b[4:end]...)
        rswap = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, bswap, c.m, c.n, c.z)
        @test rswap.status === :ok
        @test rswap.value ≈ r.value rtol=1e-14
        det = MeijerG._detect_near_logratio_jet_pair(c.a, bswap, c.m, c.n, :lower)
        @test det.status === :ok
        @test det.pair.L == 1
        @test det.b[2] - det.b[1] ≈ 1 + 1e-6 atol=1e-12
    end

    @testset "selected row values and dispatch" begin
        lower_expected = -0.36387217343218486 + 4.722879668995371im
        upper_expected = 0.0011935509858065016 + 0.00029286178110097215im

        c = lower_case()
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.value ≈ lower_expected rtol=1e-12 atol=1e-12
        @test meijerg(c.a, c.b, c.m, c.n, c.z) ≈ r.value rtol=1e-14 atol=1e-14

        c = upper_case()
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.value ≈ upper_expected rtol=1e-12 atol=1e-12
        @test meijerg(c.a, c.b, c.m, c.n, c.z) ≈ r.value rtol=1e-14 atol=1e-14

        c = lower_case(z=3e-4 * cis(pi / 8))
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.value ≈ 0.20405720745417494 + 4.531812115721854im rtol=1e-12 atol=1e-12

        c = upper_case(z=50cis(-pi / 6))
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.value ≈ 0.0008518247073055635 + 0.0010613113814142536im rtol=1e-12 atol=1e-12

        for L in (0, 1, 2, 4)
            c = lower_case(L=L)
            @test MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z).status === :ok
            c = upper_case(L=L)
            @test MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z).status === :ok
        end

        c = lower_case(L=10, epsilon=0.0)
        r = MeijerG._try_meijerg_slater_logratio_jet_exact(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected

        c = upper_case(L=10, epsilon=0.0)
        r = MeijerG._try_meijerg_slater_logratio_jet_exact(c.a, c.b, c.m, c.n, c.z)
        @test r.status === :ok
        @test r.selected

        # A representable nonzero residual remains a finite-separation case,
        # even when it is smaller than the former exact-detection tolerance.
        c = lower_case(epsilon=1e-10)
        near = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
        exact = MeijerG._try_meijerg_slater_logratio_jet_exact(c.a, c.b, c.m, c.n, c.z)
        @test near.status === :ok
        @test near.selected
        @test !exact.selected

        # Independently rounded decimal BigFloats can leave an ulp-scale
        # residual at an intended integer separation.  Snap only that
        # representation-scale residual to the exact path.
        setprecision(256) do
            a = (big"1.10", big"2.10", big"3.00")
            b = (big"0.20",)
            z = big"1.60"
            exact = MeijerG._try_meijerg_slater_logratio_jet_exact(a, b, 1, 2, z)
            near = MeijerG._try_meijerg_slater_logratio_jet_near(a, b, 1, 2, z)
            @test exact.status === :ok
            @test exact.selected
            @test !near.selected
        end
    end

    @testset "adaptive grouped-state refresh" begin
        @test MeijerG._logratio_refresh_every(:lower, 3.0 + 0.2im) == 1
        @test MeijerG._logratio_refresh_every(:lower, 3e-4 + 0im) == MeijerG._LOGRATIO_JET_REFRESH
        @test MeijerG._logratio_refresh_every(:upper, 0.5 + 0im) == 1
        @test MeijerG._logratio_refresh_every(:upper, 50.0 + 0im) == MeijerG._LOGRATIO_JET_REFRESH

        a = (0.20 + 0.25im, 0.78 - 0.15im, 1.31 + 0.19im)
        b1 = 0.52 + 0.14im
        b = (b1, b1, -0.47 + 0.28im, 1.05 - 0.21im, 1.67 + 0.16im)
        z = 3.0 + 0.96im

        r = MeijerG._try_meijerg_slater_logratio_jet_exact(a, b, 3, 2, z;
                                                           return_diagnostics=true,
                                                           working_precision=:input)
        direct = MeijerG._lower_logratio_jet_exact_series(a, b, 3, 2, z, 0;
                                                          refresh_every=1)
        fast = MeijerG._lower_logratio_jet_exact_series(a, b, 3, 2, z, 0;
                                                        refresh_every=MeijerG._LOGRATIO_JET_REFRESH)
        # |z| ≈ 3.15 is off the lower regime: the grouped sum cancels (kappa ≈ 400),
        # so the certificate declines and meijerg falls back to contour.
        @test r.status === :summation_cancellation
        @test !r.selected
        @test 1e2 < r.diagnostics.summation_kappa < 1e3
        # Exact states are refreshed in ratio form, so the recurrence path is used
        # off-regime too, and per-term refresh gives no better result.
        @test r.diagnostics.refresh_every == MeijerG._LOGRATIO_JET_REFRESH
        @test r.value == fast.value
        ref = setprecision(BigFloat, 256) do
            ComplexF64(meijerg(map(Complex{BigFloat}, a), map(Complex{BigFloat}, b), 3, 2,
                               Complex{BigFloat}(z)))
        end
        @test abs(fast.value - ref) <= 1.5 * abs(direct.value - ref)
        vc = meijerg_contour(a, b, 3, 2, z)
        @test abs(vc - ref) < 1e-13 * abs(ref)
        @test abs(vc - ref) < abs(r.value - ref)
        # Default: the off-regime exact pair goes to the double-double lane.
        rx = MeijerG._try_meijerg_slater_logratio_jet_exact(a, b, 3, 2, z;
                                                            return_diagnostics=true)
        @test rx.status === :ok && rx.selected
        @test rx.diagnostics.lane === :extended
        @test abs(rx.value - ref) < 4 * eps() * abs(ref)
        @test meijerg(a, b, 3, 2, z) == rx.value

        c = lower_case(epsilon=0.0)
        r = MeijerG._try_meijerg_slater_logratio_jet_exact(c.a, c.b, c.m, c.n, c.z;
                                                           return_diagnostics=true)
        @test r.status === :ok
        @test r.diagnostics.refresh_every == MeijerG._LOGRATIO_JET_REFRESH

        c = upper_case(epsilon=0.0)
        r = MeijerG._try_meijerg_slater_logratio_jet_exact(c.a, c.b, c.m, c.n, c.z;
                                                           return_diagnostics=true)
        @test r.status === :ok
        @test r.diagnostics.refresh_every == MeijerG._LOGRATIO_JET_REFRESH
    end

    @testset "cancellation-free pair tail certificate" begin
        directions = (0.0, 1e-6, 1e-4im, 7e-5 + 6e-5im)
        for branch in (:lower, :upper), L in (0, 1, 4), epsv in directions
            c = branch === :lower ? lower_case(L=L, epsilon=epsv) :
                upper_case(L=L, epsilon=epsv)
            exact = iszero(epsv)
            for nn in 0:12
                if branch === :lower
                    P, D = exact ?
                        MeijerG._logratio_lower_pair_state_limit(c.a, c.b, c.m, c.n,
                                                                 c.z, L, nn) :
                        MeijerG._logratio_lower_pair_state(c.a, c.b, c.m, c.n,
                                                          c.z, L, epsv, nn)
                    term = exact ?
                        (k -> MeijerG._logratio_lower_pair_term_limit(
                            c.a, c.b, c.m, c.n, c.z, L, k)) :
                        (k -> MeijerG._logratio_lower_pair_term(
                            c.a, c.b, c.m, c.n, c.z, L, epsv, k))
                else
                    P, D = exact ?
                        MeijerG._logratio_upper_pair_state_limit(c.a, c.b, c.m, c.n,
                                                                 c.z, L, nn) :
                        MeijerG._logratio_upper_pair_state(c.a, c.b, c.m, c.n,
                                                          c.z, L, epsv, nn)
                    term = exact ?
                        (k -> MeijerG._logratio_upper_pair_term_limit(
                            c.a, c.b, c.m, c.n, c.z, L, k)) :
                        (k -> MeijerG._logratio_upper_pair_term(
                            c.a, c.b, c.m, c.n, c.z, L, epsv, k))
                end
                cert = MeijerG._logratio_pair_tail_bound(
                    P, D, c.a, c.b, c.z, L, epsv, nn, branch)
                isfinite(cert.tail_bound) || continue
                # Absolute partial-tail summation is stronger than checking
                # the modulus of the corresponding signed/complex partial
                # tail.  Stop before direct closed-form terms underflow to
                # zero times an infinite logarithmic factor.
                observed = 0.0
                for k in (nn + 1):(nn + 30)
                    tk = term(k)
                    finite_value(tk) || break
                    observed += abs(tk)
                end
                @test observed <= cert.tail_bound * (1 + 1e-12)
                @test cert.theta < 1
            end
        end

        for c in (lower_case(), upper_case())
            result = MeijerG._try_meijerg_slater_logratio_jet(
                c.a, c.b, c.m, c.n, c.z; return_diagnostics=true)
            @test result.status === :ok
            @test result.diagnostics.tail_status === :accepted
            @test result.diagnostics.tail_ratio_bound < 1
            @test result.diagnostics.tail_estimate <=
                MeijerG._logratio_series_tol(c.z) *
                max(abs(result.value), one(abs(result.value)))
        end
    end

    @testset "fallback and exact interaction" begin
        c = lower_case()
        old = get(ENV, "MEIJERG_LOGRATIO_JET", nothing)
        ENV["MEIJERG_LOGRATIO_JET"] = "0"
        try
            r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z)
            @test !r.selected
            @test r.reason === :disabled
            @test finite_value(meijerg(c.a, c.b, c.m, c.n, c.z))
        finally
            old === nothing ? delete!(ENV, "MEIJERG_LOGRATIO_JET") : (ENV["MEIJERG_LOGRATIO_JET"] = old)
        end

        # The magnitude of z is not a selection criterion.  An upper case at moderate
        # |z| runs the Float64 series. One below |z| = 1 is off the convergent regime
        # and is evaluated in the double-double lane by default, while strict input
        # precision declines it on the summation certificate.
        c = upper_case(z=2.0 + 0im)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z;
                                                          return_diagnostics=true)
        @test r.selected
        @test !haskey(r.diagnostics, :lane)
        @test finite_value(meijerg(c.a, c.b, c.m, c.n, c.z))

        c = upper_case(z=0.1 + 0im)
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z;
                                                          return_diagnostics=true)
        @test r.selected
        @test r.diagnostics.lane === :extended
        rs = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z;
                                                           return_diagnostics=true,
                                                           working_precision=:input)
        # Strict input precision keeps the Float64 series only while its certificate
        # holds (here κ_sum ≈ 80), otherwise it declines it. Either way the value
        # agrees with the lane to the digits the certificate promises.
        @test rs.status === (rs.diagnostics.summation_kappa <=
                             MeijerG._LOGRATIO_JET_MAX_SUMMATION_KAPPA ? :ok :
                             :summation_cancellation)
        vstrict = meijerg_jet(c.a, c.b, c.m, c.n, c.z; working_precision=:input)
        @test finite_value(vstrict)
        @test meijerg(c.a, c.b, c.m, c.n, c.z) ≈ vstrict rtol=1e-11

        c = lower_case()
        r = MeijerG._try_meijerg_slater_logratio_jet_near(c.a, c.b, c.m, c.n, c.z; maxterms=1)
        @test r.status in (:terms_exceeded, :tail_not_accepted)
        @test finite_value(meijerg(c.a, c.b, c.m, c.n, c.z))

        exact = lower_case(epsilon=0.0)
        unified_logratio = MeijerG._try_meijerg_slater_logratio_jet(exact.a, exact.b, exact.m, exact.n, exact.z)
        exact_logratio = MeijerG._try_meijerg_slater_logratio_jet_exact(exact.a, exact.b, exact.m, exact.n, exact.z)
        exact_laurent = meijerg_laurent(exact.a, exact.b, exact.m, exact.n, exact.z)
        @test unified_logratio.status === :ok
        @test unified_logratio.selected
        @test exact_logratio.status === :ok
        @test exact_logratio.selected
        @test unified_logratio.value ≈ exact_logratio.value rtol=1e-12 atol=1e-12
        @test finite_value(exact_logratio.value)
        @test exact_logratio.value ≈ exact_laurent rtol=1e-12 atol=1e-12
        @test meijerg(exact.a, exact.b, exact.m, exact.n, exact.z) ≈ exact_logratio.value rtol=1e-12 atol=1e-12

        exact_upper_reversed = upper_case(L=1, epsilon=0.0)
        a = (exact_upper_reversed.a[2], exact_upper_reversed.a[1],
             exact_upper_reversed.a[3:end]...)
        exact_logratio = MeijerG._try_meijerg_slater_logratio_jet_exact(a, exact_upper_reversed.b,
                                                                        exact_upper_reversed.m,
                                                                        exact_upper_reversed.n,
                                                                        exact_upper_reversed.z)
        @test exact_logratio.status === :ok
        @test exact_logratio.selected
        @test meijerg(a, exact_upper_reversed.b, exact_upper_reversed.m,
                      exact_upper_reversed.n, exact_upper_reversed.z) ≈ exact_logratio.value rtol=1e-12 atol=1e-12
    end

    @testset "all-real gamma-product signs" begin
        lower_a = (0.23, 0.81, 1.36)
        lower_b = (0.51, 0.5101, -0.44, 1.09, 1.69)
        lower = MeijerG._try_meijerg_slater_logratio_jet_near(
            lower_a, lower_b, 3, 2, 2.4e-5)
        @test lower.status === :ok
        @test lower.selected
        @test lower.value ≈ meijerg_slater(lower_a, lower_b, 3, 2, 2.4e-5) rtol=1e-11

        upper_a = (0.44, 0.4401, 1.15, -0.32, 1.72)
        upper_b = (0.28, 0.89, 1.46)
        upper = MeijerG._try_meijerg_slater_logratio_jet_near(
            upper_a, upper_b, 2, 3, 34.0)
        @test upper.status === :ok
        @test upper.selected
        @test upper.value ≈ meijerg_slater(upper_a, upper_b, 2, 3, 34.0) rtol=1e-11
    end

    @testset "public/private API boundaries" begin
        @test_throws DomainError meijerg_slater((1.0,), (0.0,), 1, 1, 0.5)
        @test_throws ArgumentError meijerg_laurent((0.3,), (0.2, 1.1), 2, 1, 0.5)
        @test Base.isexported(MeijerG, :meijerg_laurent)
        @test !Base.isexported(MeijerG, :_try_meijerg_slater_logratio_jet_near)
        @test !Base.isexported(MeijerG, :_try_meijerg_slater_logratio_jet_exact)
        @test !Base.isexported(MeijerG, :_lower_logratio_jet_series)
        @test !Base.isexported(MeijerG, :_upper_logratio_jet_series)
    end
end

@testset "Summation-conditioning certificate and extended-precision lane" begin
    # Exact-pair cases off their branch's convergent regime: the grouped series
    # cancels like exp(nu*|z|^(1/nu)).  References: mpmath meijerg at 130 digits.
    cases = [
        ((0.3,), (0.2, 1.2), 2, 1, 20.0,
         big"0.1166995301016219979085603590513498095972758108483352635445102661674967803911732141289087103293602461"),
        ((), (0.2, 1.2, 0.7), 2, 0, 80.0,
         big"-0.002677243745218326127525353237575407726442827858963248195349324412420109013604634215462975032153075858"),
        ((0.2, 1.2), (0.3,), 1, 2, 0.025,
         big"2.984696898186108175763926664478959040349480677922405836058500116064110724183912604756881896378373176"),
    ]
    for (a, b, m, n, z, ref) in cases
        # Strict Float64: the certificate declines the input-precision series.
        r = MeijerG._try_meijerg_slater_logratio_jet_exact(a, b, m, n, z;
                return_diagnostics=true, working_precision=:input)
        @test r.status === :summation_cancellation
        @test !r.selected
        @test r.diagnostics.summation_kappa > MeijerG._LOGRATIO_JET_MAX_SUMMATION_KAPPA
        @test abs(meijerg_contour(a, b, m, n, z) - ref) < 1e-13 * abs(ref)
        # Default: double-double lane, or contour beyond its own kappa limit.
        rx = MeijerG._try_meijerg_slater_logratio_jet_exact(a, b, m, n, z;
                                                            return_diagnostics=true)
        @test rx.diagnostics.lane === :extended
        @test rx.selected == (rx.diagnostics.summation_kappa <=
                              MeijerG._LOGRATIO_JET_EXTENDED_MAX_KAPPA)
        v = meijerg(a, b, m, n, z)
        @test abs(v - ref) < (rx.selected ? 4 * eps() : 1e-13) * abs(ref)
        # BigFloat reruns the certified series at higher precision instead of declining.
        vb = setprecision(BigFloat, 256) do
            meijerg(map(x -> parse(BigFloat, string(x)), a),
                    map(x -> parse(BigFloat, string(x)), b), m, n,
                    parse(BigFloat, string(z)))
        end
        @test precision(vb) == 256
        @test abs(vb - ref) < big"1e-70" * abs(ref)
    end
    # The first two cases are within the lane's range, the third (kappa ~ 1e14) is not.
    # There the jet declines and meijerg returns the contour value, not the separated
    # expansion.
    @test MeijerG._try_meijerg_slater_logratio_jet_exact(cases[1][1:5]...).selected
    @test MeijerG._try_meijerg_slater_logratio_jet_exact(cases[2][1:5]...).selected
    declined = MeijerG._try_meijerg_slater_logratio_jet(cases[3][1:5]...)
    @test declined.status in MeijerG._LOGRATIO_DECLINED_STATUSES && !declined.selected
    @test meijerg(cases[3][1:5]...) == meijerg_contour(cases[3][1:5]...)
    @test !isfinite(meijerg_slater(cases[3][1:5]...))

    # Float32 inputs: doubled working precision is the Float64 path.
    v32 = meijerg((0.3f0,), (0.2f0, 1.2f0), 2, 1, 10.0f0)
    @test v32 isa Float32
    # Same exactly confluent problem in Float64 (1.2f0 - 0.2f0 rounds to 1 in Float32).
    @test abs(v32 - meijerg((Float64(0.3f0),), (Float64(0.2f0), Float64(0.2f0) + 1), 2, 1, 10.0)) <
          4 * eps(Float32) * abs(v32)

    # In-regime exact pair stays on the Float64 jet.
    r = MeijerG._try_meijerg_slater_logratio_jet_exact((0.3,), (0.2, 1.2), 2, 1, 0.5;
                                                       return_diagnostics=true)
    @test r.selected
    @test !haskey(r.diagnostics, :lane)
    @test r.diagnostics.summation_kappa < MeijerG._LOGRATIO_JET_MAX_SUMMATION_KAPPA

    # Higher-order (r = 3) jet: declined on the same certificate in strict input
    # precision, evaluated in the double-double lane by default.
    r = MeijerG._try_meijerg_slater_higher_logratio_jet((0.3,), (0.2, 1.2, 2.2), 3, 1, 5.0;
                                                        working_precision=:input)
    @test r.status === :summation_cancellation
    r = MeijerG._try_meijerg_slater_higher_logratio_jet((0.3,), (0.2, 1.2, 2.2), 3, 1, 5.0;
                                                        return_diagnostics=true)
    @test r.selected
    @test r.diagnostics.lane === :extended
    @test MeijerG._try_meijerg_slater_higher_logratio_jet((0.3,), (0.2, 1.2, 2.2), 3, 1, 0.5).selected
end
