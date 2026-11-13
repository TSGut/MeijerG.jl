using HypergeometricFunctions: pFq
import Bessels
using ClassicalOrthogonalPolynomials: laguerrel, hermiteh, jacobip
using SpecialFunctions: gamma, besselj, besselk, besseli, erf, erfc

function big_slater_reference(a, b, m, n, z)
    return setprecision(BigFloat, 512) do
        a_big = map(BigFloat, a)
        b_big = map(BigFloat, b)
        meijerg_slater(a_big, b_big, m, n, BigFloat(z))
    end
end

@testset "Supported explicit reductions" begin
    @testset "Order cancellation" begin
        for z in sample_real_points((0.15, 0.6, 1.2); seed=101, lo=0.1, hi=1.4)
            reduced = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
            reference = meijerg((), (0.25,), 1, 0, z)
            @test reduced ≈ reference rtol=1e-14
        end

        for z in sample_complex_points((0.25 + 0.15im, 0.7 - 0.2im);
                                       seed=102,
                                       real_lo=0.1,
                                       real_hi=1.0,
                                       imag_lo=-0.4,
                                       imag_hi=0.4)
            reduced = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
            reference = meijerg((), (0.25,), 1, 0, z)
            @test reduced ≈ reference rtol=1e-14
        end
    end

    @testset "Elementary and Bessel reductions" begin
        @testset "Exponential" begin
            for x in sample_real_points((0.05, 0.3, 0.9); seed=103, lo=0.02, hi=1.0)
                value = meijerg((), (), (0.0,), (), -x)
                @test value isa Float64
                @test value ≈ exp(x) rtol=1e-14
            end

            for xb in BigFloat.(sample_real_points((0.05, 0.3, 0.9); seed=104, lo=0.02, hi=1.0))
                valueb = meijerg((), (), (BigFloat(0),), (), -xb)
                @test valueb isa BigFloat
                @test valueb ≈ exp(xb) atol=big"1e-30" rtol=big"1e-30"
            end

            for zc in sample_complex_points((0.2 + 0.3im, 0.7 + 0.1im);
                                            seed=105,
                                            real_lo=0.1,
                                            real_hi=1.0,
                                            imag_lo=-0.4,
                                            imag_hi=0.4)
                valuec = meijerg((), (), (0.0,), (), zc)
                @test valuec ≈ exp(-zc) rtol=1e-14
            end
        end

        @testset "Sine / Cosine" begin
            for y in sample_real_points((0.1, 0.7, 1.4); seed=106, lo=0.05, hi=1.6)
                z = y^2 / 4
                @test sqrt(pi) * meijerg((), (), (0.5,), (0.0,), z) ≈ sin(y) rtol=1e-14
                @test sqrt(pi) * meijerg((), (), (0.0,), (0.5,), z) ≈ cos(y) rtol=1e-14
            end

            for yb in BigFloat.(sample_real_points((0.1, 0.7, 1.4); seed=107, lo=0.05, hi=1.6))
                zb = yb^2 / 4
                @test sqrt(big(pi)) * meijerg((), (), (big"0.5",), (big"0.0",), zb) ≈ sin(yb) atol=big"1e-28" rtol=big"1e-28"
                @test sqrt(big(pi)) * meijerg((), (), (big"0.0",), (big"0.5",), zb) ≈ cos(yb) atol=big"1e-28" rtol=big"1e-28"
            end
        end

        @testset "Sinh / Cosh analytic continuation" begin
            for x in sample_real_points((0.1, 0.5, 1.2); seed=108, lo=0.05, hi=2.0)
                z = -x^2 / 4 + 0im
                lhs_sinh = -im * sqrt(pi) * meijerg((), (), (0.5,), (0.0,), z)
                lhs_cosh = sqrt(pi) * meijerg((), (), (0.0,), (0.5,), z)
                @test real(lhs_sinh) ≈ sinh(x) atol=1e-12 rtol=1e-12
                @test real(lhs_cosh) ≈ cosh(x) atol=1e-12 rtol=1e-12
            end
        end

        @testset "Bessel J" begin
            ν = 0.8
            for z in sample_real_points((0.15, 1.2, 2.4); seed=109, lo=0.1, hi=2.8)
                lhs = meijerg((), (), (ν / 2,), (-ν / 2,), z)
                reference = big_slater_reference(
                    (), (ν / 2, -ν / 2), 1, 0, z)
                @test lhs ≈ reference rtol=big"2e-14"
            end
            for zb in BigFloat.(sample_real_points((0.15, 1.2, 2.4); seed=110, lo=0.1, hi=2.8))
                # Bessels.jl 0.2.8 has no BigFloat general-order method, so
                # this reduction must fall back to Slater.
                @test_throws MethodError Bessels.besselj(ν, 2 * sqrt(zb))
                @test meijerg((), (), (ν / 2,), (-ν / 2,), zb) ≈
                      meijerg_slater((), (), (ν / 2,), (-ν / 2,), zb) atol=big"1e-28" rtol=big"1e-28"
            end
            for zc in sample_complex_points((0.25 + 0.1im, 0.9 - 0.15im);
                                            seed=111,
                                            real_lo=0.1,
                                            real_hi=1.5,
                                            imag_lo=-0.5,
                                            imag_hi=0.5)
                lhsc = meijerg((), (), (ν / 2,), (-ν / 2,), zc)
                rhsc = besselj(ν, 2 * sqrt(zc))
                @test lhsc ≈ rhsc atol=1e-12 rtol=1e-12
            end
        end

        @testset "Bessel K" begin
            ν = 0.8
            for z in sample_real_points((0.15, 1.2, 2.4); seed=112, lo=0.1, hi=2.8)
                lhs = meijerg((), (), (ν / 2, -ν / 2), (), z)
                reference = big_slater_reference(
                    (), (ν / 2, -ν / 2), 2, 0, z)
                @test lhs ≈ reference rtol=big"2e-14"

                lhs_shifted = meijerg((), (), (0.3, -0.5), (), z)
                shifted_reference = setprecision(BigFloat, 512) do
                    b1 = BigFloat(0.3)
                    b2 = BigFloat(-0.5)
                    BigFloat(z)^((b1 + b2) / 2) * reference
                end
                @test lhs_shifted ≈ shifted_reference rtol=big"2e-14"
            end
            for zc in sample_complex_points((0.25 + 0.1im, 0.9 - 0.15im);
                                            seed=114,
                                            real_lo=0.1,
                                            real_hi=1.5,
                                            imag_lo=-0.5,
                                            imag_hi=0.5)
                lhsc = meijerg((), (), (ν / 2, -ν / 2), (), zc)
                rhsc = 2 * besselk(ν, 2 * sqrt(zc))
                @test lhsc ≈ rhsc atol=1e-12 rtol=1e-12
            end
        end

        @testset "Bessel I (nu=0 continuation)" begin
            ν = 0.0
            for x in sample_real_points((0.05, 0.7, 2.0); seed=115, lo=0.02, hi=3.0)
                z = -x^2 / 4
                value = meijerg((), (), (ν / 2,), (-ν / 2,), z)
                reference = big_slater_reference(
                    (), (ν / 2, -ν / 2), 1, 0, z)
                @test value ≈ reference rtol=big"2e-14"
            end
        end

        @testset "Airy Ai via Bessel K mapping" begin
            for x in sample_real_points((0.1, 0.9, 2.4); seed=117, lo=0.02, hi=3.0)
                z = x^3 / 9
                g_val = meijerg((), (), (1 / 6, -1 / 6), (), z)
                ai_api = (1 / (2 * pi)) * sqrt(x / 3) * g_val
                g_reference = big_slater_reference(
                    (), (1 / 6, -1 / 6), 2, 0, z)
                ai_reference = (inv(BigFloat(2) * big(pi)) *
                                sqrt(BigFloat(x) / 3) * g_reference)
                @test ai_api ≈ ai_reference rtol=big"2e-14"
            end
        end
    end

    @testset "Error and incomplete-gamma reductions" begin
        # The in-package erf/erfc go through the regularized incomplete gamma
        # (P, Q) at a = 1/2. SpecialFunctions is the independent reference.
        @testset "Internal _erf / _erfc vs SpecialFunctions" begin
            for x in sample_real_points((-3.5, -0.4, 0.4, 3.5, 8.0); seed=126, lo=-12.0, hi=12.0)
                @test MeijerG._erf(x) ≈ erf(x) rtol=1e-13
                @test MeijerG._erfc(x) ≈ erfc(x) rtol=1e-12
            end

            # Full complex plane, including the left half-plane that the
            # Meijer-G reductions themselves never reach (they feed sqrt(z)).
            for re in -3.0:1.0:3.0, im in -3.0:1.0:3.0
                w = complex(re, im)
                @test MeijerG._erf(w) ≈ erf(w) atol=1e-13 rtol=1e-13
                @test MeijerG._erfc(w) ≈ erfc(w) atol=1e-13 rtol=1e-13
            end

            # Tail accuracy: erfc must stay relatively accurate where erf saturates.
            for x in (5.0, 10.0, 20.0, 26.0)
                @test MeijerG._erfc(x) ≈ erfc(x) rtol=1e-12
            end

            # Small-argument branch, where w^2 underflows and 2w/sqrt(pi) takes over.
            for x in (0.0, -0.0, 1e-30, -1e-30, 1e-200, -1e-300)
                @test MeijerG._erf(x) ≈ erf(x) atol=0 rtol=1e-14
                @test MeijerG._erfc(x) ≈ erfc(x) rtol=1e-14
            end

            @test MeijerG._erf(Inf) == 1.0
            @test MeijerG._erf(-Inf) == -1.0
            @test MeijerG._erfc(Inf) == 0.0
            @test MeijerG._erfc(-Inf) == 2.0

            # Type propagation and extended precision.
            @test MeijerG._erf(1) isa Float64
            @test MeijerG._erf(big"0.5") isa BigFloat
            @test MeijerG._erf(0.5 + 0.25im) isa ComplexF64
            @test Float64(MeijerG._erf(big"0.5")) ≈ erf(0.5) rtol=1e-15
            for xb in (big"0.25", big"1.5", big"4.0", big"-2.25")
                @test MeijerG._erf(xb) + MeijerG._erfc(xb) ≈ one(xb) atol=big"1e-70"
            end
        end

        @testset "erf / erfc half-order reductions" begin
            for x in sample_real_points((0.1, 0.6, 1.5); seed=118, lo=0.05, hi=2.0)
                z = x^2
                erf_api = meijerg((1.0,), (0.5, 0.0), 1, 1, z) / sqrt(pi)
                erfc_api = meijerg((1.0,), (0.5, 0.0), 2, 0, z) / sqrt(pi)
                @test erf_api ≈ erf(x) atol=1e-12 rtol=1e-12
                @test erfc_api ≈ erfc(x) atol=1e-12 rtol=1e-12
            end

            for zb in BigFloat.(sample_real_points((0.1, 0.6, 1.5); seed=119, lo=0.05, hi=2.0)).^2
                erf_big = meijerg((BigFloat(1),), (big"0.5", big"0.0"), 1, 1, zb) / sqrt(big(pi))
                erfc_big = meijerg((BigFloat(1),), (big"0.5", big"0.0"), 2, 0, zb) / sqrt(big(pi))
                @test erf_big ≈ erf(sqrt(zb)) atol=big"1e-28" rtol=big"1e-28"
                @test erfc_big ≈ erfc(sqrt(zb)) atol=big"1e-28" rtol=big"1e-28"
            end

            for zc in sample_complex_points((0.1 + 0.1im, 0.5 + 0.2im, 1.0 - 0.15im);
                                            seed=120,
                                            real_lo=0.05,
                                            real_hi=1.2,
                                            imag_lo=-0.3,
                                            imag_hi=0.3)
                erf_complex = meijerg((1.0,), (0.5, 0.0), 1, 1, zc) / sqrt(pi)
                erfc_complex = meijerg((1.0,), (0.5, 0.0), 2, 0, zc) / sqrt(pi)
                @test erf_complex ≈ erf(sqrt(zc)) atol=1e-12 rtol=1e-12
                @test erfc_complex ≈ erfc(sqrt(zc)) atol=1e-12 rtol=1e-12
            end
        end

        @testset "Lower / upper incomplete gamma" begin
            a = 0.7
            for z in sample_real_points((0.2, 1.3, 2.8); seed=121, lo=0.1, hi=3.0)
                lower_result = meijerg((1.0,), (a, 0.0), 1, 1, z)
                upper_result = meijerg((1.0,), (a, 0.0), 2, 0, z)
                # SpecialFunctions' gamma(a, z) itself is off the 256-bit value by up to
                # ~4e-14 at some z.
                @test lower_result ≈ gamma(a) - gamma(a, z) rtol=1e-13
                @test upper_result ≈ gamma(a, z) rtol=1e-13
            end
            # G^{2,1}_{1,2}(z | 1 ; a, 0) is undefined (a_1 - b_2 = 1): no reduction, rejected
            @test_throws DomainError meijerg((1.0,), (a, 0.0), 2, 1, 0.8)
            @test_throws DomainError meijerg((1.0,), (0.5, 0.0), 2, 1, 0.8)

            ab = big"0.7"
            for zb in BigFloat.(sample_real_points((0.2, 1.3, 2.8); seed=122, lo=0.1, hi=3.0))
                lower_big = meijerg((BigFloat(1),), (ab, big"0.0"), 1, 1, zb)
                upper_big = meijerg((BigFloat(1),), (ab, big"0.0"), 2, 0, zb)
                @test lower_big ≈ gamma(ab) - gamma(ab, zb) atol=big"1e-28" rtol=big"1e-28"
                @test upper_big ≈ gamma(ab, zb) atol=big"1e-28" rtol=big"1e-28"
            end

            for zc in sample_complex_points((0.2 + 0.15im, 0.8 - 0.2im, 1.4 + 0.1im);
                                            seed=123,
                                            real_lo=0.1,
                                            real_hi=1.6,
                                            imag_lo=-0.4,
                                            imag_hi=0.4)
                lower_complex = meijerg((1.0,), (a, 0.0), 1, 1, zc)
                upper_complex = meijerg((1.0,), (a, 0.0), 2, 0, zc)
                @test lower_complex ≈ gamma(a) - gamma(a, zc) atol=1e-12 rtol=1e-12
                @test upper_complex ≈ gamma(a, zc) atol=1e-12 rtol=1e-12
            end
        end
    end

    @testset "Orthogonal-polynomial reductions" begin
        @testset "Internal scalar recurrence parity" begin
            laguerre_cases = (
                (0, 0.2, -0.7),
                (1, 0.7, 0.3),
                (6, 1.3, 2.1),
                (5, big"0.7", big"-0.4"),
            )
            for (degree, α, x) in laguerre_cases
                @test isequal(
                    MeijerG._laguerrel(degree, α, x),
                    laguerrel(degree, α, x))
            end

            jacobi_cases = (
                (0, 0.35, 0.6, -0.7),
                (1, 0.35, 0.6, 0.3),
                (6, 0.35, 0.6, 0.8),
                (5, big"0.35", big"0.6", big"-0.4"),
            )
            for (degree, α, β, x) in jacobi_cases
                @test isequal(
                    MeijerG._jacobip(degree, α, β, x),
                    jacobip(degree, α, β, x))
            end

            # ClassicalOrthogonalPolynomials does not accept complex scalar
            # arguments here, so check that extension independently through
            # the terminating hypergeometric definitions.
            n = 6
            α = 0.35
            β = 0.6
            z = 0.2 + 0.3im
            laguerre_reference =
                gamma(n + α + 1) / (gamma(α + 1) * gamma(n + 1)) *
                pFq((-n,), (α + 1,), z)
            jacobi_reference =
                gamma(n + α + 1) / (gamma(α + 1) * gamma(n + 1)) *
                pFq((-n, n + α + β + 1), (α + 1,), (1 - z) / 2)
            @test MeijerG._laguerrel(n, α, z) ≈ laguerre_reference rtol=2e-14
            @test MeijerG._jacobip(n, α, β, z) ≈ jacobi_reference rtol=2e-14
        end

        @testset "Generalized Laguerre" begin
            x_values = sample_real_points((-1.4, -0.3, 0.4, 1.2, 2.3); seed=124, lo=-1.5, hi=2.6)
            for α in (0.2, 0.7, 1.3)
                for n in 0:6
                    for x in x_values
                        lhs = gamma(big(n) + α + 1) * meijerg((n + 1.0,), (0.0, -α), 1, 0, x)
                        rhs = laguerrel(n, α, x)
                        @test lhs ≈ rhs atol=1e-14 rtol=1e-14
                    end
                end
            end

            x_values_big = BigFloat.(sample_real_points((-1.1, -0.2, 0.3, 1.0, 2.0); seed=125, lo=-1.2, hi=2.2))
            for αb in (big"0.2", big"0.7")
                for n in 0:5
                    for xb in x_values_big
                        lhsb = gamma(big(n) + αb + 1) * meijerg((big(n) + 1,), (big"0.0", -αb), 1, 0, xb)
                        rhsb = laguerrel(n, αb, xb)
                        @test lhsb ≈ rhsb atol=big"1e-28" rtol=big"1e-28"
                    end
                end
            end
        end

        @testset "Hermite" begin
            x_values = sample_real_points((-1.4, -0.8, -0.2, 0.2, 0.8, 1.4); seed=127, lo=-1.6, hi=1.6)
            for n in 0:5
                for x in x_values
                    g_even = meijerg((n + 1.0,), (0.0, 0.5), 1, 0, x^2)
                    h_even = (-1)^n * factorial(big(2n)) * sqrt(big(pi)) * g_even
                    # The Float64 reduction loses about a digit for n >= 3, |x| > 1.
                    @test h_even ≈ hermiteh(2n, x) atol=1e-13 rtol=1e-13

                    g_odd = meijerg((n + 1.0,), (0.0, -0.5), 1, 0, x^2)
                    h_odd = (-1)^n * factorial(big(2n + 1)) * x * sqrt(big(pi)) * g_odd
                    @test h_odd ≈ hermiteh(2n + 1, x) atol=1e-13 rtol=1e-14
                end
            end

            x_values_big = BigFloat.(sample_real_points((-1.2, -0.6, -0.1, 0.1, 0.6, 1.2); seed=128, lo=-1.4, hi=1.4))
            for n in 0:4
                for xb in x_values_big
                    g_even_b = meijerg((big(n) + 1,), (big"0.0", big"0.5"), 1, 0, xb^2)
                    h_even_b = (-1)^n * factorial(big(2n)) * sqrt(big(pi)) * g_even_b
                    @test h_even_b ≈ hermiteh(2n, xb) atol=big"1e-28" rtol=big"1e-28"

                    g_odd_b = meijerg((big(n) + 1,), (big"0.0", -big"0.5"), 1, 0, xb^2)
                    h_odd_b = (-1)^n * factorial(big(2n + 1)) * xb * sqrt(big(pi)) * g_odd_b
                    @test h_odd_b ≈ hermiteh(2n + 1, xb) atol=big"1e-27" rtol=big"1e-27"
                end
            end
        end

        @testset "Jacobi" begin
            α, β = 0.35, 0.6

            @testset "Float64 (+/- sample inputs)" begin
                for n in 0:6
                    for x in sample_real_points((-0.8, -0.3, 0.3, 0.8); seed=130 + n, lo=-0.85, hi=0.85)
                        z = (x - 1) / 2
                        lhs = meijerg((n + 1.0, -n - α - β), (0.0, -α), 1, 0, z)
                        rhs = jacobip(n, α, β, x) / (gamma(n + α + 1) * gamma(-n - α - β))
                        @test lhs ≈ rhs atol=1e-13 rtol=1e-13
                    end
                end
            end

            @testset "BigFloat (+/- sample inputs)" begin
                αb, βb = big"0.35", big"0.6"
                for n in 0:5
                    for xb in BigFloat.(sample_real_points((-0.7, -0.2, 0.2, 0.7); seed=140 + n, lo=-0.8, hi=0.8))
                        zb = (xb - 1) / 2
                        lhsb = meijerg((big(n) + 1, -big(n) - αb - βb), (big"0.0", -αb), 1, 0, zb)
                        rhsb = jacobip(n, αb, βb, xb) / (gamma(big(n) + αb + 1) * gamma(-big(n) - αb - βb))
                        @test lhsb ≈ rhsb atol=big"1e-28" rtol=big"1e-28"
                    end
                end
            end

            @testset "Complex (+/- sample inputs)" begin
                for n in 0:5
                    for xc in sample_complex_points((-0.35 + 0.2im, 0.35 + 0.2im, -0.35 - 0.2im, 0.35 - 0.2im);
                                                   seed=150 + n,
                                                   real_lo=-0.5,
                                                   real_hi=0.5,
                                                   imag_lo=-0.3,
                                                   imag_hi=0.3)
                        zc = (xc - 1) / 2
                        lhsc = meijerg((n + 1.0, -n - α - β), (0.0, -α), 1, 0, zc)
                        rhsc = jacobip(n, α, β, xc) / (gamma(n + α + 1) * gamma(-n - α - β))
                        @test lhsc ≈ rhsc atol=1e-12 rtol=1e-12
                    end
                end
            end

            @testset "Pattern miss fallback path" begin
                a = (3.0, -2.0)
                b = (0.1, -0.35)
                for z in (-0.6 + 0im, -0.2 + 0im, 0.2, 0.6)
                    @test meijerg(a, b, 1, 0, z) ≈ meijerg_slater(a, b, 1, 0, z) atol=1e-12 rtol=1e-13
                end
            end
        end
    end

    @testset "Logarithmic confluent-pole reduction" begin
        for z in sample_real_points((-0.7, -0.2, 0.05, 0.3, 0.9); seed=160, lo=-0.85, hi=0.95)
            result = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z)
            expected = log1p(z)
            @test result ≈ expected rtol=1e-14
        end

        for z_small in sample_real_points((1e-12, 1e-10, 1e-8); seed=161, lo=1e-12, hi=1e-8)
            result_small = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_small)
            expected_small = log1p(z_small)
            @test result_small ≈ expected_small rtol=1e-14
        end

        for zb in BigFloat.(sample_real_points((-0.7, -0.2, 0.05, 0.3, 0.9); seed=162, lo=-0.85, hi=0.95))
            result_big = meijerg((BigFloat(1), BigFloat(1)), (BigFloat(1), BigFloat(0)), 1, 2, zb)
            expected_big = log1p(zb)
            @test result_big isa BigFloat
            @test result_big ≈ expected_big atol=big"1e-30" rtol=big"1e-30"
        end

        for zc in sample_complex_points((-0.5 + 0.2im, 0.2 + 0.25im, 0.7 - 0.15im);
                                        seed=163,
                                        real_lo=-0.8,
                                        real_hi=0.8,
                                        imag_lo=-0.4,
                                        imag_hi=0.4)
            result_complex = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, zc)
            expected_complex = log1p(zc)
            @test result_complex ≈ expected_complex atol=1e-12 rtol=1e-12
        end
    end
end

@testset "Hypergeometric reduction (Public API)" begin
    @testset "Lower one-term reduction to 1F1" begin
        a1 = 2.3
        b1, b2 = 0.4, 1.7

        for z in sample_real_points((0.2, 0.7, 1.4); seed=116, lo=0.1, hi=1.6)
            lhs = meijerg((a1,), (b1, b2), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1,), (1 + b1 - b2,), z) /
                  (gamma(a1 - b1) * gamma(1 + b1 - b2))
            slater = meijerg_slater((a1,), (b1, b2), 1, 0, z)
            # Public API matches the explicit 1F1-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Lower one-term reduction to 0F1" begin
        b1, b2 = 0.35, 1.6

        for z in sample_real_points((0.15, 0.6, 1.1); seed=117, lo=0.1, hi=1.3)
            lhs = meijerg((), (b1, b2), 1, 0, z)
            rhs = z^b1 * pFq((), (1 + b1 - b2,), -z) / gamma(1 + b1 - b2)
            slater = meijerg_slater((), (b1, b2), 1, 0, z)
            # Public API matches the explicit 0F1-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Lower one-term reduction to 2F2" begin
        a1, a2 = 2.2, 2.9
        b1, b2, b3 = 0.4, 1.3, 2.1

        for z in sample_real_points((0.2, 0.8, 1.5); seed=118, lo=0.1, hi=1.7)
            lhs = meijerg((a1, a2), (b1, b2, b3), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1, 1 + b1 - a2), (1 + b1 - b2, 1 + b1 - b3), -z) /
                  (gamma(a1 - b1) * gamma(a2 - b1) * gamma(1 + b1 - b2) * gamma(1 + b1 - b3))
            slater = meijerg_slater((a1, a2), (b1, b2, b3), 1, 0, z)
            # Public API matches the explicit 2F2-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Lower one-term reduction to 3F3" begin
        a1, a2, a3 = 2.1, 2.8, 3.6
        b1, b2, b3, b4 = 0.35, 1.2, 2.05, 2.9

        for z in sample_real_points((0.2, 0.7, 1.3); seed=119, lo=0.1, hi=1.5)
            lhs = meijerg((a1, a2, a3), (b1, b2, b3, b4), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1, 1 + b1 - a2, 1 + b1 - a3),
                             (1 + b1 - b2, 1 + b1 - b3, 1 + b1 - b4), z) /
                  (gamma(a1 - b1) * gamma(a2 - b1) * gamma(a3 - b1) *
                   gamma(1 + b1 - b2) * gamma(1 + b1 - b3) * gamma(1 + b1 - b4))
            slater = meijerg_slater((a1, a2, a3), (b1, b2, b3, b4), 1, 0, z)
            # Public API matches the explicit 3F3-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 1F1" begin
        a1, a2 = 0.6, 1.8
        b1 = 0.2

        for z in sample_real_points((1.3, 2.1, 3.4); seed=120, lo=1.1, hi=3.8)
            lhs = meijerg((a1, a2), (b1,), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1,), (1 - a1 + a2,), inv(z)) /
                  (gamma(a1 - b1) * gamma(1 - a1 + a2))
            slater = meijerg_slater((a1, a2), (b1,), 0, 1, z)
            # Public API matches the explicit 1F1-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 0F1" begin
        a1, a2 = 0.7, 1.9

        for z in sample_real_points((1.2, 2.0, 3.1); seed=121, lo=1.1, hi=3.5)
            lhs = meijerg((a1, a2), (), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((), (1 - a1 + a2,), -inv(z)) / gamma(1 - a1 + a2)
            slater = meijerg_slater((a1, a2), (), 0, 1, z)
            # Public API matches the explicit 0F1-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 2F3" begin
        a1, a2, a3, a4 = 0.7, 1.8, 2.6, 3.4
        b1, b2 = 0.2, 0.9

        for z in sample_real_points((1.4, 2.2, 3.6); seed=122, lo=1.2, hi=3.9)
            lhs = meijerg((a1, a2, a3, a4), (b1, b2), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1, 1 - a1 + b2), (1 - a1 + a2, 1 - a1 + a3, 1 - a1 + a4), -inv(z)) /
                  (gamma(a1 - b1) * gamma(a1 - b2) * gamma(1 - a1 + a2) * gamma(1 - a1 + a3) * gamma(1 - a1 + a4))
            slater = meijerg_slater((a1, a2, a3, a4), (b1, b2), 0, 1, z)
            # Public API matches the explicit 2F3-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 3F4" begin
        a1, a2, a3, a4, a5 = 0.8, 1.7, 2.5, 3.3, 4.2
        b1, b2, b3 = 0.2, 0.9, 1.4

        for z in sample_real_points((1.5, 2.4, 3.8); seed=123, lo=1.2, hi=4.1)
            lhs = meijerg((a1, a2, a3, a4, a5), (b1, b2, b3), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1, 1 - a1 + b2, 1 - a1 + b3),
                                   (1 - a1 + a2, 1 - a1 + a3, 1 - a1 + a4, 1 - a1 + a5), inv(z)) /
                  (gamma(a1 - b1) * gamma(a1 - b2) * gamma(a1 - b3) *
                   gamma(1 - a1 + a2) * gamma(1 - a1 + a3) * gamma(1 - a1 + a4) * gamma(1 - a1 + a5))
            slater = meijerg_slater((a1, a2, a3, a4, a5), (b1, b2, b3), 0, 1, z)
            # Public API matches the explicit 3F4-form reduction.
            @test lhs ≈ rhs rtol=1e-14
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater rtol=1e-14
        end
    end

    @testset "Complex z consistency" begin
        for z_lower in sample_complex_points((0.2 + 0.1im, 0.45 + 0.30im, 0.9 - 0.2im);
                                             seed=124,
                                             real_lo=0.1,
                                             real_hi=1.0,
                                             imag_lo=-0.4,
                                             imag_hi=0.4)
            lhs_lower = meijerg((2.3,), (0.4, 1.7), 1, 0, z_lower)
            rhs_lower = z_lower^0.4 * pFq((1 + 0.4 - 2.3,), (1 + 0.4 - 1.7,), z_lower) /
                        (gamma(2.3 - 0.4) * gamma(1 + 0.4 - 1.7))
            slater_lower = meijerg_slater((2.3,), (0.4, 1.7), 1, 0, z_lower)
            # Public API matches lower one-term pFq reduction for complex z.
            @test lhs_lower ≈ rhs_lower rtol=1e-14
            # Public API and Slater agree on the same complex lower-mode case.
            @test lhs_lower ≈ slater_lower rtol=1e-14
        end

        for z_upper in sample_complex_points((1.3 + 0.1im, 1.8 + 0.5im, 2.8 - 0.3im);
                                             seed=125,
                                             real_lo=1.2,
                                             real_hi=3.2,
                                             imag_lo=-0.6,
                                             imag_hi=0.6)
            lhs_upper = meijerg((0.6, 1.8), (0.2,), 0, 1, z_upper)
            rhs_upper = z_upper^(0.6 - 1) * pFq((1 - 0.6 + 0.2,), (1 - 0.6 + 1.8,), inv(z_upper)) /
                        (gamma(0.6 - 0.2) * gamma(1 - 0.6 + 1.8))
            slater_upper = meijerg_slater((0.6, 1.8), (0.2,), 0, 1, z_upper)
            # Public API matches upper one-term pFq reduction for complex z.
            @test lhs_upper ≈ rhs_upper rtol=1e-14
            # Public API and Slater agree on the same complex upper-mode case.
            @test lhs_upper ≈ slater_upper rtol=1e-14
        end
    end
end
