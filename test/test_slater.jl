using HypergeometricFunctions: pFq
include("test_utilities.jl")

@testset "Slater evaluator" begin
    @testset "Fallback behavior" begin
        a = (0.3, 0.8)
        b = (0.2, 1.1)
        for z in sample_real_points((0.2, 0.7, 1.1); seed=201, lo=0.1, hi=1.3)
            # Split and Slater APIs agree.
            @test meijerg(a, b, 1, 1, z) ≈ meijerg_slater(a, b, 1, 1, z) atol=1e-12 rtol=1e-14
        end
    end

    @testset "Pure Slater evaluator" begin
        a = (0.3, 0.8)
        b = (0.2, 1.1)
        for z in sample_real_points((0.2, 0.7, 1.1); seed=202, lo=0.1, hi=1.3)
            # Slater matches public evaluator.
            @test meijerg_slater(a, b, 1, 1, z) ≈ meijerg(a, b, 1, 1, z) atol=1e-12 rtol=1e-14
        end
    end

    @testset "Domain error on integer parameter differences" begin
        # Integer difference raises domain error.
        @test_throws DomainError meijerg_slater((1.0,), (0.0,), 1, 1, 0.5)
    end
end

@testset "Slater sanity checks" begin
    @testset "Sinh / Cosh (analytic continuation)" begin
        for x in sample_real_points((0.1, 0.5, 1.2); seed=300, lo=0.05, hi=2.0)
            # use a complex z for negative real bases to avoid real-power domain errors
            z = -x^2 / 4 + 0im
            g_sin = meijerg_slater((), (0.5, 0.0), 1, 0, z)
            lhs_sinh = -im * sqrt(pi) * g_sin
            @test real(lhs_sinh) ≈ sinh(x) atol=1e-12 rtol=1e-12

            g_cos = meijerg_slater((), (0.0, 0.5), 1, 0, z)
            lhs_cosh = sqrt(pi) * g_cos
            @test real(lhs_cosh) ≈ cosh(x) atol=1e-12 rtol=1e-12
        end
    end

    @testset "erf / erfc / erfcx" begin
        for x in sample_real_points((0.1, 0.6, 1.5); seed=301, lo=0.05, hi=2.0)
            z = x^2
            erf_slater = (1 / sqrt(pi)) * meijerg_slater((1.0,), (0.5, 0.0), 1, 1, z)
            @test erf_slater ≈ erf(x) atol=1e-12 rtol=1e-12

            # meijerg_slater rejects the direct upper-gamma pairing used by the
            # naive G->erfc mapping (a[j]-b[k] a positive integer).  Use the
            # relation erfc = 1 - erf for the slater-based sanity check instead.
            erfc_slater = 1 - erf_slater
            @test erfc_slater ≈ erfc(x) atol=1e-12 rtol=1e-12

            erfcx_slater = exp(z) * erfc_slater
            @test erfcx_slater ≈ exp(x^2) * erfc(x) atol=1e-12 rtol=1e-12
        end
    end

    @testset "Modified Bessel I (nu=0)" begin
        ν = 0.0
        for x in sample_real_points((0.05, 0.7, 2.0); seed=302, lo=0.02, hi=3.0)
            # use complex z for negative-real bases (analytic continuation)
            z = -x^2 / 4 + 0im
            g_val = meijerg_slater((), (ν / 2, -ν / 2), 1, 0, z)
            # For ν=0 the analytic-continuation prefactor i^{-ν}=1 and J_0(i x)=I_0(x).
            @test g_val ≈ besseli(ν, x) atol=1e-11 rtol=1e-12
        end
    end

    @testset "Fresnel via erf mapping" begin
        for x in sample_real_points((0.1, 0.6, 1.1); seed=303, lo=0.05, hi=2.0)
            zarg = exp(1im * pi / 4) * sqrt(pi / 2) * x
            # compute erf via slater mapping then form C+iS
            erf_slater = (1 / sqrt(pi)) * meijerg_slater((1.0,), (0.5, 0.0), 1, 1, zarg^2)
            cs_slater = (1 + 1im) / 2 * erf_slater
            cs_expected = (1 + 1im) / 2 * erf(zarg)
            @test cs_slater ≈ cs_expected atol=1e-11 rtol=1e-12
        end
    end

    @testset "Airy Ai via Bessel-K mapping" begin
        for x in sample_real_points((0.1, 0.9, 2.4); seed=304, lo=0.02, hi=3.0)
            z = x^3 / 9
            g_val = meijerg_slater((), (1/6, -1/6), 2, 0, z)
            ai_slater = (1 / (2 * pi)) * sqrt(x / 3) * g_val
            ai_expected = (1 / pi) * sqrt(x / 3) * besselk(1/3, (2 / 3) * x^(3/2))
            @test ai_slater ≈ ai_expected atol=1e-11 rtol=1e-12
        end
    end
end

@testset "Slater expansion vs special-function reductions" begin
    @testset "Exponential" begin
        # G_{0,1}^{1,0}(z | - ; 0) = exp(-z)
        for z in sample_real_points((0.3, 1.0, 2.5); seed=203, lo=0.1, hi=2.8)
            @test meijerg_slater((), (0.0,), 1, 0, z) ≈ exp(-z) atol=1e-12 rtol=1e-14
        end

        for zc in sample_complex_points((0.2 + 0.1im, 0.6 + 0.4im, 1.0 - 0.2im);
                                        seed=204,
                                        real_lo=0.1,
                                        real_hi=1.2,
                                        imag_lo=-0.4,
                                        imag_hi=0.4)
            # Exponential identity also holds for complex input.
            @test meijerg_slater((), (0.0,), 1, 0, zc) ≈ exp(-zc) atol=1e-11 rtol=1e-14
        end
    end

    @testset "Sine" begin
        # G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2√z)/√π
        for z in sample_real_points((0.1, 0.25, 1.0); seed=205, lo=0.05, hi=1.2)
            @test meijerg_slater((), (0.5, 0.0), 1, 0, z) ≈ sin(2 * sqrt(z)) / sqrt(pi) atol=1e-11 rtol=1e-14
        end
    end

    @testset "Cosine" begin
        # G_{0,2}^{1,0}(z | - ; 0, 1/2) = cos(2√z)/√π
        for z in sample_real_points((0.1, 0.25, 1.0); seed=206, lo=0.05, hi=1.2)
            @test meijerg_slater((), (0.0, 0.5), 1, 0, z) ≈ cos(2 * sqrt(z)) / sqrt(pi) atol=1e-11 rtol=1e-14
        end
    end

    @testset "Bessel K" begin
        # G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2√z)  (non-integer ν avoids confluent poles)
        ν = 0.8
        for z in sample_real_points((0.5, 1.2, 3.0); seed=207, lo=0.2, hi=3.3)
            @test meijerg_slater((), (ν / 2, -ν / 2), 2, 0, z) ≈ 2 * besselk(ν, 2 * sqrt(z)) atol=1e-11 rtol=1e-14
        end
    end

    @testset "Order reduction" begin
        # G_{1,2}^{1,1}(z | 1 ; 0.25, 1) reduces to G_{0,1}^{1,0}(z | - ; 0.25)
        # via the a[1] = b[2] = 1 cancellation
        for z in sample_real_points((0.3, 0.6, 0.9); seed=208, lo=0.1, hi=1.1)
            @test meijerg_slater((1.0,), (0.25, 1.0), 1, 1, z) ≈ meijerg_slater((), (0.25,), 1, 0, z) atol=1e-12 rtol=1e-14
        end
    end
end

@testset "Hypergeometric reduction (Slater)" begin
    @testset "Lower one-term reduction to 1F1" begin
        a1 = 2.3
        b1, b2 = 0.4, 1.7

        for z in sample_real_points((0.2, 0.7, 1.4); seed=209, lo=0.1, hi=1.6)
            lhs = meijerg_slater((a1,), (b1, b2), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1,), (1 + b1 - b2,), z) /
                  (gamma(a1 - b1) * gamma(1 + b1 - b2))
            # Slater lower expansion matches explicit 1F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Lower one-term reduction to 0F1" begin
        b1, b2 = 0.35, 1.6

        for z in sample_real_points((0.15, 0.6, 1.1); seed=210, lo=0.1, hi=1.3)
            lhs = meijerg_slater((), (b1, b2), 1, 0, z)
            rhs = z^b1 * pFq((), (1 + b1 - b2,), -z) / gamma(1 + b1 - b2)
            # Slater lower expansion matches explicit 0F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Lower one-term reduction to 2F2" begin
        a1, a2 = 2.2, 2.9
        b1, b2, b3 = 0.4, 1.3, 2.1

        for z in sample_real_points((0.2, 0.8, 1.5); seed=211, lo=0.1, hi=1.7)
            lhs = meijerg_slater((a1, a2), (b1, b2, b3), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1, 1 + b1 - a2), (1 + b1 - b2, 1 + b1 - b3), -z) /
                  (gamma(a1 - b1) * gamma(a2 - b1) * gamma(1 + b1 - b2) * gamma(1 + b1 - b3))
            # Slater lower expansion matches explicit 2F2-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Lower one-term reduction to 3F3" begin
        a1, a2, a3 = 2.1, 2.8, 3.6
        b1, b2, b3, b4 = 0.35, 1.2, 2.05, 2.9

        for z in sample_real_points((0.2, 0.7, 1.3); seed=212, lo=0.1, hi=1.5)
            lhs = meijerg_slater((a1, a2, a3), (b1, b2, b3, b4), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1, 1 + b1 - a2, 1 + b1 - a3),
                             (1 + b1 - b2, 1 + b1 - b3, 1 + b1 - b4), z) /
                  (gamma(a1 - b1) * gamma(a2 - b1) * gamma(a3 - b1) *
                   gamma(1 + b1 - b2) * gamma(1 + b1 - b3) * gamma(1 + b1 - b4))
            # Slater lower expansion matches explicit 3F3-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 1F1" begin
        a1, a2 = 0.6, 1.8
        b1 = 0.2

        for z in sample_real_points((1.3, 2.1, 3.4); seed=213, lo=1.1, hi=3.8)
            lhs = meijerg_slater((a1, a2), (b1,), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1,), (1 - a1 + a2,), inv(z)) /
                  (gamma(a1 - b1) * gamma(1 - a1 + a2))
            # Slater upper expansion matches explicit 1F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 0F1" begin
        a1, a2 = 0.7, 1.9

        for z in sample_real_points((1.2, 2.0, 3.1); seed=214, lo=1.1, hi=3.5)
            lhs = meijerg_slater((a1, a2), (), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((), (1 - a1 + a2,), -inv(z)) / gamma(1 - a1 + a2)
            # Slater upper expansion matches explicit 0F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 2F3" begin
        a1, a2, a3, a4 = 0.7, 1.8, 2.6, 3.4
        b1, b2 = 0.2, 0.9

        for z in sample_real_points((1.4, 2.2, 3.6); seed=215, lo=1.2, hi=3.9)
            lhs = meijerg_slater((a1, a2, a3, a4), (b1, b2), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1, 1 - a1 + b2), (1 - a1 + a2, 1 - a1 + a3, 1 - a1 + a4), -inv(z)) /
                  (gamma(a1 - b1) * gamma(a1 - b2) * gamma(1 - a1 + a2) * gamma(1 - a1 + a3) * gamma(1 - a1 + a4))
            # Slater upper expansion matches explicit 2F3-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Upper one-term reduction to 3F4" begin
        a1, a2, a3, a4, a5 = 0.8, 1.7, 2.5, 3.3, 4.2
        b1, b2, b3 = 0.2, 0.9, 1.4

        for z in sample_real_points((1.5, 2.4, 3.8); seed=216, lo=1.2, hi=4.1)
            lhs = meijerg_slater((a1, a2, a3, a4, a5), (b1, b2, b3), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1, 1 - a1 + b2, 1 - a1 + b3),
                                   (1 - a1 + a2, 1 - a1 + a3, 1 - a1 + a4, 1 - a1 + a5), inv(z)) /
                  (gamma(a1 - b1) * gamma(a1 - b2) * gamma(a1 - b3) *
                   gamma(1 - a1 + a2) * gamma(1 - a1 + a3) * gamma(1 - a1 + a4) * gamma(1 - a1 + a5))
            # Slater upper expansion matches explicit 3F4-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-14
        end
    end

    @testset "Complex z consistency" begin
        for z_lower in sample_complex_points((0.2 + 0.1im, 0.45 + 0.30im, 0.9 - 0.2im);
                                             seed=217,
                                             real_lo=0.1,
                                             real_hi=1.0,
                                             imag_lo=-0.4,
                                             imag_hi=0.4)
            lhs_lower = meijerg_slater((2.3,), (0.4, 1.7), 1, 0, z_lower)
            rhs_lower = z_lower^0.4 * pFq((1 + 0.4 - 2.3,), (1 + 0.4 - 1.7,), z_lower) /
                        (gamma(2.3 - 0.4) * gamma(1 + 0.4 - 1.7))
            # Lower one-term pFq reduction remains valid for complex z.
            @test lhs_lower ≈ rhs_lower atol=1e-11 rtol=1e-14
        end

        for z_upper in sample_complex_points((1.3 + 0.1im, 1.8 + 0.5im, 2.8 - 0.3im);
                                             seed=218,
                                             real_lo=1.2,
                                             real_hi=3.2,
                                             imag_lo=-0.6,
                                             imag_hi=0.6)
            lhs_upper = meijerg_slater((0.6, 1.8), (0.2,), 0, 1, z_upper)
            rhs_upper = z_upper^(0.6 - 1) * pFq((1 - 0.6 + 0.2,), (1 - 0.6 + 1.8,), inv(z_upper)) /
                        (gamma(0.6 - 0.2) * gamma(1 - 0.6 + 1.8))
            # Upper one-term pFq reduction remains valid for complex z.
            @test lhs_upper ≈ rhs_upper atol=1e-11 rtol=1e-14
        end
    end
    @testset "Slater NaN-avoidance" begin
        @testset "Skip pFq when gamma-ratio zero" begin
            val = meijerg_slater((1/3, 2/3), (1.0, 2.0), 1, 1, 0.5 + 0.5im)
            @test isfinite(real(val)) && isfinite(imag(val))
            @test abs(val) ≈ 0.0 atol=1e-14
        end
    end
end


