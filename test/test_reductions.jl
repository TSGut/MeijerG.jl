using HypergeometricFunctions: pFq
using ClassicalOrthogonalPolynomials: laguerrel, hermiteh
using SpecialFunctions: gamma, besselj, besselk
include("test_utilities.jl")

@testset "Public Meijer G API" begin
    @testset "Order reduction mapping" begin
        for z in sample_real_points((0.15, 0.6, 1.2); seed=101, lo=0.1, hi=1.4)
            reduced = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
            reference = meijerg((), (0.25,), 1, 0, z)
            # Reduced and direct forms agree.
            @test reduced ≈ reference atol=1e-12 rtol=1e-12
        end
    end

    @testset "Elementary explicit mappings" begin
        for x in sample_real_points((0.05, 0.3, 0.9); seed=102, lo=0.02, hi=1.0)
            # Exponential mapping evaluates correctly.
            @test meijerg((), (), (0.0,), (), -x) ≈ exp(x) atol=1e-12 rtol=1e-12
        end

        for y in sample_real_points((0.1, 0.7, 1.4); seed=103, lo=0.05, hi=1.6)
            # Sine mapping evaluates correctly.
            @test sqrt(pi) * meijerg((), (), (0.5,), (0.0,), y^2 / 4) ≈ sin(y) atol=1e-11 rtol=1e-11
            # Cosine mapping evaluates correctly.
            @test sqrt(pi) * meijerg((), (), (0.0,), (0.5,), y^2 / 4) ≈ cos(y) atol=1e-11 rtol=1e-11
        end

        ν = 0.8
        for z in sample_real_points((0.15, 1.2, 2.4); seed=104, lo=0.1, hi=2.8)
            lhs = meijerg((), (), (ν / 2, -ν / 2), (), z)
            rhs = 2 * besselk(ν, 2 * sqrt(z))
            # Bessel-K mapping evaluates correctly.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11

            lhs_j = meijerg((), (), (ν / 2,), (-ν / 2,), z)
            rhs_j = besselj(ν, 2 * sqrt(z))
            # Bessel-J mapping evaluates correctly.
            @test lhs_j ≈ rhs_j atol=1e-11 rtol=1e-11
        end
    end

    @testset "Type stability and BigFloat support" begin
        for x in sample_real_points((0.02, 0.2, 0.8); seed=105, lo=0.01, hi=0.9)
            v_float = meijerg((), (), (0.0,), (), -x)
            # Float input returns Float output.
            @test v_float isa Float64
            # Float value matches exponential reference.
            @test v_float ≈ exp(x) atol=1e-12 rtol=1e-12
        end

        for x_big in BigFloat.(sample_real_points((0.02, 0.2, 0.8); seed=106, lo=0.01, hi=0.9))
            v_big = meijerg((), (), (BigFloat(0),), (), -x_big)
            # BigFloat input returns BigFloat output.
            @test v_big isa BigFloat
            # BigFloat value matches exponential reference.
            @test v_big ≈ exp(x_big) atol=big"1e-30" rtol=big"1e-30"
        end

        for z_complex in sample_complex_points((0.2 + 0.3im, 0.7 + 0.1im);
                                               seed=107,
                                               real_lo=0.1,
                                               real_hi=1.0,
                                               imag_lo=-0.4,
                                               imag_hi=0.4)
            v_complex = meijerg((), (0.0,), 1, 0, z_complex)
            # Complex input returns complex output.
            @test v_complex isa ComplexF64
            # Complex exponential identity is preserved.
            @test v_complex ≈ exp(-z_complex) atol=1e-11 rtol=1e-11
        end
    end

    @testset "Logarithmic special cases" begin
        # log(1+z) via G_{2,2}^{1,2}(z | 1, 1 ; 1, 0) = log(1+z)/z
        # Reference: Gradshteyn & Ryzhik 9.353, DLMF 15.8.2
        for z_val in sample_real_points((0.05, 0.3, 0.9); seed=108, lo=0.02, hi=0.95)
            result = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_val)
            expected = log1p(z_val) / z_val
            # Log mapping matches reference value.
            @test result ≈ expected atol=1e-12 rtol=1e-12
        end

        for z_small in sample_real_points((1e-12, 1e-10, 1e-8); seed=109, lo=1e-12, hi=1e-8)
            result_small = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_small)
            expected_small = log1p(z_small) / z_small
            # Small-z limit is numerically stable.
            @test result_small ≈ expected_small atol=1e-10 rtol=1e-10
        end

        for z_big in BigFloat.(sample_real_points((0.05, 0.3, 0.9); seed=110, lo=0.02, hi=0.95))
            result_big = meijerg((BigFloat(1), BigFloat(1)), (BigFloat(1), BigFloat(0)), 1, 2, z_big)
            # BigFloat log mapping preserves type.
            @test result_big isa BigFloat
            expected_big = log1p(z_big) / z_big
            # BigFloat log mapping matches reference.
            @test result_big ≈ expected_big atol=big"1e-30" rtol=big"1e-30"
        end
    end

    @testset "Incomplete gamma special cases" begin
        a = 0.7

        for z_val in sample_real_points((0.2, 1.3, 2.8); seed=111, lo=0.1, hi=3.0)
            lower_result = meijerg((1.0,), (a, 0.0), 1, 1, z_val)
            lower_expected = gamma(a) - gamma(a, z_val)
            # Lower incomplete gamma mapping matches reference value.
            @test lower_result ≈ lower_expected atol=1e-12 rtol=1e-12

            upper_result = meijerg((1.0,), (a, 0.0), 2, 1, z_val)
            upper_expected = gamma(a, z_val)
            # Upper incomplete gamma mapping matches reference value.
            @test upper_result ≈ upper_expected atol=1e-12 rtol=1e-12
        end
    end
end

@testset "Elementary identities" begin
    @testset "Exponential and trigonometric reductions" begin
        for x in sample_real_points((0.05, 0.3, 0.9); seed=112, lo=0.02, hi=1.0)
            # Exponential identity holds.
            @test meijerg((), (), (0,), (), -x) ≈ exp(x) atol=1e-12 rtol=1e-12
        end

        for y in sample_real_points((0.1, 0.7, 1.4); seed=113, lo=0.05, hi=1.6)
            # Sine identity holds.
            @test sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2 / 4) ≈ sin(y) atol=1e-11 rtol=1e-11
            # Cosine identity holds.
            @test sqrt(pi) * meijerg((), (), (0.0,), (0.5,), y^2 / 4) ≈ cos(y) atol=1e-11 rtol=1e-11
        end
    end

    @testset "Bessel-K reduction" begin
        ν = 0.8
        for z in sample_real_points((0.15, 1.2, 2.4); seed=114, lo=0.1, hi=2.8)
            lhs = meijerg((), (), (ν / 2, -ν / 2), (), z)
            rhs = 2 * besselk(ν, 2 * sqrt(z))
            # Bessel-K identity holds.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end
end

@testset "Reduction behavior" begin
    @testset "Order reduction equivalence" begin
        for z in sample_real_points((0.15, 0.6, 1.2); seed=115, lo=0.1, hi=1.4)
            lhs = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
            rhs = meijerg((), (0.25,), 1, 0, z)
            # Cancellation reduction matches direct form.
            @test lhs ≈ rhs atol=1e-12 rtol=1e-12
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
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
        end
    end

    @testset "Lower one-term reduction to 0F1" begin
        b1, b2 = 0.35, 1.6

        for z in sample_real_points((0.15, 0.6, 1.1); seed=117, lo=0.1, hi=1.3)
            lhs = meijerg((), (b1, b2), 1, 0, z)
            rhs = z^b1 * pFq((), (1 + b1 - b2,), -z) / gamma(1 + b1 - b2)
            slater = meijerg_slater((), (b1, b2), 1, 0, z)
            # Public API matches the explicit 0F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
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
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
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
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
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
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
        end
    end

    @testset "Upper one-term reduction to 0F1" begin
        a1, a2 = 0.7, 1.9

        for z in sample_real_points((1.2, 2.0, 3.1); seed=121, lo=1.1, hi=3.5)
            lhs = meijerg((a1, a2), (), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((), (1 - a1 + a2,), -inv(z)) / gamma(1 - a1 + a2)
            slater = meijerg_slater((a1, a2), (), 0, 1, z)
            # Public API matches the explicit 0F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
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
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
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
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same pFq case.
            @test lhs ≈ slater atol=1e-12 rtol=1e-12
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
            @test lhs_lower ≈ rhs_lower atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same complex lower-mode case.
            @test lhs_lower ≈ slater_lower atol=1e-12 rtol=1e-12
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
            @test lhs_upper ≈ rhs_upper atol=1e-11 rtol=1e-11
            # Public API and Slater agree on the same complex upper-mode case.
            @test lhs_upper ≈ slater_upper atol=1e-12 rtol=1e-12
        end
    end
end

@testset "Polynomial and orthogonal-polynomial reductions" begin
    @testset "Generalized Laguerre polynomial reductions" begin
        # Via 1F1 reduction and L_n^(α)(x) = ((α+1)_n / n!) 1F1(-n; α+1; x):
        # G_{1,2}^{1,0}(x | n+1 ; 0, -α) = L_n^(α)(x) / Γ(n+α+1)
        x_values = sample_real_points((0.05, 0.4, 1.2, 2.3); seed=126, lo=0.02, hi=2.6)
        for α in (0.2, 0.7, 1.3)
            for n in 0:6
                for x in x_values
                    lhs = gamma(n + α + 1) * meijerg((n + 1.0,), (0.0, -α), 1, 0, x)
                    rhs = laguerrel(n, α, x)
                    @test lhs ≈ rhs atol=1e-11 rtol=1e-11
                end
            end
        end
    end

    @testset "Hermite polynomial reductions" begin
        # DLMF 13.6(v) / confluent-hypergeometric identities:
        # H_{2n}(x)   = (-1)^n (2n)!/n! * 1F1(-n; 1/2; x^2)
        # H_{2n+1}(x) = (-1)^n (2n+1)!/n! * 2x * 1F1(-n; 3/2; x^2)
        # Combined with lower one-term Slater reduction for 1F1.
        x_values = sample_real_points((0.05, 0.2, 0.8, 1.4); seed=127, lo=0.02, hi=1.6)
        for n in 0:5
            for x in x_values
                g_even = meijerg((n + 1.0,), (0.0, 0.5), 1, 0, x^2)
                h_even = (-1)^n * factorial(big(2n)) * sqrt(big(pi)) * g_even
                @test h_even ≈ hermiteh(2n, x) atol=1e-10 rtol=1e-10

                g_odd = meijerg((n + 1.0,), (0.0, -0.5), 1, 0, x^2)
                h_odd = (-1)^n * factorial(big(2n + 1)) * x * sqrt(big(pi)) * g_odd
                @test h_odd ≈ hermiteh(2n + 1, x) atol=1e-10 rtol=1e-10
            end
        end
    end
end