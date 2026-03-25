using HypergeometricFunctions: pFq

@testset "Slater evaluator" begin
    @testset "Fallback behavior" begin
        a = (0.3, 0.8)
        b = (0.2, 1.1)
        z = 0.7
        # Split and Slater APIs agree.
        @test meijerg(a, b, 1, 1, z) ≈ meijerg_slater(a, b, 1, 1, z) atol=1e-12 rtol=1e-12
    end

    @testset "Pure Slater evaluator" begin
        a = (0.3, 0.8)
        b = (0.2, 1.1)
        z = 0.7
        # Slater matches public evaluator.
        @test meijerg_slater(a, b, 1, 1, z) ≈ meijerg(a, b, 1, 1, z) atol=1e-12 rtol=1e-12
    end

    @testset "Domain error on integer parameter differences" begin
        # Integer difference raises domain error.
        @test_throws DomainError meijerg_slater((1.0,), (0.0,), 1, 1, 0.5)
    end
end

@testset "Slater expansion vs special-function reductions" begin
    @testset "Exponential" begin
        # G_{0,1}^{1,0}(z | - ; 0) = exp(-z)
        for z in (0.3, 1.0, 2.5)
            @test meijerg_slater((), (0.0,), 1, 0, z) ≈ exp(-z) atol=1e-12 rtol=1e-12
        end

        zc = 0.6 + 0.4im
        # Exponential identity also holds for complex input.
        @test meijerg_slater((), (0.0,), 1, 0, zc) ≈ exp(-zc) atol=1e-11 rtol=1e-11
    end

    @testset "Sine" begin
        # G_{0,2}^{1,0}(z | - ; 1/2, 0) = sin(2√z)/√π
        for z in (0.1, 0.25, 1.0)
            @test meijerg_slater((), (0.5, 0.0), 1, 0, z) ≈ sin(2 * sqrt(z)) / sqrt(pi) atol=1e-11 rtol=1e-11
        end
    end

    @testset "Cosine" begin
        # G_{0,2}^{1,0}(z | - ; 0, 1/2) = cos(2√z)/√π
        for z in (0.1, 0.25, 1.0)
            @test meijerg_slater((), (0.0, 0.5), 1, 0, z) ≈ cos(2 * sqrt(z)) / sqrt(pi) atol=1e-11 rtol=1e-11
        end
    end

    @testset "Bessel K" begin
        # G_{0,2}^{2,0}(z | - ; ν/2, -ν/2) = 2K_ν(2√z)  (non-integer ν avoids confluent poles)
        ν = 0.8
        for z in (0.5, 1.2, 3.0)
            @test meijerg_slater((), (ν / 2, -ν / 2), 2, 0, z) ≈ 2 * besselk(ν, 2 * sqrt(z)) atol=1e-11 rtol=1e-11
        end
    end

    @testset "Order reduction" begin
        # G_{1,2}^{1,1}(z | 1 ; 0.25, 1) reduces to G_{0,1}^{1,0}(z | - ; 0.25)
        # via the a[1] = b[2] = 1 cancellation
        for z in (0.3, 0.6, 0.9)
            @test meijerg_slater((1.0,), (0.25, 1.0), 1, 1, z) ≈ meijerg_slater((), (0.25,), 1, 0, z) atol=1e-12 rtol=1e-12
        end
    end
end

@testset "Hypergeometric reduction (Slater)" begin
    @testset "Lower one-term reduction to 1F1" begin
        a1 = 2.3
        b1, b2 = 0.4, 1.7

        for z in (0.2, 0.7, 1.4)
            lhs = meijerg_slater((a1,), (b1, b2), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1,), (1 + b1 - b2,), z) /
                  (gamma(a1 - b1) * gamma(1 + b1 - b2))
            # Slater lower expansion matches explicit 1F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Lower one-term reduction to 0F1" begin
        b1, b2 = 0.35, 1.6

        for z in (0.15, 0.6, 1.1)
            lhs = meijerg_slater((), (b1, b2), 1, 0, z)
            rhs = z^b1 * pFq((), (1 + b1 - b2,), -z) / gamma(1 + b1 - b2)
            # Slater lower expansion matches explicit 0F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Lower one-term reduction to 2F2" begin
        a1, a2 = 2.2, 2.9
        b1, b2, b3 = 0.4, 1.3, 2.1

        for z in (0.2, 0.8, 1.5)
            lhs = meijerg_slater((a1, a2), (b1, b2, b3), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1, 1 + b1 - a2), (1 + b1 - b2, 1 + b1 - b3), -z) /
                  (gamma(a1 - b1) * gamma(a2 - b1) * gamma(1 + b1 - b2) * gamma(1 + b1 - b3))
            # Slater lower expansion matches explicit 2F2-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Lower one-term reduction to 3F3" begin
        a1, a2, a3 = 2.1, 2.8, 3.6
        b1, b2, b3, b4 = 0.35, 1.2, 2.05, 2.9

        for z in (0.2, 0.7, 1.3)
            lhs = meijerg_slater((a1, a2, a3), (b1, b2, b3, b4), 1, 0, z)
            rhs = z^b1 * pFq((1 + b1 - a1, 1 + b1 - a2, 1 + b1 - a3),
                             (1 + b1 - b2, 1 + b1 - b3, 1 + b1 - b4), z) /
                  (gamma(a1 - b1) * gamma(a2 - b1) * gamma(a3 - b1) *
                   gamma(1 + b1 - b2) * gamma(1 + b1 - b3) * gamma(1 + b1 - b4))
            # Slater lower expansion matches explicit 3F3-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Upper one-term reduction to 1F1" begin
        a1, a2 = 0.6, 1.8
        b1 = 0.2

        for z in (1.3, 2.1, 3.4)
            lhs = meijerg_slater((a1, a2), (b1,), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1,), (1 - a1 + a2,), inv(z)) /
                  (gamma(a1 - b1) * gamma(1 - a1 + a2))
            # Slater upper expansion matches explicit 1F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Upper one-term reduction to 0F1" begin
        a1, a2 = 0.7, 1.9

        for z in (1.2, 2.0, 3.1)
            lhs = meijerg_slater((a1, a2), (), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((), (1 - a1 + a2,), -inv(z)) / gamma(1 - a1 + a2)
            # Slater upper expansion matches explicit 0F1-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Upper one-term reduction to 2F3" begin
        a1, a2, a3, a4 = 0.7, 1.8, 2.6, 3.4
        b1, b2 = 0.2, 0.9

        for z in (1.4, 2.2, 3.6)
            lhs = meijerg_slater((a1, a2, a3, a4), (b1, b2), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1, 1 - a1 + b2), (1 - a1 + a2, 1 - a1 + a3, 1 - a1 + a4), -inv(z)) /
                  (gamma(a1 - b1) * gamma(a1 - b2) * gamma(1 - a1 + a2) * gamma(1 - a1 + a3) * gamma(1 - a1 + a4))
            # Slater upper expansion matches explicit 2F3-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Upper one-term reduction to 3F4" begin
        a1, a2, a3, a4, a5 = 0.8, 1.7, 2.5, 3.3, 4.2
        b1, b2, b3 = 0.2, 0.9, 1.4

        for z in (1.5, 2.4, 3.8)
            lhs = meijerg_slater((a1, a2, a3, a4, a5), (b1, b2, b3), 0, 1, z)
            rhs = z^(a1 - 1) * pFq((1 - a1 + b1, 1 - a1 + b2, 1 - a1 + b3),
                                   (1 - a1 + a2, 1 - a1 + a3, 1 - a1 + a4, 1 - a1 + a5), inv(z)) /
                  (gamma(a1 - b1) * gamma(a1 - b2) * gamma(a1 - b3) *
                   gamma(1 - a1 + a2) * gamma(1 - a1 + a3) * gamma(1 - a1 + a4) * gamma(1 - a1 + a5))
            # Slater upper expansion matches explicit 3F4-form reduction.
            @test lhs ≈ rhs atol=1e-11 rtol=1e-11
        end
    end

    @testset "Complex z consistency" begin
        z_lower = 0.45 + 0.30im
        lhs_lower = meijerg_slater((2.3,), (0.4, 1.7), 1, 0, z_lower)
        rhs_lower = z_lower^0.4 * pFq((1 + 0.4 - 2.3,), (1 + 0.4 - 1.7,), z_lower) /
                    (gamma(2.3 - 0.4) * gamma(1 + 0.4 - 1.7))
        # Lower one-term pFq reduction remains valid for complex z.
        @test lhs_lower ≈ rhs_lower atol=1e-11 rtol=1e-11

        z_upper = 1.8 + 0.5im
        lhs_upper = meijerg_slater((0.6, 1.8), (0.2,), 0, 1, z_upper)
        rhs_upper = z_upper^(0.6 - 1) * pFq((1 - 0.6 + 0.2,), (1 - 0.6 + 1.8,), inv(z_upper)) /
                    (gamma(0.6 - 0.2) * gamma(1 - 0.6 + 1.8))
        # Upper one-term pFq reduction remains valid for complex z.
        @test lhs_upper ≈ rhs_upper atol=1e-11 rtol=1e-11
    end
end
