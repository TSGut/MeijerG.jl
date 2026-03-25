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
