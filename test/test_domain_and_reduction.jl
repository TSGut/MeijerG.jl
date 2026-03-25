@testset "Domain and singularity behavior" begin
    @testset "Inputs" begin
        @test_throws ArgumentError meijerg((0.2,), (0.1,), 2, 0, 0.5)
        @test_throws ArgumentError meijerg((0.2,), (0.1,), 0, 2, 0.5)
        @test_throws DomainError meijerg((0.2,), (0.1,), 1, 1, 0.0)
    end

    @testset "Confluent poles and logarithmic terms" begin
        @test_throws DomainError meijerg((1.0,), (0.0,), 1, 1, 0.5)

        z_in = 0.5
        z_out = 2.0
        # log(1+z) = G_{2,2}^{1,2}(z | 1,1 ; 1,0)
        lhs_in = meijerg((1.0, 1.0), (), (1.0,), (0.0,), z_in)
        lhs_out = meijerg((1.0, 1.0), (), (1.0,), (0.0,), z_out)
        @test lhs_in ≈ log1p(z_in) atol=1e-10 rtol=1e-10
        @test lhs_out ≈ log1p(z_out) atol=1e-10 rtol=1e-10
    end

    @testset "Order reduction equivalence" begin
        z = 0.6
        lhs = meijerg_reduce((1.0,), (0.25, 1.0), 1, 1, z)
        rhs = meijerg((), (0.25,), 1, 0, z)
        @test lhs ≈ rhs atol=1e-12 rtol=1e-12
    end

    @testset "Near |z| = 1 branch-switch" begin
        a = (0.3, 0.8)
        b = (0.2, 1.1)
        inside = meijerg(a, b, 1, 1, 0.999)
        outside = meijerg(a, b, 1, 1, 1.001)
        @test isfinite(inside)
        @test isfinite(outside)
    end
end
