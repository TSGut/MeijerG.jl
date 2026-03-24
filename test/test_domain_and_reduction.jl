@testset "Domain, singularity, and reduction behavior" begin
    @testset "Inputs" begin
        @test_throws ArgumentError meijerg((0.2,), (0.1,), 2, 0, 0.5)
        @test_throws ArgumentError meijerg((0.2,), (0.1,), 0, 2, 0.5)
        @test_throws DomainError meijerg((0.2,), (0.1,), 1, 1, 0.0)
    end

    @testset "Confluent/simple-pole guards" begin
        @test_throws DomainError meijerg((0.2,), (0.1, 0.1), 2, 1, 0.5)
        @test_throws DomainError meijerg((0.2, 0.2), (0.1,), 1, 2, 2.0)
    end

    @testset "Order reduction equivalence" begin
        z = 0.6
        lhs = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
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
