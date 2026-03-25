@testset "Domain behavior" begin
    @testset "Inputs" begin
        # Invalid m index raises argument error.
        @test_throws ArgumentError meijerg((0.2,), (0.1,), 2, 0, 0.5)
        # Invalid n index raises argument error.
        @test_throws ArgumentError meijerg((0.2,), (0.1,), 0, 2, 0.5)
        # Zero input is outside supported domain.
        @test_throws DomainError meijerg((0.2,), (0.1,), 1, 1, 0.0)
    end

    @testset "Near |z| = 1 branch-switch" begin
        a = (0.3, 0.8)
        b = (0.2, 1.1)
        inside = meijerg(a, b, 1, 1, 0.999)
        outside = meijerg(a, b, 1, 1, 1.001)
        # Value remains finite just inside boundary.
        @test isfinite(inside)
        # Value remains finite just outside boundary.
        @test isfinite(outside)
    end
end
