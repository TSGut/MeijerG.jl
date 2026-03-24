@testset "Elementary identities" begin
    @testset "Exponential and trigonometric reductions" begin
        x = 0.3
        @test meijerg((), (), (0,), (), -x) ≈ exp(x) atol=1e-12 rtol=1e-12

        y = 0.7
        @test sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2 / 4) ≈ sin(y) atol=1e-11 rtol=1e-11
        @test sqrt(pi) * meijerg((), (), (0.0,), (0.5,), y^2 / 4) ≈ cos(y) atol=1e-11 rtol=1e-11
    end

    @testset "Bessel-K reduction" begin
        ν = 0.8
        z = 1.2
        lhs = meijerg((), (), (ν / 2, -ν / 2), (), z)
        rhs = 2 * besselk(ν, 2 * sqrt(z))
        @test lhs ≈ rhs atol=1e-11 rtol=1e-11
    end
end
