@testset "Explicit reduction API" begin
    @testset "Fallback behavior" begin
        a = (0.3, 0.8)
        b = (0.2, 1.1)
        z = 0.7
        @test meijerg_reduce(a, b, 1, 1, z) ≈ meijerg(a, b, 1, 1, z) atol=1e-12 rtol=1e-12
    end

    @testset "Order reduction mapping" begin
        z = 0.6
        reduced = meijerg_reduce((1.0,), (0.25, 1.0), 1, 1, z)
        reference = meijerg((), (0.25,), 1, 0, z)
        @test reduced ≈ reference atol=1e-12 rtol=1e-12
    end

    @testset "Elementary explicit mappings" begin
        x = 0.3
        @test meijerg_reduce((), (), (0.0,), (), -x) ≈ exp(x) atol=1e-12 rtol=1e-12

        y = 0.7
        @test sqrt(pi) * meijerg_reduce((), (), (0.5,), (0.0,), y^2 / 4) ≈ sin(y) atol=1e-11 rtol=1e-11
        @test sqrt(pi) * meijerg_reduce((), (), (0.0,), (0.5,), y^2 / 4) ≈ cos(y) atol=1e-11 rtol=1e-11

        ν = 0.8
        z = 1.2
        lhs = meijerg_reduce((), (), (ν / 2, -ν / 2), (), z)
        rhs = 2 * besselk(ν, 2 * sqrt(z))
        @test lhs ≈ rhs atol=1e-11 rtol=1e-11
    end

    @testset "Type stability and BigFloat support" begin
        v_float = meijerg_reduce((), (), (0.0,), (), -0.2)
        @test v_float isa Float64

        x_big = BigFloat("0.2")
        v_big = meijerg_reduce((), (), (BigFloat(0),), (), -x_big)
        @test v_big isa BigFloat
        @test v_big ≈ exp(x_big) atol=big"1e-30" rtol=big"1e-30"
    end
end
