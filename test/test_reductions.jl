@testset "Public Meijer G API" begin
    @testset "Order reduction mapping" begin
        z = 0.6
        reduced = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
        reference = meijerg((), (0.25,), 1, 0, z)
        # Reduced and direct forms agree.
        @test reduced ≈ reference atol=1e-12 rtol=1e-12
    end

    @testset "Elementary explicit mappings" begin
        x = 0.3
        # Exponential mapping evaluates correctly.
        @test meijerg((), (), (0.0,), (), -x) ≈ exp(x) atol=1e-12 rtol=1e-12

        y = 0.7
        # Sine mapping evaluates correctly.
        @test sqrt(pi) * meijerg((), (), (0.5,), (0.0,), y^2 / 4) ≈ sin(y) atol=1e-11 rtol=1e-11
        # Cosine mapping evaluates correctly.
        @test sqrt(pi) * meijerg((), (), (0.0,), (0.5,), y^2 / 4) ≈ cos(y) atol=1e-11 rtol=1e-11

        ν = 0.8
        z = 1.2
        lhs = meijerg((), (), (ν / 2, -ν / 2), (), z)
        rhs = 2 * besselk(ν, 2 * sqrt(z))
        # Bessel-K mapping evaluates correctly.
        @test lhs ≈ rhs atol=1e-11 rtol=1e-11
    end

    @testset "Type stability and BigFloat support" begin
        v_float = meijerg((), (), (0.0,), (), -0.2)
        # Float input returns Float output.
        @test v_float isa Float64

        x_big = BigFloat("0.2")
        v_big = meijerg((), (), (BigFloat(0),), (), -x_big)
        # BigFloat input returns BigFloat output.
        @test v_big isa BigFloat
        # BigFloat value matches exponential reference.
        @test v_big ≈ exp(x_big) atol=big"1e-30" rtol=big"1e-30"
    end

    @testset "Logarithmic special cases" begin
        # log(1+z) via G_{2,2}^{1,2}(z | 1, 1 ; 1, 0) = log(1+z)/z
        # Reference: Gradshteyn & Ryzhik 9.353, DLMF 15.8.2
        z_val = 0.3
        result = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_val)
        expected = log1p(z_val) / z_val
        # Log mapping matches reference value.
        @test result ≈ expected atol=1e-12 rtol=1e-12

        z_small = 1e-10
        result_small = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_small)
        expected_small = log1p(z_small) / z_small
        # Small-z limit is numerically stable.
        @test result_small ≈ expected_small atol=1e-10 rtol=1e-10

        z_big = BigFloat("0.3")
        result_big = meijerg((BigFloat(1), BigFloat(1)), (BigFloat(1), BigFloat(0)), 1, 2, z_big)
        # BigFloat log mapping preserves type.
        @test result_big isa BigFloat
        expected_big = log1p(z_big) / z_big
        # BigFloat log mapping matches reference.
        @test result_big ≈ expected_big atol=big"1e-30" rtol=big"1e-30"
    end
end

@testset "Elementary identities" begin
    @testset "Exponential and trigonometric reductions" begin
        x = 0.3
        # Exponential identity holds.
        @test meijerg((), (), (0,), (), -x) ≈ exp(x) atol=1e-12 rtol=1e-12

        y = 0.7
        # Sine identity holds.
        @test sqrt(pi) * meijerg((), (), (0.5,), (0,), y^2 / 4) ≈ sin(y) atol=1e-11 rtol=1e-11
        # Cosine identity holds.
        @test sqrt(pi) * meijerg((), (), (0.0,), (0.5,), y^2 / 4) ≈ cos(y) atol=1e-11 rtol=1e-11
    end

    @testset "Bessel-K reduction" begin
        ν = 0.8
        z = 1.2
        lhs = meijerg((), (), (ν / 2, -ν / 2), (), z)
        rhs = 2 * besselk(ν, 2 * sqrt(z))
        # Bessel-K identity holds.
        @test lhs ≈ rhs atol=1e-11 rtol=1e-11
    end
end

@testset "Reduction behavior moved from domain" begin
    @testset "Order reduction equivalence" begin
        z = 0.6
        lhs = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
        rhs = meijerg((), (0.25,), 1, 0, z)
        # Cancellation reduction matches direct form.
        @test lhs ≈ rhs atol=1e-12 rtol=1e-12
    end

    @testset "log(1+z) via confluent expansion" begin
        z_in = 0.5
        z_out = 2.0
        lhs_in = meijerg((1.0, 1.0), (), (1.0,), (0.0,), z_in)
        lhs_out = meijerg((1.0, 1.0), (), (1.0,), (0.0,), z_out)
        # In-domain confluent case matches log reference.
        @test lhs_in ≈ log1p(z_in) / z_in atol=1e-10 rtol=1e-10
        # Out-of-domain confluent case matches log reference.
        @test lhs_out ≈ log1p(z_out) / z_out atol=1e-10 rtol=1e-10
    end
end
