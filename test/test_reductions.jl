using HypergeometricFunctions: pFq

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

@testset "Reduction behavior" begin
    @testset "Order reduction equivalence" begin
        z = 0.6
        lhs = meijerg((1.0,), (0.25, 1.0), 1, 1, z)
        rhs = meijerg((), (0.25,), 1, 0, z)
        # Cancellation reduction matches direct form.
        @test lhs ≈ rhs atol=1e-12 rtol=1e-12
    end

end

@testset "Hypergeometric reduction (Public API)" begin
    @testset "Lower one-term reduction to 1F1" begin
        a1 = 2.3
        b1, b2 = 0.4, 1.7

        for z in (0.2, 0.7, 1.4)
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

        for z in (0.15, 0.6, 1.1)
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

        for z in (0.2, 0.8, 1.5)
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

        for z in (0.2, 0.7, 1.3)
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

        for z in (1.3, 2.1, 3.4)
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

        for z in (1.2, 2.0, 3.1)
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

        for z in (1.4, 2.2, 3.6)
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

        for z in (1.5, 2.4, 3.8)
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
end
