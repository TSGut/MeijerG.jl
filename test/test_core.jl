@testset "Core polyalgorithm" begin
    relerr(x, y) = abs(x - y) / max(abs(y), eps(typeof(abs(y))))

    @testset "Slater path: low confluence risk, all input types" begin
        a, b, m, n, z = (0.5,), (0.2, 0.9), 2, 1, 0.5
        # High-precision reference via BigFloat.
        ref = setprecision(512) do
            meijerg(map(BigFloat, a), map(BigFloat, b), m, n, BigFloat(z))
        end

        # Float64: result is real, matches reference closely.
        result_f64 = meijerg(a, b, m, n, z)
        @test result_f64 isa Float64
        @test isfinite(result_f64)
        @test relerr(BigFloat(result_f64), ref) < big"1e-12"

        # Complex{Float64}: real input path returns complex type, real part matches.
        result_c64 = meijerg(a, b, m, n, z + 0im)
        @test result_c64 isa Complex{Float64}
        @test abs(result_c64 - result_f64) < 1e-12

        # BigFloat: result is BigFloat and matches high-precision reference.
        result_bf = meijerg(map(BigFloat, a), map(BigFloat, b), m, n, BigFloat(z))
        @test result_bf isa BigFloat
        @test isfinite(result_bf)
        @test relerr(result_bf, ref) < big"1e-30"

        # Complex{BigFloat}: result is complex, real part matches BigFloat result.
        result_cbf = meijerg(map(BigFloat, a), map(BigFloat, b), m, n,
                             BigFloat(z) + zero(BigFloat) * im)
        @test result_cbf isa Complex{BigFloat}
        @test abs(result_cbf - result_bf) < big"1e-30"
    end

    @testset "Loop contour lower: integer b-parameter difference, all input types" begin
        a, b, m, n = (0.5,), (0.01, 1.01), 2, 1

        for z in (0.3, 0.5, 0.8)
            # Float64 input.
            result_f64 = meijerg(a, b, m, n, z)
            @test result_f64 isa Number
            @test isfinite(result_f64)

            # Complex{Float64}: passes a complex z, returns complex type.
            result_c64 = meijerg(a, b, m, n, z + 0im)
            @test result_c64 isa Complex{Float64}
            @test abs(result_c64 - result_f64) < 1e-8 * max(abs(result_f64), 1e-10)

            # BigFloat: higher precision value (possibly complex due to contour path).
            result_bf = meijerg(map(BigFloat, a), map(BigFloat, b), m, n, BigFloat(z))
            @test result_bf isa Number
            @test isfinite(result_bf)

            # Complex{BigFloat}: should match BigFloat-path value.
            result_cbf = meijerg(map(BigFloat, a), map(BigFloat, b), m, n,
                                 BigFloat(z) + zero(BigFloat) * im)
            @test result_cbf isa Complex{BigFloat}
            @test abs(result_cbf - result_bf) < big"1e-20"
        end
    end

    @testset "Loop contour upper: integer a-parameter difference, all input types" begin
        a, b, m, n = (0.01, 1.01), (0.5,), 1, 2

        for z in (1.2, 1.5, 2.0)
            # Float64 input.
            result_f64 = meijerg(a, b, m, n, z)
            @test result_f64 isa Number
            @test isfinite(result_f64)

            # Complex{Float64}: complex z, complex result consistent with real result.
            result_c64 = meijerg(a, b, m, n, z + 0im)
            @test result_c64 isa Complex{Float64}
            @test abs(result_c64 - result_f64) < 1e-8 * max(abs(result_f64), 1e-10)

            # BigFloat: higher precision value (possibly complex due to contour path).
            result_bf = meijerg(map(BigFloat, a), map(BigFloat, b), m, n, BigFloat(z))
            @test result_bf isa Number
            @test isfinite(result_bf)

            # Complex{BigFloat}: should match BigFloat-path value.
            result_cbf = meijerg(map(BigFloat, a), map(BigFloat, b), m, n,
                                 BigFloat(z) + zero(BigFloat) * im)
            @test result_cbf isa Complex{BigFloat}
            @test abs(result_cbf - result_bf) < big"1e-20"
        end
    end
end

@testset "Non-finite Slater results fall back to contour" begin
    # Balanced (p == q) confluent geometries are outside the log-ratio jet
    # scope and make raw Slater return Inf/NaN. meijerg should not return that.
    # References: mpmath at 30 digits.
    cases = (
        (((0.5, 0.5), (0, 0), 2, 2, 0.5), 11.649494770833713188),
        (((0.5, 0.5), (0, 1), 2, 2, 0.5), 2.6615984032139112661),
        (((0.5, 1.5), (0.2, 0.1), 1, 2, 2.0), -8.5841561929703619733),
        (((1.329, -0.407), (0.872, 1.872), 2, 1, 0.46946936645877857), 0.037179718359677448818),
    )
    for (args, ref) in cases
        v = meijerg(args...)
        @test isfinite(v)
        @test v ≈ ref rtol=1e-8
    end
end

@testset "High-precision checks" begin
    @testset "BigFloat identities" begin
        setprecision(256) do
            x = big"0.3"
            value = meijerg((), (), (big"0.0",), (), -x)
            # BigFloat exponential result keeps BigFloat type.
            @test value isa BigFloat
            # BigFloat exponential value matches reference.
            @test value ≈ exp(x) rtol=big"1e-70"

            y = big"0.7"
            s = sqrt(big(pi)) * meijerg((), (), (big"0.5",), (big"0.0",), y^2 / 4)
            # BigFloat sine identity is satisfied.
            @test s ≈ sin(y) rtol=big"1e-70"
        end
    end

    @testset "Precision convergence against high-precision internal reference" begin
        a = (big"0.25", big"0.75", big"1.25")
        b = (big"0.1", big"0.6", big"1.6")
        z = big"0.7"

        reference = setprecision(1024) do
            meijerg(a, b, 1, 1, z)
        end

        coarse = setprecision(128) do
            meijerg(a, b, 1, 1, z)
        end

        fine = setprecision(512) do
            meijerg(a, b, 1, 1, z)
        end

        relerr(x, y) = abs(x - y) / max(abs(y), eps(typeof(abs(y))))
        # Higher precision improves relative error.
        @test relerr(fine, reference) < relerr(coarse, reference)
        # Require near-full 512-bit accuracy (eps(512 bits) ≈ 1.5e-154).
        @test relerr(fine, reference) < big"1e-150"
    end
end
