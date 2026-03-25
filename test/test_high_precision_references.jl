@testset "High-precision and external-reference checks" begin
    @testset "BigFloat identities" begin
        setprecision(256) do
            x = big"0.3"
            value = meijerg((), (), (big"0.0",), (), -x)
            # BigFloat exponential result keeps BigFloat type.
            @test value isa BigFloat
            # BigFloat exponential value matches reference.
            @test value ≈ exp(x) atol=big"1e-60" rtol=big"1e-60"

            y = big"0.7"
            s = sqrt(big(pi)) * meijerg((), (), (big"0.5",), (big"0.0",), y^2 / 4)
            # BigFloat sine identity is satisfied.
            @test s ≈ sin(y) atol=big"1e-50" rtol=big"1e-50"
        end
    end

    @testset "Reference values documented by mpmath" begin
        ref1 = 271.46290321152464592 - 703.03330399954820169im
        val1 = meijerg((), (1 + 1im,), (1,), (1,), 3 + 4im)
        # Complex reference 1 matches mpmath.
        @test val1 ≈ ref1 atol=1e-10 rtol=1e-10

        ref2 = -1.338096165935754898687431
        val2 = meijerg((-3.0,), (-0.5,), (-1.0,), (-2.5,), -0.5)
        # Real reference 2 matches mpmath.
        @test val2 ≈ ref2 atol=1e-12 rtol=1e-12

        ref3 = -(pi + 4) / (4pi)
        val3 = meijerg((-3.0,), (-0.5,), (-1.0,), (-2.5,), -1.0)
        # Real reference 3 matches closed form.
        @test val3 ≈ ref3 atol=1e-12 rtol=1e-12
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
        # Fine precision reaches target tolerance.
        @test relerr(fine, reference) < big"1e-40"
    end
end
