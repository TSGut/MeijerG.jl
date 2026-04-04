include("test_utilities.jl")

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

    @testset "Relevant example quality near |z| = 1" begin
        a_str = ("0.45", "0.95", "1.55", "2.05")
        b_str = ("0.12", "0.72", "1.22", "1.82")
        a64 = Float64.(parse.(Ref(BigFloat), a_str))
        b64 = Float64.(parse.(Ref(BigFloat), b_str))
        z_literals = ("0.98", "0.99", "0.995", "0.9975", "0.999", "0.9995",
                      "1.0005", "1.001", "1.0025", "1.005", "1.01", "1.02")

        for z_literal in z_literals
            z64 = Float64(parse(BigFloat, z_literal))
            ref = setprecision(1024) do
                a_ref = parse.(Ref(BigFloat), a_str)
                b_ref = parse.(Ref(BigFloat), b_str)
                meijerg_slater(a_ref, b_ref, 1, 1, parse(BigFloat, z_literal))
            end

            guarded = meijerg_slater(a64, b64, 1, 1, z64)
            err_guarded = abs(BigFloat(guarded) - ref) / max(abs(ref), eps(BigFloat))
            @test err_guarded < 1e-12

            direct_big = meijerg_slater(BigFloat.(a64), BigFloat.(b64), 1, 1, BigFloat(z64))
            agree_err = abs(BigFloat(guarded) - direct_big) / max(abs(direct_big), eps(BigFloat))
            # Boundary guard should match explicit BigFloat Slater to rounding accuracy.
            @test agree_err < 1e-12
        end
    end

    @testset "Explicit complex branch selection" begin
        # Real input that would require a complex branch should throw,
        # mirroring Julia's standard real/complex behavior.
        @test_throws DomainError meijerg((), (1, 0), 2, 0, -1.0)
        # Explicit complex input should evaluate on the complex branch.
        @test meijerg((), (1, 0), 2, 0, -1.0 + 0im) isa Complex
    end

    @testset "A compact support test" begin
        # G_{1,1}^{1,0}(x | a ; b) = x^b (1-x)^(a-b-1) / Γ(a-b) on 0 < x < 1,
        # and vanishes for real x > 1 when Re(a) > Re(b).
        a, b = 2.5, 0.5

        for x in sample_real_points((0.25, 0.75); seed=201, lo=0.05, hi=0.95)
            expected = x^b * (1 - x)^(a - b - 1) / gamma(a - b)
            @test meijerg((a,), (b,), 1, 0, x) ≈ expected rtol=1e-14
            @test meijerg_slater((a,), (b,), 1, 0, x) ≈ expected rtol=1e-14
        end

        for x in sample_real_points((1.1, 1.5, 2.0); seed=202, lo=1.01, hi=2.5)
            @test iszero(meijerg((a,), (b,), 1, 0, x))
            @test iszero(meijerg_slater((a,), (b,), 1, 0, x))
        end

        ab, bb = big"2.5", big"0.5"
        for xb in BigFloat.(sample_real_points((1.1, 1.5, 2.0); seed=203, lo=1.01, hi=2.5))
            @test iszero(meijerg((ab,), (bb,), 1, 0, xb))
            @test iszero(meijerg_slater((ab,), (bb,), 1, 0, xb))
        end
    end
end
