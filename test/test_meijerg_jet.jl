using Test
using MeijerG

# Public standalone jet evaluator and its working-precision modes.
# References: mpmath meijerg at 130 digits (branch-matched series), as in
# test_logratio_jet.jl and test_extended_lanes.jl.

@testset "meijerg_jet" begin
    exact_cases = [
        ((0.3,), (0.2, 1.2), 2, 1, 20.0,
         big"0.1166995301016219979085603590513498095972758108483352635445102661674967803911732141289087103293602461"),
        ((), (0.2, 1.2, 0.7), 2, 0, 80.0,
         big"-0.002677243745218326127525353237575407726442827858963248195349324412420109013604634215462975032153075858"),
    ]

    @testset "off the convergent regime" begin
        for (a, b, m, n, z, ref) in exact_cases
            vauto = meijerg_jet(a, b, m, n, z)
            @test vauto isa Float64
            @test vauto == meijerg(a, b, m, n, z)
            @test abs(vauto - ref) < 4 * eps() * abs(ref)
            @test meijerg_jet(a, b, m, n, z; working_precision=:extended) == vauto
            # The uncertified input-precision sum loses digits to cancellation.
            vin = meijerg_jet(a, b, m, n, z; working_precision=:input)
            @test isfinite(vin)
            @test abs(vin - ref) > 1e-13 * abs(ref)
        end
        # Finite separation 1e-4, lower branch at z = 20.
        a, b, m, n, z = (0.3,), (0.2, 1.2001), 2, 1, 20.0
        ref = big"0.116703223857659018894238408036303642577908279"
        @test abs(meijerg_jet(a, b, m, n, z) - ref) < 4 * eps() * ref
        @test abs(meijerg_jet(a, b, m, n, z; working_precision=:extended) - ref) < 4 * eps() * ref
        @test abs(meijerg_jet(a, b, m, n, z; working_precision=:input) - ref) > 1e-8 * ref
        # Higher-order (r = 3) cluster.
        a, b, m, n, z = (0.3,), (0.2, 1.2, 2.2), 3, 1, 5.0
        v = meijerg_jet(a, b, m, n, z)
        @test v == meijerg(a, b, m, n, z)
        @test meijerg_jet(a, b, m, n, z; working_precision=:extended) == v
        @test isfinite(meijerg_jet(a, b, m, n, z; working_precision=:input))
    end

    @testset "inside the convergent regime" begin
        # :auto and :input keep the Float64 series, :extended sums in double-double.
        a, b, m, n = (0.3,), (0.2, 1.2), 2, 1
        for z in (0.05, 0.4, 0.9)
            ref = setprecision(BigFloat, 256) do
                meijerg((big"0.3",), (big"0.2", big"1.2"), m, n, big(z))
            end
            vauto = meijerg_jet(a, b, m, n, z)
            @test meijerg_jet(a, b, m, n, z; working_precision=:input) == vauto
            @test abs(vauto - ref) < 1e-14 * abs(ref)
            vext = meijerg_jet(a, b, m, n, z; working_precision=:extended)
            @test abs(vext - ref) < 2 * eps() * abs(ref)
        end
        # Upper exact pair.
        for z in (2.5, 20.0)
            ref = setprecision(BigFloat, 256) do
                meijerg((big"0.3", big"-0.7"), (big"0.2",), 1, 2, big(z))
            end
            vext = meijerg_jet((0.3, -0.7), (0.2,), 1, 2, z; working_precision=:extended)
            @test abs(vext - ref) < 2 * eps() * abs(ref)
        end
        r = meijerg_jet((0.3,), (0.2, 1.2, 2.2), 3, 1, 0.5)
        @test abs(meijerg_jet((0.3,), (0.2, 1.2, 2.2), 3, 1, 0.5;
                              working_precision=:extended) - r) < 1e-14 * abs(r)
    end

    @testset "BigFloat and Float32" begin
        a, b, m, n, z, ref = exact_cases[1]
        big_in(x) = parse(BigFloat, string(x))
        setprecision(BigFloat, 256) do
            ab, bb, zb = map(big_in, a), map(big_in, b), big_in(z)
            for wp in (:auto, :extended)
                v = meijerg_jet(ab, bb, m, n, zb; working_precision=wp)
                @test precision(v) == 256
                @test abs(v - ref) < big"1e-70" * abs(ref)
            end
            # Inside the regime :extended doubles the precision, :auto keeps it.
            r = setprecision(BigFloat, 1024) do
                meijerg(map(big_in, a), map(big_in, b), m, n, big_in(0.5))
            end
            vauto = meijerg_jet(ab, bb, m, n, big_in(0.5))
            vext = meijerg_jet(ab, bb, m, n, big_in(0.5); working_precision=:extended)
            @test precision(vext) == 256
            @test abs(vext - r) <= abs(vauto - r)
            @test abs(vext - r) < 4 * eps(BigFloat) * abs(r)
        end
        for wp in (:auto, :extended)
            @test meijerg_jet((0.3f0,), (0.2f0, 1.2001f0), 2, 1, 20.0f0;
                              working_precision=wp) isa Float32
        end
        # Summed in Float32 alone, the series off its regime is not finite.
        @test_throws ErrorException meijerg_jet((0.3f0,), (0.2f0, 1.2001f0), 2, 1, 20.0f0;
                                                working_precision=:input)
        @test meijerg_jet((0.3f0,), (0.2f0, 1.2001f0), 2, 1, 0.5f0;
                          working_precision=:input) isa Float32
        @test abs(meijerg_jet((0.3f0,), (0.2f0, 1.2001f0), 2, 1, 20.0f0) -
                  meijerg((0.3f0,), (0.2f0, 1.2001f0), 2, 1, 20.0f0)) <
              4 * eps(Float32) * abs(meijerg((0.3f0,), (0.2f0, 1.2001f0), 2, 1, 20.0f0))
    end

    @testset "calling conventions and errors" begin
        @test meijerg_jet([0.3], [], [0.2, 1.2], [], 0.4) ==
              meijerg_jet((0.3,), (0.2, 1.2), 2, 1, 0.4)
        @test meijerg_jet([0.3], [0.2, 1.2], 2, 1, 0.4) ==
              meijerg_jet((0.3,), (0.2, 1.2), 2, 1, 0.4)
        @test meijerg_jet((0.3,), (0.2, 1.2), 2, 1, 0.4 + 0im) isa ComplexF64
        @test_throws ArgumentError meijerg_jet((0.3,), (0.2, 0.7), 2, 1, 0.5)
        @test_throws ArgumentError meijerg_jet((0.3, 0.4), (0.2, 1.2), 1, 1, 0.5)
        @test_throws ArgumentError meijerg_jet((big"0.3",), (big"0.2", big"1.2", big"2.2"),
                                               3, 1, big"5.0")
        @test_throws ArgumentError meijerg_jet((0.3,), (0.2, 1.2), 2, 1, 0.5;
                                               working_precision=:float64)
        @test_throws DomainError meijerg_jet((0.3,), (0.2, 1.2), 2, 1, -0.5)
        @test_throws DomainError meijerg_jet((0.3,), (0.2, 1.2), 2, 1, 0.0)
    end
end
