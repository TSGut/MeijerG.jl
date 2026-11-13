@testset "API consistency and structural identities" begin
    @testset "Index form matches split form" begin
        value1 = meijerg((0.25, 1.75), (0.5, 1.25), 1, 1, 0.5)
        value2 = meijerg((0.25,), (1.75,), (0.5,), (1.25,), 0.5)
        # Index and split signatures agree.
        @test value1 ≈ value2 rtol=1e-15
    end

    @testset "Vector and tuple arguments agree" begin
        a = [0.25, 1.75]
        b = [0.5, 1.25]
        # Vector and tuple inputs produce same value.
        @test meijerg(a, b, 1, 1, 0.5) ≈ meijerg(tuple(a...), tuple(b...), 1, 1, 0.5) rtol=1e-15
    end

    @testset "Inversion symmetry" begin
        a = (0.25, 0.75)
        b = (0.1, 0.6)
        z = 0.5
        lhs = meijerg(a, b, 1, 1, z)
        rhs = meijerg(Tuple(1 .- b), Tuple(1 .- a), 1, 1, inv(z))
        # Inversion identity is satisfied.
        @test lhs ≈ rhs rtol=2e-15
    end

    @testset "Stress cases: split/index consistency" begin
        cases = [
            ((0.23, 1.43), (0.19, 1.91), 0.13),
            ((0.31, 1.57), (0.22, 2.01), 0.21),
            ((0.44, 1.68), (0.35, 2.12), 0.29),
            ((0.52, 1.79), (0.41, 2.28), 0.37),
            ((0.63, 1.83), (0.54, 2.35), 0.45),
            ((0.74, 1.94), (0.62, 2.47), 0.53),
            ((0.81, 2.08), (0.73, 2.54), 0.61),
            ((0.92, 2.17), (0.85, 2.68), 0.69),
            ((1.03, 2.26), (0.96, 2.74), 0.77),
            ((1.11, 2.39), (1.08, 2.87), 0.85),
        ]

        for (a, b, z) in cases
            v_index = meijerg(a, b, 1, 1, z)
            v_split = meijerg((a[1],), (a[2],), (b[1],), (b[2],), z)
            # Stress case keeps split/index consistency.
            @test v_index ≈ v_split rtol=1e-15
        end
    end
end

@testset "Return type follows the inputs on every path" begin
    # (label, a, b, m, n, z): one case per evaluation path of meijerg
    cases = (
        ("reduction",           (),         (0, 1//2),     1, 0, 1//2),
        ("slater",              (3//10, 4//5), (1//5, 11//10), 1, 1, 7//10),
        ("near-confluent jet",  (3//10,),   (1//5, 1//5 + 1//10^7), 2, 1, 1//2),
        ("exact-confluent jet", (3//10,),   (1//5, 6//5),  2, 1, 1//2),
        ("contour fallback",    (1//2, 1//2), (0, 1),      2, 2, 1//2),
    )
    for (label, a, b, m, n, z) in cases, T in (Float32, Float64, BigFloat)
        label == "contour fallback" && T === BigFloat && continue
        v = meijerg(T.(a), T.(b), m, n, T(z))
        @test v isa T
        vc = meijerg(T.(a), T.(b), m, n, Complex{T}(z))
        @test vc isa Complex{T}
    end
    # integer input promotes to Float64 (no integer zero from exact cancellations)
    @test meijerg((1, 2), (0, 1), 1, 1, 2) isa Float64
    for f in (meijerg, meijerg_slater, meijerg_contour)
        @test f((-1,), (0, 2), 1, 1, 1) isa Float64
    end
    @test meijerg_laurent((-1,), (0, 1), 2, 1, 1) isa Float64
    @test meijerg_laurent((0.3,), (0.2, 1.2), 2, 1, 0.5) isa Float64
    @test meijerg_contour((0.3,), (0.2, 1.1), 1, 1, 0.7) isa Float64
end
