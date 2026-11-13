using Test
using MeijerG

relerr_higher(x, y) = abs(x - y) / max(abs(y), eps(Float64))

@testset "Higher-order confluence internals" begin
    @testset "standalone Laurent exact multiplicity" begin
        b1 = 0.5
        cases = (
            (id=:lower_r3, a=(0.3, 0.7, 1.4),
             b=(b1, b1 + 1, b1 + 2, -0.4, 1.1, 1.7),
             m=3, n=1, z=3e-4 + 0im,
             expected=-0.0002611044611187978 + 0im, r=3),
            (id=:lower_r4, a=(0.3, 0.7, 1.4, 0.9),
             b=(b1, b1 + 1, b1 + 2, b1 + 3, -0.4, 1.1, 1.7, 2.1),
             m=4, n=1, z=3e-4 + 0im,
             expected=6.366499430334363e-5 + 0im, r=4),
            (id=:upper_r3, a=(0.4, 1.4, 2.4, 0.7, 1.1, 1.6),
             b=(0.3, -0.5, 0.9),
             m=1, n=3, z=50.0 + 0im,
             expected=-62.4043178028642 + 0im, r=3),
            (id=:upper_r4, a=(0.4, 1.4, 2.4, 3.4, 0.7, 1.1, 1.6, 2.1),
             b=(0.3, -0.5, 0.9, 1.3),
             m=1, n=4, z=50.0 + 0im,
             expected=110.54693696056081 + 0im, r=4),
        )

        for c in cases
            r = MeijerG._try_meijerg_slater_higher_laurent_exact(c.a, c.b, c.m,
                                                                 c.n, c.z)
            @test r.status === :ok
            @test r.max_multiplicity == c.r
            @test relerr_higher(r.value, c.expected) < 5e-11
            @test relerr_higher(meijerg_laurent(c.a, c.b, c.m, c.n, c.z),
                                c.expected) < 5e-11
            raw_status = try
                v = meijerg_slater(c.a, c.b, c.m, c.n, c.z)
                isfinite(real(v)) && isfinite(imag(v))
            catch
                false
            end
            @test !raw_status
        end
    end

    @testset "standalone Laurent unsupported cases" begin
        b1 = 0.5
        a = (0.3, 0.7, 1.4)
        b5 = (ntuple(k -> b1 + (k - 1), 5)..., -0.4)
        r = MeijerG._try_meijerg_slater_higher_laurent_exact(a, b5, 5, 1,
                                                            3e-4 + 0im)
        @test r.status === :unsupported_multiplicity
        @test r.max_multiplicity == 5
        @test_throws ArgumentError meijerg_laurent(a, b5, 5, 1, 3e-4 + 0im)

        @test_throws ArgumentError meijerg_laurent((0.3, 0.7),
                                                  (0.2, 0.2 + 1e-4),
                                                  2, 1, 0.5)
    end

    @testset "higher log-ratio exact and near" begin
        b1 = 0.5
        exact_cases = (
            (a=(0.3, 0.7, 1.4),
             b=(b1, b1 + 1, b1 + 2, -0.4, 1.1),
             m=3, n=1, z=3e-4 + 0im, r=3),
            (a=(0.3, 0.7, 1.4, 0.9),
             b=(b1, b1 + 1, b1 + 2, b1 + 3, -0.4, 1.1),
             m=4, n=1, z=3e-4 + 0im, r=4),
            (a=(0.4, 1.4, 2.4, 0.7, 1.1, 1.6),
             b=(0.3, -0.5, 0.9),
             m=1, n=3, z=50.0 + 0im, r=3),
            (a=(0.4, 1.4, 2.4, 3.4, 0.7, 1.1, 1.6, 2.1),
             b=(0.3, -0.5, 0.9, 1.3),
             m=1, n=4, z=50.0 + 0im, r=4),
        )
        for c in exact_cases
            jet = MeijerG._try_meijerg_slater_higher_logratio_jet(c.a, c.b,
                                                                  c.m, c.n, c.z)
            laurent = meijerg_laurent(c.a, c.b, c.m, c.n, c.z)
            @test jet.status === :ok
            @test jet.selected
            @test relerr_higher(jet.value, laurent) < 5e-11
            @test relerr_higher(meijerg(c.a, c.b, c.m, c.n, c.z), jet.value) < 5e-11
        end

        near_cases = (
            (a=(0.3, 0.7, 1.4),
             b=(b1, b1 + 1 + 1e-6, b1 + 2 + 2e-6, -0.4, 1.1),
             m=3, n=1, z=3e-4 + 0im,
             expected=0.001519363169323765 + 0im, r=3),
            (a=(0.3, 0.7, 1.4, 0.9),
             b=(b1, b1 + 1 + 1e-6, b1 + 2 + 2e-6,
                b1 + 3 + 3e-6, -0.4, 1.1),
             m=4, n=1, z=3e-4 + 0im,
             expected=0.0013701208761540918 + 0im, r=4),
            (a=(0.4, 1.4 + 1e-6, 2.4 + 2e-6, 0.7, 1.1, 1.6),
             b=(0.3, -0.5, 0.9),
             m=1, n=3, z=50.0 + 0im,
             expected=-62.40255566116208 + 0im, r=3),
            (a=(0.4, 1.4 + 1e-6, 2.4 + 2e-6, 3.4 + 3e-6,
                0.7, 1.1, 1.6, 2.1),
             b=(0.3, -0.5, 0.9, 1.3),
             m=1, n=4, z=50.0 + 0im,
             expected=110.54106437019861 + 0im, r=4),
        )
        for c in near_cases
            jet = MeijerG._try_meijerg_slater_higher_logratio_jet(c.a, c.b,
                                                                  c.m, c.n, c.z)
            @test jet.status === :ok
            @test jet.selected
            @test relerr_higher(jet.value, c.expected) < 5e-11
            @test relerr_higher(meijerg(c.a, c.b, c.m, c.n, c.z), jet.value) < 5e-11
        end

        # A tightly spaced upper r=4 cluster remains well conditioned in the
        # Opitz representation all the way down to the Float64 snap boundary.
        tight_epsilon = 1e-13
        tight_a = (0.4, 1.4 + tight_epsilon, 2.4 + 2tight_epsilon,
                   3.4 + 3tight_epsilon, 0.7, 1.1, 1.6, 2.1)
        tight_b = (0.3, -0.5, 0.9, 1.3)
        tight_jet = MeijerG._try_meijerg_slater_higher_logratio_jet(
            tight_a, tight_b, 1, 4, 50.0)
        @test tight_jet.status === :ok
        @test tight_jet.selected
        @test relerr_higher(tight_jet.value, 110.54693695997373) < 1e-13

        # Keep a nearby reciprocal-gamma zero outside the logarithmic jet.
        # This is the real r=4, L=0 tube case that previously stalled at 2e-6.
        a = (0.291, 0.963, 1.355, -0.23299999999999998, 1.809)
        b = (0.523, 0.526, 0.529, 0.532, -0.421, 1.161,
             1.9829999999999999, 2.537)
        jet = MeijerG._try_meijerg_slater_higher_logratio_jet(a, b, 4, 2,
                                                              2.73e-5)
        @test jet.status === :ok
        @test jet.selected
        @test relerr_higher(jet.value, 0.008545740357663327) < 5e-13

        b5 = (ntuple(k -> b1 + (k - 1) + (k - 1) * 1e-6, 5)..., -0.4)
        jet = MeijerG._try_meijerg_slater_higher_logratio_jet((0.3, 0.7, 1.4),
                                                              b5, 5, 1,
                                                              3e-4 + 0im)
        @test jet.status === :unsupported_multiplicity
        @test !jet.selected
    end
end

@testset "BigFloat exact higher-order clusters use the Laurent evaluator" begin
    # The higher-order jet is Float32/Float64-only. BigFloat exact r = 3 clusters must
    # still reach full precision (references: mpmath at 100 digits).
    cases = (
        ((big"0.3",), (big"0.2", big"1.2", big"2.2"), 3, 1, big"0.5",
         big"0.58939848492806939162743636025093643817672788820373760314859353153001072947725659"),
        # needs more than the default 40 pole levels at 256 bits
        ((big"0.3", big"0.65"), (big"0.2", big"1.2", big"2.2"), 3, 2, big"0.7",
         big"1.21348917696077602044028550139719466358965720112035777984660710550570143105491142"),
    )
    for (a, b, m, n, z, ref) in cases
        jet = MeijerG._try_meijerg_slater_higher_logratio_jet(a, b, m, n, z)
        @test jet.status === :unsupported_precision
        v = meijerg(a, b, m, n, z)
        @test v isa BigFloat
        @test abs(v - ref) <= 1e3 * eps(BigFloat) * abs(ref)
        @test meijerg_laurent(a, b, m, n, z) ≈ ref rtol = 1e3 * eps(BigFloat)
    end
end
