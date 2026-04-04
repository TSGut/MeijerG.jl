using Test
using MeijerG

include("test_utilities.jl")

@testset "Main evaluator vs contour evaluator" begin
    relerr(x, y) = abs(x - y) / max(abs(y), eps(typeof(abs(y))), 1e-16)

    @testset "Non-confluent lower mode (Float64, ComplexF64)" begin
        a, b, m, n = (0.35,), (0.2, 1.1), 1, 1

        for z in sample_real_points((0.2, 0.6, 0.9); seed=401, lo=0.15, hi=1.0, count=2)
            main_val = meijerg(a, b, m, n, z)
            contour_val = meijerg_contour(a, b, m, n, z;
                                          contour_case=:auto,
                                          rtol=1e-9,
                                          atol=1e-11,
                                          max_truncation_steps=6)
            @test relerr(main_val, contour_val) < 1e-7
        end

        for z in sample_complex_points((0.3 + 0.2im, 0.8 - 0.15im);
                                       seed=402,
                                       real_lo=0.2,
                                       real_hi=0.9,
                                       imag_lo=-0.3,
                                       imag_hi=0.3,
                                       count=2)
            main_val = meijerg(a, b, m, n, z)
            contour_val = meijerg_contour(a, b, m, n, z;
                                          contour_case=:auto,
                                          rtol=1e-9,
                                          atol=1e-11,
                                          max_truncation_steps=6)
            @test relerr(main_val, contour_val) < 1e-6
        end
    end

    @testset "Non-confluent upper mode (Float64, ComplexF64)" begin
        a, b, m, n = (0.4, 1.3), (0.2,), 1, 1

        for z in sample_real_points((1.2, 1.8, 2.4); seed=403, lo=1.1, hi=2.6, count=2)
            main_val = meijerg(a, b, m, n, z)
            contour_val = meijerg_contour(a, b, m, n, z;
                                          contour_case=:auto,
                                          rtol=1e-9,
                                          atol=1e-11,
                                          max_truncation_steps=6)
            @test relerr(main_val, contour_val) < 1e-7
        end

        for z in sample_complex_points((1.4 + 0.25im, 2.0 - 0.2im);
                                       seed=404,
                                       real_lo=1.2,
                                       real_hi=2.3,
                                       imag_lo=-0.3,
                                       imag_hi=0.3,
                                       count=2)
            main_val = meijerg(a, b, m, n, z)
            contour_val = meijerg_contour(a, b, m, n, z;
                                          contour_case=:auto,
                                          rtol=1e-9,
                                          atol=1e-11,
                                          max_truncation_steps=6)
            @test relerr(main_val, contour_val) < 1e-6
        end
    end

    @testset "Confluent lower mode" begin
        # Integer-separated active lower parameters force confluent behavior.
        a, b, m, n = (0.5,), (0.01, 1.01), 2, 1

        for z in (0.3, 0.8)
            main_val = meijerg(a, b, m, n, z)
            contour_val = meijerg_contour(a, b, m, n, z;
                                          contour_case=:auto,
                                          rtol=1e-8,
                                          atol=1e-10,
                                          max_truncation_steps=6)
            @test relerr(main_val, contour_val) < 2e-5
        end

        z_bf = BigFloat("0.5")
        a_bf = map(BigFloat, a)
        b_bf = map(BigFloat, b)
        main_bf = meijerg(a_bf, b_bf, m, n, z_bf)
        contour_bf = meijerg_contour(a_bf, b_bf, m, n, z_bf;
                                     contour_case=:auto,
                                     rtol=big"1e-24",
                                     atol=big"1e-30",
                                     max_truncation_steps=6)
        @test relerr(main_bf, contour_bf) < big"1e-10"
    end

    @testset "Confluent upper mode" begin
        # Integer-separated active upper parameters force confluent behavior.
        a, b, m, n = (0.01, 1.01), (0.5,), 1, 2

        for z in (1.3, 2.0)
            main_val = meijerg(a, b, m, n, z)
            contour_val = meijerg_contour(a, b, m, n, z;
                                          contour_case=:auto,
                                          rtol=1e-8,
                                          atol=1e-10,
                                          max_truncation_steps=6)
            @test relerr(main_val, contour_val) < 2e-5
        end

        z_bf = BigFloat("1.6")
        a_bf = map(BigFloat, a)
        b_bf = map(BigFloat, b)
        main_bf = meijerg(a_bf, b_bf, m, n, z_bf)
        contour_bf = meijerg_contour(a_bf, b_bf, m, n, z_bf;
                                     contour_case=:auto,
                                     rtol=big"1e-24",
                                     atol=big"1e-30",
                                     max_truncation_steps=6)
        @test relerr(main_bf, contour_bf) < big"1e-10"
    end
end
