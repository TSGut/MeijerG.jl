using Test
using MeijerG


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

    @testset "Explicit loop contour residue continuations" begin
        # Regression cases for DLMF loop cases ii/iii where no vertical strip
        # exists. These are checked against raw Slater when nonconfluent and
        # against mpmath-supported reference values at exact confluence.
        a_zero = (0.19 + 0.31im, 0.74 - 0.27im, 1.36 + 0.16im)
        b_zero = (-0.42 + 0.18im, 0.58 + 0.18im, 0.69 + 0.24im,
                  1.28 - 0.11im, 1.77 + 0.21im)
        z_zero = 3e-5 + 0im

        b_zero_perturbed = (-0.42 + 0.18im, 0.58 + 0.18im + 1e-4,
                            0.69 + 0.24im, 1.28 - 0.11im, 1.77 + 0.21im)
        lower_raw = meijerg_slater(a_zero, b_zero_perturbed, 2, 2, z_zero)
        lower_loop = meijerg_contour(a_zero, b_zero_perturbed, 2, 2, z_zero;
                                     contour_case=:ii, rtol=1e-9, atol=1e-11)
        @test relerr(lower_loop, lower_raw) < 1e-12

        lower_ref = -2.821531876170197 - 1.109576498726359im
        lower_confluent = meijerg_contour(a_zero, b_zero, 2, 2, z_zero;
                                          contour_case=:ii, rtol=1e-8, atol=1e-10)
        @test relerr(lower_confluent, lower_ref) < 1e-8

        a_upper = (-0.22 + 0.18im, 1.78 + 0.18im,
                   2.41 - 0.22im, 3.03 + 0.31im)
        b_upper = (0.26 + 0.31im, 1.12 - 0.16im)
        z_upper = 50.0 + 0im

        a_upper_perturbed = (-0.22 + 0.18im, 1.78 + 0.18im + 1e-4,
                             2.41 - 0.22im, 3.03 + 0.31im)
        upper_raw = meijerg_slater(a_upper_perturbed, b_upper, 1, 2, z_upper)
        upper_loop = meijerg_contour(a_upper_perturbed, b_upper, 1, 2, z_upper;
                                     contour_case=:iii, rtol=1e-9, atol=1e-11)
        @test relerr(upper_loop, upper_raw) < 1e-12

        upper_ref = -29.171421299993224 - 56.326573936372796im
        upper_confluent = meijerg_contour(a_upper, b_upper, 1, 2, z_upper;
                                          contour_case=:iii, rtol=1e-8, atol=1e-10)
        @test relerr(upper_confluent, upper_ref) < 1e-8
    end

    @testset "Loop contour circles do not overlap nearby nonmember poles" begin
        a = (0.20 + 0.25im, 0.78 - 0.15im, 1.31 + 0.19im)
        b1 = 0.52 + 0.14im
        btail = (-0.47 + 0.28im, 1.05 - 0.21im, 1.67 + 0.16im)
        z = 2.5e-5 * cis(-0.09)

        # The following test is a regression test for a bug in older versions:
        # These separations previously sat between the clustering tolerance and
        # the old absolute radius floor, so separate circles could enclose the
        # neighboring active pole and return stable but wrong values.
        for epsv in (1e-5, 3e-6, 1e-6, 3e-8)
            b = (b1, b1 + epsv, btail...)
            reference = MeijerG._lower_logratio_jet_series(a, b, 3, 2, z, 0, epsv)
            @test reference.status === :ok

            contour = meijerg_contour(a, b, 3, 2, z;
                                      contour_case=:ii,
                                      rtol=1e-9,
                                      atol=1e-11,
                                      max_truncation_steps=8)
            @test relerr(contour, reference.value) < 1e-9
        end
    end

    @testset "Vertical contour only where its integral converges" begin
        # (m, n, p, q) = (4, 2, 5, 8): a separating strip exists, but 2(m + n) < p + q, so the
        # vertical integrand grows along the line. Case (i) is not selected and the loop contour
        # evaluates the integral.
        a = (0.18 + 0.24im, 0.83 - 0.16im, 1.24 + 0.21im, -0.33 + 0.28im, 1.71 - 0.18im)
        b = (0.42 + 0.13im, 0.421 + 0.13im, 0.422 + 0.13im, 1.38 - 0.29im,
             -0.52 + 0.31im, 1.09 - 0.23im, 1.87 + 0.19im, 2.41 - 0.14im)
        @test MeijerG._has_vertical_strip(a, b, 4, 2)
        @test !MeijerG._vertical_contour_converges(a, b, 4, 2, 200.0)
        @test MeijerG._resolve_contour_case(:auto, a, b, 4, 2, 200.0) === :ii
        value = meijerg_contour(a, b, 4, 2, 200.0)
        @test isfinite(value)
        @test relerr(value, meijerg(a, b, 4, 2, 200.0)) < 1e-5

        # The mirrored upper geometry (3, 4, 8, 5) converges for |arg z| < π/2 only.
        au = (0.51 + 0.17im, 0.511 + 0.17im, 0.512 + 0.17im, -0.42 + 0.27im,
              1.19 - 0.22im, -0.74 + 0.33im, 1.83 + 0.2im, 2.27 - 0.15im)
        bu = (0.28 - 0.18im, 0.94 + 0.25im, 1.46 - 0.16im, -0.37 + 0.3im, 1.88 - 0.2im)
        @test MeijerG._resolve_contour_case(:auto, au, bu, 3, 4, 0.01) === :i
        @test MeijerG._resolve_contour_case(:auto, au, bu, 3, 4, 0.01 * cis(0.3pi)) === :i
        @test MeijerG._resolve_contour_case(:auto, au, bu, 3, 4, 0.01 * cis(0.9pi)) === :iii

        # An explicit vertical request on the divergent geometry still throws.
        @test_throws DomainError meijerg_contour(a, b, 4, 2, 200.0; contour_case=:i)
    end
end
