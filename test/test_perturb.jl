@testset "Perturbative evaluator" begin
    relerr(x, y) = abs(x - y) / max(abs(y), eps(typeof(abs(y))))

    @testset "Mapped-log perturbation" begin
        @testset "Logarithmic confluent mapping" begin
            # log(1+z) via G_{2,2}^{1,2}(z | 1, 1 ; 1, 0) = log(1+z)/z
            z_val = 0.3
            result = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_val)
            expected = log1p(z_val) / z_val
            # Mapped-log value matches the analytic reference.
            @test result ≈ expected atol=1e-12 rtol=1e-12

            z_small = 1e-10
            result_small = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_small)
            expected_small = log1p(z_small) / z_small
            # Near-zero evaluation remains stable against the same reference.
            @test result_small ≈ expected_small atol=1e-10 rtol=1e-10

            z_big = BigFloat("0.3")
            result_big = meijerg((BigFloat(1), BigFloat(1)), (BigFloat(1), BigFloat(0)), 1, 2, z_big)
            expected_big = log1p(z_big) / z_big
            # BigFloat input keeps BigFloat output.
            @test result_big isa BigFloat
            # BigFloat mapped-log value matches the analytic reference.
            @test result_big ≈ expected_big atol=big"1e-30" rtol=big"1e-30"

            z_complex = 0.3 + 0.2im
            result_complex = meijerg((1.0, 1.0), (1.0, 0.0), 1, 2, z_complex)
            expected_complex = log1p(z_complex) / z_complex
            # Complex mapped-log value matches the analytic reference.
            @test result_complex ≈ expected_complex atol=1e-11 rtol=1e-11
        end

        @testset "Confluent expansion consistency" begin
            z_in = 0.5
            z_out = 2.0
            lhs_in = meijerg((1.0, 1.0), (), (1.0,), (0.0,), z_in)
            lhs_out = meijerg((1.0, 1.0), (), (1.0,), (0.0,), z_out)
            # In-domain split-form confluent case matches log reference.
            @test lhs_in ≈ log1p(z_in) / z_in atol=1e-10 rtol=1e-10
            # Out-of-domain split-form confluent case matches log reference.
            @test lhs_out ≈ log1p(z_out) / z_out atol=1e-10 rtol=1e-10
        end

        @testset "Confluent rejection when not applicable" begin
            # Pairing violation (a[j]-b[k] positive integer) remains invalid.
            @test_throws DomainError meijerg((1.0,), (0.0,), 1, 1, 0.5)
        end
    end

    @testset "Lower-mode perturbation stress" begin
        # Lower-mode confluent families stay finite and agree with high-precision reference.
        cases = [
            ((0.35,), (0.10, 1.10, 1.85), 2, 1, 0.20),
            ((0.45,), (0.20, 1.20, 2.05), 2, 1, 0.35),
            ((0.55,), (0.30, 1.30, 2.25), 2, 1, 0.50),
            ((0.75,), (0.50, 1.50, 2.65), 2, 1, 0.90),
        ]

        for (a, b, m, n, z) in cases
            value64 = meijerg(a, b, m, n, z)
            reference = setprecision(512) do
                meijerg(map(BigFloat, a), map(BigFloat, b), m, n, BigFloat(z))
            end

            # Float64 evaluation should remain finite on stress inputs.
            @test isfinite(value64)
            # Float64 result should track high-precision reference closely.
            @test relerr(BigFloat(value64), reference) < big"1e-12"
        end
    end

    @testset "Upper-mode perturbation stress" begin
        # Upper-mode confluent families stay finite and agree with high-precision reference.
        cases = [
            ((1.10, 2.10, 3.00), (0.20,), 1, 2, 1.60),
            ((1.50, 2.50, 3.60), (0.35,), 1, 2, 2.70),
            ((1.10, 2.10, 2.90), (0.20,), 1, 2, 1.20),
            ((1.40, 2.40, 3.50), (0.50,), 1, 2, 2.00),
            ((1.50, 2.50, 3.70), (0.60,), 1, 2, 2.40),
        ]

        for (a, b, m, n, z) in cases
            value64 = meijerg(a, b, m, n, z)
            reference = setprecision(512) do
                meijerg(map(BigFloat, a), map(BigFloat, b), m, n, BigFloat(z))
            end

            # Float64 evaluation should remain finite on stress inputs.
            @test isfinite(value64)
            # Float64 result should track high-precision reference closely.
            @test relerr(BigFloat(value64), reference) < big"1e-12"
        end
    end
end