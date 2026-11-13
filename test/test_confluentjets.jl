module TestConfluentJets

using MeijerG.ConfluentJets
using Gamma
using Test

@testset "ConfluentJets" begin

@testset "Jet arithmetic" begin
    x = affine_jet(2.0, 1.0, 5)
    @test coefficients(log(x))[1:4] ≈
          [log(2.0), 1 / 2, -1 / 8, 1 / 24]
    @test coefficients(exp(log(x))) ≈ coefficients(x)
    @test coefficients(x * inv(x)) ≈ coefficients(one(x))
    @test coefficient(compose_scale(x, 3), 1) == 3
    @test coefficients(affine_log_jet(2.0, 1.0, 5))[1:4] ==
          [log(2.0), 1 / 2, -1 / 8, 1 / 24]
end

@testset "Gamma-family jets" begin
    x = affine_jet(0.7 + 0.2im, 1.0 + 0im, 9)
    lg = loggamma(x)
    @test coefficient(lg, 0) == loggamma(0.7 + 0.2im)
    @test coefficient(lg, 1) == digamma(0.7 + 0.2im)
    @test coefficients(gamma(x)) ≈ coefficients(exp(lg))
    @test coefficients(digamma(x))[1] == digamma(0.7 + 0.2im)
    @test polygamma(1, x) isa Jet

    pole = affine_jet(-2.0 + 0im, 1.0 + 0im, 8)
    reciprocal = reciprocal_gamma(pole)
    @test reciprocal isa Jet
    @test iszero(coefficient(reciprocal, 0))
end

@testset "Opitz divided differences" begin
    logjet = log(affine_jet(2.0 + 0im, 1.0 + 0im, 10))
    nodes = ComplexF64[-1e-4, 2e-4, 4e-4]
    value = exp_divided_difference(logjet, nodes)
    @test value ≈ 0.0 atol=1e-12

    factor = affine_jet(1.0 + 0im, 0.5 + 0im, 10)
    product_value =
        exp_product_divided_difference(factor, logjet, nodes)
    @test isfinite(product_value)
end

@testset "Jet-valued 2F1 series" begin
    N = 2
    a = affine_jet(1.0, 1.0, N)
    b = constant_jet(2.0, N)
    c = constant_jet(3.0, N)
    z = 0.2
    term = one(a)
    total = zero(a)
    derivative_reference = 0.0
    scalar_term = 1.0
    harmonic_increment = 0.0
    for k in 0:39
        total += term
        derivative_reference += scalar_term * harmonic_increment
        term *= ((a + k) * (b + k) / (c + k)) * (z / (k + 1))
        scalar_term *= ((1.0 + k) * (2.0 + k) / (3.0 + k)) *
                       (z / (k + 1))
        harmonic_increment += inv(1.0 + k)
    end
    @test coefficient(total, 1) ≈ derivative_reference rtol=2e-14
end

@testset "Complex BigFloat" begin
    setprecision(512) do
        z = affine_jet(Complex{BigFloat}(big"-3.4", big"0.2"),
                       Complex{BigFloat}(1), 19)
        value = loggamma(z)
        @test value isa Jet{Complex{BigFloat}}
        @test all(isfinite, real.(coefficients(value)))
        @test all(isfinite, imag.(coefficients(value)))
    end
end
end

end # module TestConfluentJets
