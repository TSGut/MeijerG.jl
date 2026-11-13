using Test
using MeijerG

# Doubled-precision lanes: the private double-double type with the generic jet
# machinery, and the finite-separation, higher-order and separated-Slater lanes.
# References: mpmath meijerg at 60 digits (branch-matched series), with the
# Float64 parameter values spelled out as decimals.

@testset "Extended-precision lanes" begin
    DD = MeijerG._DD
    CDD = MeijerG._CDD
    tobig(z::CDD) = Complex{BigFloat}(BigFloat(real(z)), BigFloat(imag(z)))
    relerr(x::CDD, r) = Float64(abs(tobig(x) - r) / abs(r))
    lg = MeijerG.loggamma
    psi = MeijerG.digamma

    @testset "double-double primitives and generic jets" begin
        setprecision(BigFloat, 256) do
            for (x, y) in ((0.3 + 0.7im, -2.25 + 0.1im), (12.5 - 3.0im, 0.5im))
                X = CDD(x); Y = CDD(y)
                bx = Complex{BigFloat}(x); by = Complex{BigFloat}(y)
                @test relerr(X * Y, bx * by) < 1e-30
                @test relerr(X / Y, bx / by) < 1e-30
                @test relerr(exp(X), exp(bx)) < 1e-30
                @test relerr(log(X), log(bx)) < 1e-30
                @test relerr(log1p(Y / 16), log1p(by / 16)) < 1e-30
                @test relerr(expm1(Y / 16), expm1(by / 16)) < 1e-30
                @test relerr(psi(Y), psi(by)) < 1e-30
                @test relerr(sinpi(Y), sinpi(by)) < 1e-30
                @test relerr(X^-3, bx^-3) < 1e-30
                # log-gamma is formed modulo 2πi
                d = tobig(lg(X)) - lg(bx)
                @test abs(d - 2im * big(pi) * round(imag(d) / (2 * big(pi)))) < 1e-29
            end
            for x in (0.3 + 0.2im, -2.6 + 0.0im, 7.0 + 1.0im)
                jd = lg(MeijerG.ConfluentJets.affine_jet(CDD(x), CDD(1), 8))
                jb = lg(MeijerG.ConfluentJets.affine_jet(Complex{BigFloat}(x),
                                                         Complex{BigFloat}(1), 8))
                @test maximum(Float64(abs(tobig(jd.coeffs[k]) - jb.coeffs[k]) /
                                      abs(jb.coeffs[k])) for k in 2:9) < 1e-28
            end
            for (x, e) in ((0.3 + 0.2im, 1e-4), (-2.6 + 0.0im, 1e-9))
                dd = MeijerG._dd_loggamma_divdiff(CDD(x), CDD(e), 1)
                bx = Complex{BigFloat}(x); be = Complex{BigFloat}(e)
                @test relerr(dd, (lg(bx + be) - lg(bx)) / be) < 1e-28
            end
        end
        # Opitz divided difference: the jet route (used for double-double) agrees
        # with the matrix-exponential route (Float64).
        L = MeijerG.ConfluentJets.Jet(ComplexF64[0.3 + 0.1im, 1.1, -0.4, 0.2, 0.05, -0.01, 0.002])
        deltas = ComplexF64[-1e-3, 0.0, 2e-3]
        v64 = MeijerG.ConfluentJets.exp_divided_difference(L, deltas)
        vdd = MeijerG.ConfluentJets.exp_divided_difference(
            MeijerG.ConfluentJets.Jet(CDD.(L.coeffs)), CDD.(deltas))
        @test abs(v64 - ComplexF64(vdd)) < 1e-13 * abs(v64)
        @test DD(0.1) + 1 isa DD
        @test round(Int, DD(2.6)) == 3
        @test ceil(Int, DD(2.1)) == 3
        @test convert(DD, 1 // 3) isa DD
        @test Float64(DD(pi)) == Float64(pi)
        @test iszero(zero(CDD)) && isfinite(one(CDD))
    end

    @testset "finite-separation lane off the convergent regime" begin
        # lower: b2 = b1 + 1 + 1e-4 at z = 20 (series 1)
        a, b, m, n, z = (0.3,), (0.2, 1.2001), 2, 1, 20.0
        ref = big"0.116703223857659018894238408036303642577908279"
        r = MeijerG._try_meijerg_slater_logratio_jet_near(a, b, m, n, z; return_diagnostics=true)
        @test r.status === :ok && r.selected
        @test r.diagnostics.lane === :extended
        @test abs(r.value - ref) < 1e-14 * abs(ref)
        @test abs(meijerg(a, b, m, n, z) - ref) < 1e-14 * abs(ref)
        rs = MeijerG._try_meijerg_slater_logratio_jet_near(a, b, m, n, z; return_diagnostics=true,
                                                           working_precision=:input)
        @test rs.status === :summation_cancellation
        @test rs.diagnostics.summation_kappa > MeijerG._LOGRATIO_JET_MAX_SUMMATION_KAPPA
        # A declined series goes to contour, not to the separated expansion, which is
        # wrong here.
        @test abs(meijerg_contour(a, b, m, n, z) - ref) < 1e-12 * abs(ref)
        @test !(abs(meijerg_slater(a, b, m, n, z) - ref) < 1e-6 * abs(ref))
        jet = MeijerG._try_meijerg_slater_logratio_jet(a, b, m, n, z; working_precision=:input)
        @test jet.status in MeijerG._LOGRATIO_DECLINED_STATUSES
        # Non-canonical order is canonicalized. Float32 input uses the Float64 path.
        @test meijerg(a, (b[2], b[1]), m, n, z) ≈ meijerg(a, b, m, n, z) rtol=1e-14
        a32 = Float32.(a); b32 = Float32.(b)
        @test meijerg(a32, b32, m, n, Float32(z)) ≈
              meijerg(Float64.(a32), Float64.(b32), m, n, z) rtol=1e-6
        # upper: a2 = a1 - 1 + 1e-4 at z = 0.1 (series 2)
        a, b, m, n, z = (1.2, 0.2001), (0.3,), 1, 2, 0.1
        ref = big"4.49118212466876414146221643152695134120445603"
        r = MeijerG._try_meijerg_slater_logratio_jet_near(a, b, m, n, z; return_diagnostics=true)
        @test r.status === :ok && r.selected
        @test r.diagnostics.lane === :extended
        @test abs(meijerg(a, b, m, n, z) - ref) < 1e-14 * abs(ref)
        @test abs(meijerg((a[2], a[1]), b, m, n, z) - ref) < 1e-14 * abs(ref)
    end

    @testset "separated Slater certificate and lane" begin
        a, b, m, n = (0.3,), (0.2, 0.7, 1.9), 2, 1
        ref50 = big"0.0287400054125126953093808139095445920821706416"
        s = MeijerG._try_meijerg_slater_certified(a, b, m, n, 50.0)
        @test s.status === :ok && s.lane === :extended
        @test s.kappa > MeijerG._LOGRATIO_JET_MAX_SUMMATION_KAPPA
        @test abs(meijerg(a, b, m, n, 50.0) - ref50) < 1e-14 * abs(ref50)
        ss = MeijerG._try_meijerg_slater_certified(a, b, m, n, 50.0; working_precision=:input)
        @test ss.status === :summation_cancellation
        @test !(abs(meijerg_slater(a, b, m, n, 50.0) - ref50) < 1e-13 * abs(ref50))
        @test abs(meijerg_contour(a, b, m, n, 50.0) - ref50) < 1e-12 * abs(ref50)
        # Exact confluence off the convergent regime: the separated expansion is singular
        # in any precision, so the lane must report it (not throw) and `meijerg` must
        # continue to contour integration (reference: mpmath at 40 digits).
        a5, b5 = (0.3,), (0.2, 0.2, 0.2, 0.2, 0.2)
        s5 = MeijerG._try_meijerg_slater_certified(a5, b5, 5, 1, 30.0)
        @test s5.status === :nonfinite && s5.lane === :extended
        @test abs(meijerg(a5, b5, 5, 1, 30.0) - 0.12683318981312108) < 1e-12 * 0.12683318981312108
        # In the convergent regime the native Float64 expansion is certified and used.
        s = MeijerG._try_meijerg_slater_certified(a, b, m, n, 0.05)
        @test s.status === :ok && s.lane === :native
        @test s.kappa <= MeijerG._LOGRATIO_JET_MAX_SUMMATION_KAPPA
        @test meijerg(a, b, m, n, 0.05) == meijerg_slater(a, b, m, n, 0.05)
        # Beyond the lane's own limit, contour integration.
        a3, b3 = (0.3, 0.9), (0.2, 0.7, 1.9)
        ref100 = big"0.0591957012776794760366521348356606702417122646"
        s = MeijerG._try_meijerg_slater_certified(a3, b3, 3, 1, 100.0)
        @test s.status === :summation_cancellation
        @test s.kappa > MeijerG._LOGRATIO_JET_EXTENDED_MAX_KAPPA
        @test abs(meijerg(a3, b3, 3, 1, 100.0) - ref100) < 1e-12 * abs(ref100)
        # Upper branch at small |z|.
        au, bu = (0.2, 0.7, 1.9), (0.3,)
        refu = big"0.322861839840010358426246190560869180835246929"
        s = MeijerG._try_meijerg_slater_certified(au, bu, 1, 2, 0.005)
        @test s.status === :ok && s.lane === :extended
        @test abs(meijerg(au, bu, 1, 2, 0.005) - refu) < 1e-14 * abs(refu)
        # Float32 input runs on the Float64 path, BigFloat input by precision escalation.
        @test meijerg(Float32.(a), Float32.(b), m, n, 50.0f0) isa Float32
        @test abs(meijerg(Float32.(a), Float32.(b), m, n, 50.0f0) - ref50) < 1e-6 * abs(ref50)
        sb = MeijerG._try_meijerg_slater_certified(big.(a), big.(b), m, n, big(50.0))
        @test sb.status === :ok && sb.lane === :escalated
        @test abs(meijerg(big.(a), big.(b), m, n, big(50.0)) - ref50) < 1e-14 * abs(ref50)
        # The certificate itself: inter-family cancellation measured at z = 20.
        v, κ = MeijerG._meijerg_slater_certified(a, b, m, n, 20.0)
        @test 1e6 < κ < 1e7
    end

    @testset "higher-order lane and complete cluster enumeration" begin
        a, b, m, n = (0.3,), (0.2, 1.2, 2.2), 3, 1
        ref20 = big"0.191002138678678707297650321498557947047919041"
        r = MeijerG._try_meijerg_slater_higher_logratio_jet(a, b, m, n, 20.0; return_diagnostics=true)
        @test r.status === :ok && r.selected
        @test r.diagnostics.lane === :extended
        @test abs(r.value - ref20) < 1e-14 * abs(ref20)
        @test abs(meijerg(a, b, m, n, 20.0) - ref20) < 1e-14 * abs(ref20)
        # Strict Float64: the series now runs over complete clusters only (no internal
        # error from a missing cluster member) and is declined on the certificate.
        rs = MeijerG._try_meijerg_slater_higher_logratio_jet(a, b, m, n, 20.0;
                                                             return_diagnostics=true,
                                                             working_precision=:input)
        @test rs.status in (:summation_cancellation, :terms_exceeded)
        @test rs.status in MeijerG._LOGRATIO_DECLINED_STATUSES
        @test abs(meijerg_contour(a, b, m, n, 20.0) - ref20) < 1e-12 * abs(ref20)
        # Beyond the lane's own limit (κ_sum ≈ 1e17 at z = 200) the jet declines and
        # meijerg returns the contour value, not the separated expansion.
        rl = MeijerG._try_meijerg_slater_higher_logratio_jet(a, b, m, n, 200.0;
                                                             return_diagnostics=true)
        @test rl.status === :summation_cancellation && !rl.selected
        @test rl.diagnostics.summation_kappa > MeijerG._LOGRATIO_JET_EXTENDED_MAX_KAPPA
        ref200 = setprecision(BigFloat, 256) do
            Float64(meijerg(map(big, a), map(big, b), m, n, big(200.0)))
        end
        @test meijerg(a, b, m, n, 200.0) == meijerg_contour(a, b, m, n, 200.0)
        @test abs(meijerg(a, b, m, n, 200.0) - ref200) < 1e-12 * abs(ref200)
        d = MeijerG._higher_logratio_branch_sum(a, b, m, n, 20.0, :lower;
                                                maxell=MeijerG._higher_logratio_maxell(20.0))
        @test d.status in (:ok, :terms_exceeded)
        # In-regime behaviour is unchanged: Float64 series, no lane.
        ri = MeijerG._try_meijerg_slater_higher_logratio_jet(a, b, m, n, 0.5; return_diagnostics=true)
        @test ri.status === :ok && ri.selected
        @test !haskey(ri.diagnostics, :lane)
        @test MeijerG._higher_logratio_maxell_extended(20.0, 2, :lower) >
              MeijerG._higher_logratio_maxell_extended(0.5, 2, :lower)
    end
end
