using Test
using MeijerG

@testset "Private active confluence cluster detection" begin
    @testset "lower exact and near shifted family clusters" begin
        b0 = -0.31 + 0.11im
        a = (0.42 + 0.23im, 1.16 - 0.18im, 1.77 + 0.19im)
        b = (b0, b0 + 1, b0 + 2, 1.37 + 0.22im, 1.91 - 0.17im)
        s = MeijerG._active_confluence_summary(a, b, 3, 1, 3e-4 + 0im)
        @test s.status === :ok
        @test s.branch === :lower
        @test s.max_multiplicity == 3
        @test length(s.clusters) == 1
        @test only(s.clusters).kind === :exact
        @test only(s.clusters).families == (1, 2, 3)
        @test only(s.clusters).ells == (2, 1, 0)

        bnear = (b0, b0 + 1 + 2e-4, b0 + 2 - 3e-4,
                 1.37 + 0.22im, 1.91 - 0.17im)
        s = MeijerG._active_confluence_summary(a, bnear, 3, 1, 3e-4 + 0im)
        @test s.status === :ok
        @test only(s.clusters).kind === :near
        @test only(s.clusters).multiplicity == 3
        @test only(s.clusters).diameter ≈ 5e-4 atol=1e-15
    end

    @testset "upper exact and near shifted family clusters" begin
        a0 = -0.22 + 0.18im
        a = (a0, a0 - 1, a0 - 2, 1.41 - 0.22im, 2.03 + 0.31im)
        b = (0.26 + 0.31im, 1.12 - 0.16im, -0.44 + 0.27im)
        s = MeijerG._active_confluence_summary(a, b, 1, 3, 50.0 + 0im)
        @test s.status === :ok
        @test s.branch === :upper
        @test only(s.clusters).kind === :exact
        @test only(s.clusters).multiplicity == 3
        @test only(s.clusters).families == (1, 2, 3)

        anear = (a0, a0 - 1 + 1e-4, a0 - 2 - 2e-4,
                 1.41 - 0.22im, 2.03 + 0.31im)
        s = MeijerG._active_confluence_summary(anear, b, 1, 3, 50.0 + 0im)
        @test s.status === :ok
        @test only(s.clusters).kind === :near
        @test only(s.clusters).multiplicity == 3
    end

    @testset "arbitrary multiplicity detection is not capped" begin
        b0 = 0.18 - 0.07im
        a = (0.42 + 0.23im, 1.16 - 0.18im, 1.77 + 0.19im)
        b = (ntuple(k -> b0 + (k - 1), 5)..., 3.2 + 0.4im)
        s = MeijerG._active_confluence_summary(a, b, 5, 1, 3e-4 + 0im)
        @test s.status === :ok
        @test s.max_multiplicity == 5
        @test only(s.clusters).multiplicity == 5
        @test only(s.clusters).kind === :exact
    end

    @testset "safe and unsupported geometries" begin
        a = (0.42 + 0.23im, 1.16 - 0.18im)
        b = (-0.31 + 0.11im, 0.48 - 0.23im, 1.37 + 0.22im)
        s = MeijerG._active_confluence_summary(a, b, 2, 1, 3e-4 + 0im)
        @test s.status === :not_near_confluent
        @test isempty(s.clusters)

        s = MeijerG._active_confluence_summary((0.3, 0.8), (0.2, 0.2 + 1e-6),
                                               2, 1, 0.5 + 0im)
        @test s.status === :unsupported_finite_p_eq_q

        b0 = -0.31 + 0.11im
        bden = 1 + b0
        bbad = (b0, b0 + 1e-6, bden)
        s = MeijerG._active_confluence_summary(a, bbad, 2, 1, 3e-4 + 0im)
        @test s.status === :unsupported_local_structure
    end

    @testset "multiple independent active clusters are visible" begin
        a = (0.42 + 0.23im, 1.16 - 0.18im, 1.77 + 0.19im)
        b1 = -0.31 + 0.11im
        b2 = 0.73 - 0.24im
        b = (b1, b1 + 1e-5, b2, b2 + 1 + 2e-5, 1.91 - 0.17im)
        s = MeijerG._active_confluence_summary(a, b, 4, 1, 3e-4 + 0im)
        @test s.status === :unsupported_multiple_clusters
        @test length(s.clusters) == 2
        @test all(c -> c.multiplicity == 2, s.clusters)
    end

    @testset "branch override and pole-location helper" begin
        a0 = 0.2 + 0.1im
        a = (a0, a0 - 1 + 1e-4, 1.5 + 0.2im)
        b = (0.3 - 0.1im, 1.2 + 0.4im)
        s = MeijerG._active_confluence_summary(a, b, 1, 2, 20.0 + 0im; branch=:upper)
        @test s.status === :ok
        @test s.branch === :upper
        @test only(s.clusters).multiplicity == 2

        pole_clusters = MeijerG._detect_active_confluence_clusters(a, b, 1, 2, :upper;
                                                                   maxell=3,
                                                                   near_radius=1e-3)
        @test !isempty(pole_clusters)
        @test any(c -> c.multiplicity == 2 && c.kind === :near, pole_clusters)
    end
end
