using Test
using ShipShapes
using ShipShapes: wigley_sdf, wigley_volume
using StaticArrays

@testset "ShipShapes" begin

    @testset "Wigley analytic SDF" begin
        L, B, T = 2.5, 0.25, 0.156

        # Sign: midship-keel interior point is inside
        @test wigley_sdf(SVector(0.0, 0.0, -T/2), L, B, T) < 0

        # Sign: well outside in y is outside
        @test wigley_sdf(SVector(0.0, B, -T/2), L, B, T) > 0

        # Sign: above the waterline is outside
        @test wigley_sdf(SVector(0.0, 0.0, 0.1), L, B, T) > 0

        # Sign: below the keel is outside
        @test wigley_sdf(SVector(0.0, 0.0, -T - 0.1), L, B, T) > 0

        # Sign: ahead of the bow is outside
        @test wigley_sdf(SVector(L, 0.0, -T/2), L, B, T) > 0

        # Surface: midship max-beam, half-draught — should be exactly zero
        half_beam_mid = (B/2) * 1 * (1 - (-T/2/T)^2)
        @test isapprox(
            wigley_sdf(SVector(0.0, half_beam_mid, -T/2), L, B, T),
            0.0;
            atol = 1e-12,
        )

        # Bow / stern have zero half-beam — surface passes through y=0
        @test isapprox(wigley_sdf(SVector( L/2, 0.0, -T/2), L, B, T), 0.0; atol=1e-12)
        @test isapprox(wigley_sdf(SVector(-L/2, 0.0, -T/2), L, B, T), 0.0; atol=1e-12)
    end

    @testset "Wigley volume via Monte Carlo" begin
        # Sample N points in a box bounding the hull, count interior fraction
        L, B, T = 2.5, 0.25, 0.156
        V_analytic = wigley_volume(L, B, T)
        @test isapprox(V_analytic, 4 * L * B * T / 9; rtol=1e-12)

        # Box volume
        boxL, boxB, boxT = 1.2*L, 1.2*B, 1.2*T
        V_box = boxL * boxB * boxT

        # Deterministic seed via plain LCG to keep this test reproducible.
        N = 200_000
        rng_state = UInt64(2025_05_20)
        @inline function next_unit()
            rng_state = rng_state * 6364136223846793005 + 1442695040888963407
            return (rng_state >> 11) / (1 << 53)
        end
        inside = 0
        for _ in 1:N
            x = (next_unit() - 0.5) * boxL
            y = (next_unit() - 0.5) * boxB
            z = (next_unit() - 1.0) * boxT  # z in [-boxT, 0]
            inside += wigley_sdf(SVector(x, y, z), L, B, T) < 0 ? 1 : 0
        end
        V_mc = inside / N * V_box

        # MC error scales as 1/√N ~ 0.22% on volume fraction at N=2e5;
        # generous 5% tolerance keeps the test stable run-to-run.
        @test isapprox(V_mc, V_analytic; rtol = 0.05)
    end

    @testset "Wigley wraps as AutoBody (interface only)" begin
        # We don't import WaterLily here directly to avoid the test runtime
        # depending on KernelAbstractions startup. Just smoke-check that the
        # constructor returns *something* and is callable.
        hull = Wigley(L=2.5, B=0.25, T=0.156)
        @test hull !== nothing
        # AutoBody.sdf should evaluate at a known interior point
        @test hull.sdf(SVector(0.0, 0.0, -0.05), 0.0) < 0
        @test hull.sdf(SVector(0.0, 0.5, -0.05), 0.0) > 0
    end

    @testset "TabulatedHull: sample analytic Wigley and re-evaluate" begin
        L, B, T = 2.5, 0.25, 0.156
        ana = (x, t) -> wigley_sdf(x, L, B, T)

        # Sample box: 1.4× analytic extent
        ox, oy, oz = -0.7 * L, -0.7 * B, -1.4 * T
        spx = 1.4 * L / 99
        spy = 1.4 * B / 49
        spz = 1.4 * T / 49
        table = sample_sdf(ana, (ox, oy, oz), (spx, spy, spz), (100, 50, 50);
                           T = Float64)

        # At sample-aligned points the value should be bit-identical
        for (ig, jg, kg) in [(20, 25, 25), (50, 10, 15), (80, 40, 35)]
            x = SVector(ox + ig * spx, oy + jg * spy, oz + kg * spz)
            @test isapprox(table(x, 0.0), ana(x, 0.0); atol=1e-12)
        end

        # At interpolation points the value should match analytic within
        # second-order interpolation error.
        for (ig, jg, kg) in [(20.5, 25.5, 25.5), (50.5, 10.5, 15.5)]
            x = SVector(ox + ig * spx, oy + jg * spy, oz + kg * spz)
            d_tab = table(x, 0.0)
            d_ana = ana(x, 0.0)
            @test isapprox(d_tab, d_ana; atol = 0.02)  # ~half a cell tolerance
        end

        # Outside the sample box → clamped to "far outside"
        x_far = SVector(10.0, 0.0, 0.0)
        @test table(x_far, 0.0) > 0
    end

    @testset "TabulatedHull: tabulated_sdf wraps as AutoBody" begin
        L, B, T = 2.5, 0.25, 0.156
        ana = (x, t) -> wigley_sdf(x, L, B, T)
        table = sample_sdf(ana, (-1.5, -0.2, -0.25),
                           (3.0/49, 0.4/19, 0.3/19), (50, 20, 20))
        hull = tabulated_sdf(table)
        @test hull !== nothing
        # Interior point — table interpolates close to analytic
        x_in = SVector(0.0, 0.0, -0.05)
        @test hull.sdf(x_in, 0.0) ≈ ana(x_in, 0.0) atol=0.02
    end

    @testset "Containership hull SDF + volume" begin
        L, B, T = 5.0, 0.7, 0.3
        # Inside the parallel midbody, half-beam = B/2 everywhere
        @test ShipShapes.containership_sdf(SVector(0.0, 0.0, -T/2), L, B, T) < 0
        @test ShipShapes.containership_sdf(SVector(0.0, 0.4, -T/2), L, B, T) > 0
        # At bow (s=1) half-beam → 0
        @test ShipShapes.containership_sdf(SVector(L/2, 0.0, -T/2), L, B, T) ≈ 0 atol=1e-9
        @test ShipShapes.containership_sdf(SVector(-L/2, 0.0, -T/2), L, B, T) ≈ 0 atol=1e-9
        # Cb ≈ (1 + par_frac)/2 = 0.75 for par_frac=0.5
        V = containership_volume(L, B, T, 0.5)
        @test V / (L * B * T) ≈ 0.75 atol=1e-9
        # Higher par_frac → higher Cb
        @test containership_volume(L, B, T, 0.8) > containership_volume(L, B, T, 0.5)
    end

    @testset "Containership body runs in a WaterLily Simulation" begin
        using WaterLily
        L_c = 30f0; B_c = 7f0; T_c = 4f0
        hull_xc = 20f0; hull_yc = 16f0; hull_zc = 16f0
        hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
        hull = Containership(; L=L_c, B=B_c, T=T_c, par_frac=0.5, map=hull_map)
        sim = WaterLily.Simulation((48, 32, 32), (1f0, 0f0, 0f0), L_c;
            T=Float32, body=hull, Δt=0.25f0, ϵ=1, U=1f0)
        for _ in 1:5
            WaterLily.mom_step!(sim.flow, sim.pois)
        end
        @test isfinite(maximum(abs, sim.flow.u))
        @test maximum(abs, sim.flow.u) < 5f0
    end

    @testset "Wigley body runs in a WaterLily Simulation" begin
        # End-to-end: build a Simulation with a Wigley body and confirm
        # a few mom_step!s don't blow up. Smoke-level proof that the
        # SDF is BDIM-compatible.
        using WaterLily
        L_c = 24f0; B_c = 6f0; T_c = 4f0
        hull_xc = 16f0; hull_yc = 16f0; hull_zc = 16f0
        hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
        hull = Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
        sim = WaterLily.Simulation((48, 32, 32), (1f0, 0f0, 0f0), L_c;
            T = Float32, body = hull, Δt = 0.25f0, ϵ = 1, U = 1f0)
        for _ in 1:5
            WaterLily.mom_step!(sim.flow, sim.pois)
        end
        @test isfinite(maximum(abs, sim.flow.u))
        @test maximum(abs, sim.flow.u) < 5f0   # well-bounded
    end

    @testset "wigley_sdf gradient ≈ 1 (first-order Eikonal)" begin
        # The Eikonal-normalised SDF should have |∇φ| ≈ 1 to first
        # order; second-order curvature terms can push it a few %
        # above 1. Check several interior points and require |∇φ| < 1.2.
        L, B, T = 3.0, 0.4, 0.3
        h = 1e-3
        pts = (SVector(0.0, 0.05, -0.05),    # close to midship, near waterline
               SVector(0.0, 0.05, -0.15),    # midship, mid-depth
               SVector(0.5, 0.05, -0.10))    # 1/3 from bow
        for p in pts
            fx = (ShipShapes.wigley_sdf(SVector(p[1]+h, p[2], p[3]), L, B, T) -
                  ShipShapes.wigley_sdf(SVector(p[1]-h, p[2], p[3]), L, B, T)) / (2h)
            fy = (ShipShapes.wigley_sdf(SVector(p[1], p[2]+h, p[3]), L, B, T) -
                  ShipShapes.wigley_sdf(SVector(p[1], p[2]-h, p[3]), L, B, T)) / (2h)
            fz = (ShipShapes.wigley_sdf(SVector(p[1], p[2], p[3]+h), L, B, T) -
                  ShipShapes.wigley_sdf(SVector(p[1], p[2], p[3]-h), L, B, T)) / (2h)
            g = sqrt(fx^2 + fy^2 + fz^2)
            @test 0.85 ≤ g ≤ 1.20
        end
    end

end
