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

end
