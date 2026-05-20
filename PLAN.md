# PLAN 5 — ShipShapes.jl

**Repo:** `pankgeorg/ShipShapes.jl` (private).
**Depends on:** only WaterLily core — needs nothing from PLAN 1.
This package can ship *first* and is the simplest target.

## Scope

Signed distance functions for benchmark ship hulls, exposed as
`AbstractBody`s. Four hulls in priority order:

1. **Wigley parabolic hull** — closed-form analytic SDF. Toy but useful;
   has decades of experimental data for resistance and wave-making.
2. **KCS** (KRISO container ship) — most-used modern ship benchmark.
   Hull offsets tabulated; SDF from a mesh.
3. **DTC** (Duisburg test case) — what OpenFOAM ships in
   `tutorials/incompressibleVoF/DTCHull*`. Matches our integration tests.
4. **DTMB 5415** — US Navy combatant; canonical for resistance and
   manoeuvring studies.

Also bundles utilities:
- Rigid-body kinematics helpers (heave, pitch, sway, yaw at prescribed
  motion or 6DoF).
- Hull-fixed reference frame coordinate transforms.
- A "double-body" mode that mirrors the hull about the waterline (for
  pre-Phase-2 single-phase tests without VoF).

## Non-goals

- **No hull-form design / parametric morphing.** That's a research field.
- **No** appendages (rudders, bilge keels) packaged with the hulls.
  Provide as separate `AbstractBody`s users can compose with `+`/`-`.
- **No** structural deformation. Hulls are rigid.
- **No** stability/equilibrium solver. Floating equilibrium is
  ShipFlow.jl's problem.

## API

```julia
using WaterLily, ShipShapes
using StaticArrays

# Analytic Wigley — no I/O needed
hull = ShipShapes.Wigley(L = 2.5, B = 0.25, T = 0.156)

# Tabulated hull — loads packaged geometry data
hull = ShipShapes.KCS()                     # default scale (Lpp = 230 m model)
hull = ShipShapes.DTC(scale = 1/59.4)       # model scale used in DTC paper
hull = ShipShapes.DTMB5415()

# Compose with rigid motion
moving = ShipShapes.moving_hull(hull;
    velocity = (t)-> SVector(U, 0, 0),
    rotation = (t)-> RotZ(0.0),
)

# Double-body mode for pre-VoF tests
sub = ShipShapes.double_body(hull, waterline_z = 0.0)

Simulation(dims, uBC, L; body = moving, ν, ...)
```

All hulls implement the `AbstractBody` interface: `measure(body, x, t) = (d, n, V)`.

## Algorithms

- Wigley analytic: `y(x,z) = (B/2)(1 - (2x/L)²)(1 - (z/T)²)`. SDF in closed
  form via projection onto the surface.
- Tabulated hulls: hull offset tables → triangulated mesh → SDF via
  point-to-triangle distance with BVH acceleration. Use
  `Meshing.jl` + a hand-rolled BVH (or `RegionTrees.jl` /
  `Distances.jl` if a clean dep exists).
- Geometry data: bundled as compressed offset tables (`hulls/*.csv` or
  packed JLD2). For KCS and DTC, sourced from the SIMMAN workshop and
  the el Moctar et al. 2012 paper respectively — both publicly
  redistributable, document the provenance.

## Validation

### Layer 1 — geometric (fast, per-PR CI)

- **Wigley SDF correctness.** Sample 10⁴ random points; check
  `|d - d_analytic| < 1e-6 L`. Check `∇d = n` via AD (`ForwardDiff`)
  and confirm `|n| ≈ 1`.
- **Volume / displaced mass.** Numerical integration of indicator
  function over a fine grid should match the published hull
  displacement within ±0.5% for each hull at the design draught.
  Published values:
  - Wigley: analytic
  - KCS: 52030 m³ (full scale)
  - DTC: 11377 m³
  - DTMB 5415: 8424 m³
- **Wetted surface area.** Within ±1% of published values.
- **Sign and continuity of SDF.** No interior point has positive
  distance; no exterior point negative. Crossings clean (no
  discontinuities > 2 Δx in a sweep).

### Layer 2 — comparison with OpenFOAM-meshed hull (nightly CI)

For DTC specifically (since we have the OpenFOAM tutorial in-tree):

1. Run `snappyHexMesh` in the `DTCHull` tutorial container to generate
   the body-fitted mesh.
2. Extract the body surface as STL.
3. Sample our SDF on the OpenFOAM mesh nodes and compare wetted-surface
   geometry: location of stem, stern, midship section, design waterline.
   All within ±0.5% of the OpenFOAM mesh extents.

This isn't a flow test — purely geometric — but it catches data-loading
mistakes that L1 alone would miss.

### Layer 3 — release-blocking flow validation

(via ShipFlow.jl, using bare WaterLily — no VoF, no propeller.)

Drag coefficient `C_T` at zero Froude number (double-body) on a
KCS hull at model scale, compared to the SIMMAN workshop towing-tank
data. Pass: within ±20% on first cut — bare WaterLily without
turbulence isn't expected to match closely; this is mostly a
"does the geometry produce a sane flow" check.

## Performance budget

- Wigley SDF: < 50 ns per `measure` call (closed form).
- Tabulated hulls: < 1 µs per `measure` call with BVH.
- `measure!(flow, body)` over a 256³ grid: < 1 second on CPU,
  < 100 ms on GPU.

If tabulated-hull SDF is the per-step bottleneck, precompute a
coarse-grid SDF cache + trilinear interpolation in a band around the
body.

## Harness

- Layer 1 tests in `test/runtests.jl`. Pure geometry; no fluid solver
  invoked. Runs in seconds.
- Layer 2 tests in `test/openfoam/`, nightly. Needs the OpenFOAM
  Docker image to mesh DTC; output is cached so we don't re-run unless
  the SDF code changes.
- Layer 3 is hosted in ShipFlow.jl as part of the integration suite.

## Risks & open questions

- **Hull data licensing.** SIMMAN benchmark hulls are public-domain for
  research use; verify before committing the offset files. If a hull's
  redistribution is restricted, ship a fetcher script instead of the
  data itself.
- **SDF accuracy in concave regions** (transom, bulbous bow). BVH
  point-to-triangle gives correct *signed* distance only with a
  consistent inside/outside test — use generalized winding number
  (Jacobson et al. 2013) rather than naïve normal voting.
- **Coordinate frame conventions.** OpenFOAM, WaterLily, and ship
  hydrodynamics all use different axis conventions. Document and
  enforce a single one (e.g., x forward, y starboard, z up) at the
  package boundary; convert at I/O.

## Milestones

| # | Goal                                            | Done when                              |
|---|-------------------------------------------------|----------------------------------------|
| 1 | Wigley analytic SDF + L1 tests                  | Layer 1 passes                         |
| 2 | KCS tabulated SDF + L1 + L2                     | volume ±0.5%, OpenFOAM extents match   |
| 3 | DTC SDF (matches our OpenFOAM tutorial)         | feeds VoF.jl Phase-3 integration test  |
| 4 | DTMB 5415                                       | L1 passes                              |
| 5 | Rigid-motion helpers + double-body mode         | KCS drag ±20% in ShipFlow.jl           |
