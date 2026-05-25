module ShipShapes

using WaterLily
using StaticArrays

export Wigley, wigley_volume, Containership, containership_volume,
       TabulatedHull, sample_sdf, tabulated_sdf

"""
    wigley_sdf(p, L, B, T)

Approximate signed distance from point `p::SVector{3}` (in hull-fixed coords,
midship at origin, waterline at `z=0`, keel at `z=-T`) to the Wigley
parabolic hull defined by

```
y_surface(x,z) = ±(B/2) * (1 - (2x/L)²) * (1 - (z/T)²),  x ∈ [-L/2, L/2], z ∈ [-T, 0]
```

The returned value is a pseudo-SDF: signed (negative inside), zero on the
surface, and consistent with the geometric normal direction. Magnitude
is corrected by `WaterLily.AutoBody`'s gradient normalization, so this
is suitable for use with `AutoBody(wigley_sdf_closure)`.
"""
@inline function wigley_sdf(p, L, B, T)
    x, y, z = p[1], p[2], p[3]
    xc = clamp(x, -L/2, L/2)
    zc = clamp(z, -T, 0)
    fx = 1 - (2xc/L)^2          # 1 at midship, 0 at bow/stern
    fz = 1 - (zc/T)^2           # 1 at waterline, 0 at keel
    half_beam = (B/2) * fx * fz
    in_box = (-L/2 ≤ x ≤ L/2) & (-T ≤ z ≤ 0)
    if in_box
        # Eikonal-normalised SDF: divide the y-axis distance to the
        # surface by |∇half_beam_extended|, where the gradient lives in
        # all three directions because the surface curves in x and z.
        # ∂half_beam/∂x = -2B·x/L² · fz,  ∂half_beam/∂z = -B·z/T² · fx.
        dhdx = -2 * B * xc / (L^2) * fz
        dhdz = -B * zc / (T^2) * fx
        gnorm = sqrt(1 + dhdx^2 + dhdz^2)
        return (abs(y) - half_beam) / gnorm
    else
        y_surf = clamp(y, -half_beam, half_beam)
        dx = x - xc
        dz = z - zc
        dy = y - y_surf
        return sqrt(dx^2 + dy^2 + dz^2)
    end
end

"""
    Wigley(; L, B, T, map=(x,t)->x)

Wigley parabolic hull as a `WaterLily.AutoBody`. Coordinate frame:
midship at the origin, x along the length, y across the beam, z vertical
(positive up, waterline at z=0, keel at z=-T).

Pass `map` to translate / rotate / animate the hull in the world frame.

# Example
```julia
hull = Wigley(L=2.5, B=0.25, T=0.156)
sim = Simulation((256,64,32), (1.0,0,0), 2.5; body=hull, ν=1e-5)
```
"""
function Wigley(; L::Real, B::Real, T::Real, map=(x,t)->x)
    L_, B_, T_ = float(L), float(B), float(T)
    sdf(x, t) = wigley_sdf(x, L_, B_, T_)
    AutoBody(sdf, map)
end

"""
    wigley_volume(L, B, T)

Analytic displaced-water volume of the Wigley hull: `4 L B T / 9`.
"""
wigley_volume(L, B, T) = 4 * L * B * T / 9

"""
    containership_sdf(p, L, B, T, par_frac=0.5)

Approximate signed distance to a containership-style hull: parallel
mid-body of length `par_frac × L` at full beam B and full draught T,
linearly tapered to a point at bow and stern. Vertical sides (Cb
much higher than Wigley's 0.44 — typically 0.65-0.70).

Coordinate frame: midship at origin, x = length, y = beam, z =
vertical (waterline at z=0, keel at z=-T). Used as a stand-in for
the Duisburg Test Case (DTC) hull until proper offset data is
available.
"""
@inline function containership_sdf(p, L, B, T, par_frac=0.5)
    x, y, z = p[1], p[2], p[3]
    xc = clamp(x, -L/2, L/2)
    zc = clamp(z, -T, 0)
    # Plan-view half-beam: full B/2 over the parallel mid-body
    # |2x/L| ≤ par_frac, linearly tapered outside.
    s = abs(2*xc / L)
    if s ≤ par_frac
        half_beam_x = B/2
    else
        # taper from full beam at s=par_frac to 0 at s=1
        half_beam_x = (B/2) * (1 - s) / (1 - par_frac)
    end
    in_box = (-L/2 ≤ x ≤ L/2) & (-T ≤ z ≤ 0)
    if in_box
        # SDF inside: distance to the rectangular cross-section
        # bounded by ±half_beam_x in y and 0..-T in z. Vertical
        # sides mean the section is a rectangle in y-z.
        return abs(y) - half_beam_x
    else
        y_surf = clamp(y, -half_beam_x, half_beam_x)
        dx = x - xc
        dz = z - zc
        dy = y - y_surf
        return sqrt(dx^2 + dy^2 + dz^2)
    end
end

"""
    Containership(; L, B, T, par_frac=0.5, map=(x,t)->x)

Containership-style hull (parallel mid-body + linear-tapered ends,
vertical sides) as a `WaterLily.AutoBody`. Approximation for the
Duisburg Test Case (DTC) family. Cb ≈ `(1 + par_frac) / 2`, so the
default `par_frac=0.5` gives Cb ≈ 0.75.

```julia
hull = Containership(L=355, B=51, T=14.5)
```
"""
function Containership(; L::Real, B::Real, T::Real, par_frac::Real=0.5,
                         map=(x,t)->x)
    L_, B_, T_, p_ = float(L), float(B), float(T), float(par_frac)
    sdf(x, t) = containership_sdf(x, L_, B_, T_, p_)
    AutoBody(sdf, map)
end

"""
    containership_volume(L, B, T, par_frac=0.5)

Analytic displaced-water volume of the containership stand-in:
parallel midbody of length `par_frac·L` × `B` × `T` plus two
triangular-prism end-caps of length `(1−par_frac)·L/2` × half-beam
average × `T`.
"""
function containership_volume(L, B, T, par_frac=0.5)
    L_par = par_frac * L
    L_end = (1 - par_frac) * L / 2     # each end
    # midbody = L_par × B × T; each end is a triangular prism with
    # base L_end and tapering half-beam from B/2 to 0 → volume
    # per end = (1/2) · L_end · (B/2) · 2 · T = L_end · B/2 · T.
    L_par * B * T + 2 * (L_end * B / 2 * T)
end

# ----------------------------------------------------------------------------
# Tabulated SDF — for hulls without a closed-form distance (DTC, KCS, …)
# ----------------------------------------------------------------------------

"""
    TabulatedHull{T,N}

A signed-distance function sampled on a uniform 3D grid. Use this when
you have a hull whose surface is defined by a triangle mesh or a set of
sample points, but where evaluating the SDF analytically is expensive.

# Fields
- `grid` — `Array{T,3}` of distance values, with `grid[i,j,k]` =
  SDF at world position `origin + spacing .* (i-1, j-1, k-1)`.
- `origin` — `SVector{3,T}`, world-frame coordinates of `grid[1,1,1]`.
- `spacing` — `SVector{3,T}`, grid step in each axis.

Trilinear interpolation is used for queries between sample points.
Outside the sampled box, the SDF returns a positive value (clamped
to "far outside") — *do not* place a body so that its surface goes
through the SDF box boundary.
"""
struct TabulatedHull{T, A<:AbstractArray{T,3}}
    grid    :: A
    origin  :: SVector{3,T}
    spacing :: SVector{3,T}
end

function TabulatedHull(grid::AbstractArray{T,3},
                       origin, spacing) where T
    o = SVector{3,T}(origin)
    s = SVector{3,T}(spacing)
    TabulatedHull{T, typeof(grid)}(grid, o, s)
end

# Trilinear interpolation, returning a value with the same type as the
# input coordinate (so AutoBody / ForwardDiff Duals stay differentiable).
@inline function (h::TabulatedHull{T})(p, _t = 0) where T
    nx, ny, nz = size(h.grid)
    Tin = promote_type(eltype(p), T)
    # World → grid index (continuous)
    gx = (p[1] - h.origin[1]) / h.spacing[1]
    gy = (p[2] - h.origin[2]) / h.spacing[2]
    gz = (p[3] - h.origin[3]) / h.spacing[3]
    # Outside the tabulation box: return a continuous distance-to-box plus
    # a one-cell margin so the SDF stays C0 across the box face. Returning a
    # constant (the old behaviour) made BDIM's ∇sdf finite-difference see a
    # spurious surface at the box boundary.
    if gx < 0 || gy < 0 || gz < 0 || gx > nx-1 || gy > ny-1 || gz > nz-1
        dx = max(zero(Tin), Tin(-gx), Tin(gx - (nx - 1))) * Tin(h.spacing[1])
        dy = max(zero(Tin), Tin(-gy), Tin(gy - (ny - 1))) * Tin(h.spacing[2])
        dz = max(zero(Tin), Tin(-gz), Tin(gz - (nz - 1))) * Tin(h.spacing[3])
        d_box = sqrt(dx*dx + dy*dy + dz*dz)
        # Anchor at the (positive) value on the nearest face — the
        # tabulation is supposed to have valid SDF values out to its
        # boundary, so reading the corner approximates the face SDF.
        i = clamp(round(Int, gx) + 1, 1, nx)
        j = clamp(round(Int, gy) + 1, 1, ny)
        k = clamp(round(Int, gz) + 1, 1, nz)
        @inbounds face_val = Tin(h.grid[i, j, k])
        return face_val + d_box
    end
    i = clamp(floor(Int, gx), 0, nx-2)
    j = clamp(floor(Int, gy), 0, ny-2)
    k = clamp(floor(Int, gz), 0, nz-2)
    fx = gx - i; fy = gy - j; fz = gz - k
    i += 1; j += 1; k += 1                    # 1-indexed
    g = h.grid
    @inbounds c000 = Tin(g[i,   j,   k  ])
    @inbounds c100 = Tin(g[i+1, j,   k  ])
    @inbounds c010 = Tin(g[i,   j+1, k  ])
    @inbounds c110 = Tin(g[i+1, j+1, k  ])
    @inbounds c001 = Tin(g[i,   j,   k+1])
    @inbounds c101 = Tin(g[i+1, j,   k+1])
    @inbounds c011 = Tin(g[i,   j+1, k+1])
    @inbounds c111 = Tin(g[i+1, j+1, k+1])
    c00 = c000 * (1 - fx) + c100 * fx
    c01 = c001 * (1 - fx) + c101 * fx
    c10 = c010 * (1 - fx) + c110 * fx
    c11 = c011 * (1 - fx) + c111 * fx
    c0  = c00  * (1 - fy) + c10  * fy
    c1  = c01  * (1 - fy) + c11  * fy
    return c0 * (1 - fz) + c1 * fz
end

"""
    sample_sdf(sdf::Function, origin, spacing, dims; T=Float32)

Tabulate an analytic SDF on a uniform grid. Returns a `TabulatedHull`.
`origin`, `spacing` are `NTuple{3}`; `dims` is `NTuple{3,Int}` giving the
sample count along each axis.

```julia
hull_analytic = Wigley(L=L, B=B, T=T)
hull_table = sample_sdf((x,t) -> ShipShapes.wigley_sdf(x, L, B, T),
                        (-L/2 - L/10, -B - 0.1, -T - 0.1),  # origin
                        (1.2L/63, 2.2B/31, 1.2T/31),         # spacing
                        (64, 32, 32))
```
"""
function sample_sdf(sdf, origin, spacing, dims::NTuple{3,Int}; T::Type=Float32)
    nx, ny, nz = dims
    o = SVector{3,T}(origin)
    s = SVector{3,T}(spacing)
    grid = Array{T,3}(undef, nx, ny, nz)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        x = SVector{3,T}(o[1] + s[1] * (i-1),
                         o[2] + s[2] * (j-1),
                         o[3] + s[3] * (k-1))
        grid[i, j, k] = T(sdf(x, T(0)))
    end
    return TabulatedHull(grid, o, s)
end

"""
    tabulated_sdf(table::TabulatedHull; map=(x,t)->x)

Convert a `TabulatedHull` into a `WaterLily.AutoBody` ready for use in
a Simulation. `map` is the standard WaterLily coordinate-mapping
function.
"""
function tabulated_sdf(table::TabulatedHull; map=(x,t)->x)
    WaterLily.AutoBody((x, t) -> table(x, t), map)
end

end # module
