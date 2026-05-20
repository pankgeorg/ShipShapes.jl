module ShipShapes

using WaterLily
using StaticArrays

export Wigley, wigley_volume

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
    half_beam = (B/2) * (1 - (2xc/L)^2) * (1 - (zc/T)^2)
    in_box = (-L/2 ≤ x ≤ L/2) & (-T ≤ z ≤ 0)
    if in_box
        return abs(y) - half_beam
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

end # module
