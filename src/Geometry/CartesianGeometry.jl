"""
  StructureDomain1D

Defines a rectangular domain for 1D structural problems on a conforming 2D domain.
The domain is defined by its length `L` and an optional starting point `x₀` in 2D space.

Variables:
- L::Float64 = 1.0: Length of the 1D structure domain
- x₀::Vector{Float64} = [0.0, 1.0]: Starting point of the structure domain in 2D space
"""
@with_kw struct StructureDomain1D
    L::Float64 = 1.0
    x₀::Vector{Float64} = [0.0, 1.0]
end

"""
  DampingZone1D

Defines a 1D damping zone for absorbing outgoing waves, with length `L` and a damping coefficient `σ`.
The damping is applied in the region defined by `x₀` to `x₀ + L` along the x-axis.

Variables:
- L::Float64 = 0.5: Length of the damping zone
- σ::Float64 = 1.0: Damping coefficient
- x₀::Vector{Float64} = [0.0, 1.0]: Starting point of the damping zone in 2D space
"""
@with_kw struct DampingZone1D
    L::Float64 = 0.5
    σ::Float64 = 1.0
    x₀::Vector{Float64} = [0.0, 1.0]
end

"""
  TankDomain2D 

Defines a rectangular domain for 2D problems, with dimensions `L` and `H` [0,L]x[0,H], discretized into `nx` by `ny` elements. 
The `map` function allows for coordinate transformations if needed.

Variables:
- L::Float64 = 4.0: Length of the tank domain
- H::Float64 = 1.0: Height of the tank domain
- nx::Int = 16: Number of elements in the x-direction
- ny::Int = 2: Number of elements in the y-direction
- map::Function = (x, y) -> (x, y): Identity mapping function for coordinates
- structure_domains::Vector{StructureDomain1D} = Vector{StructureDomain1D}(): List of structure domains within the tank
- damping_zones::Vector{DampingZone1D} = Vector{DampingZone1D}(): List of damping zones within the
"""
@with_kw struct TankDomain2D
    L::Float64 = 4.0
    H::Float64 = 1.0
    nx::Int = 16
    ny::Int = 2
    map::Function = x->x
    structure_domains::Vector{StructureDomain1D} = Vector{StructureDomain1D}()
    damping_zones::Vector{DampingZone1D} = Vector{DampingZone1D}()
end

"""
  build_model(domain::TankDomain2D)

Build a `CartesianDiscreteModel` from the specifications in `domain`.
The model is constructed on the rectangular domain defined by `L` and `H`, 
discretized into `nx` by `ny` elements.
"""
function build_model(domain::TankDomain2D)
    x_min, x_max = 0.0, domain.L
    y_min, y_max = 0.0, domain.H
    model = CartesianDiscreteModel((x_min, x_max, y_min, y_max), (domain.nx, domain.ny),map=domain.map)
    return model
end


# ─────────────────────────────────────────────────────────────
# Surface masking
# ─────────────────────────────────────────────────────────────

"""
    _centroid(xs) -> VectorValue

Compute the centroid of a cell given its node coordinates.
"""
function _centroid(xs)
    n = length(xs)
    (1 / n) * sum(xs)
end

"""
    surface_mask(zone::StructureDomain1D) -> Function

Return a closure `(xs) -> Bool` that tests whether the centroid of
a surface cell lies within the structure domain.

The zone spans `[x₀[1], x₀[1]+L]` at `y = x₀[2]`.
"""
function surface_mask(zone::StructureDomain1D)
    x_lo = zone.x₀[1]
    x_hi = zone.x₀[1] + zone.L
    y_ref = zone.x₀[2]
    return function (xs)
        c = _centroid(xs)
        (x_lo <= c[1] <= x_hi) && (c[2] ≈ y_ref)
    end
end

"""
    surface_mask(zone::DampingZone1D) -> Function

Return a closure `(xs) -> Bool` that tests whether the centroid of
a surface cell lies within the damping zone.

The zone spans `[x₀[1], x₀[1]+L]` at `y = x₀[2]`.
"""
function surface_mask(zone::DampingZone1D)
    x_lo = zone.x₀[1]
    x_hi = zone.x₀[1] + zone.L
    y_ref = zone.x₀[2]
    return function (xs)
        c = _centroid(xs)
        (x_lo <= c[1] <= x_hi) && (c[2] ≈ y_ref)
    end
end

"""
    surface_masks(domain::TankDomain2D)

Build mask closures for every structure domain and damping zone in `domain`.

Returns `(structure_masks, damping_masks)` where each is a `Vector{Function}`,
ordered to match `domain.structure_domains` and `domain.damping_zones`.
"""
function surface_masks(domain::TankDomain2D)
    smasks = [surface_mask(s) for s in domain.structure_domains]
    dmasks = [surface_mask(d) for d in domain.damping_zones]
    return smasks, dmasks
end

