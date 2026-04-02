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
    domain_symbol::Symbol = :Γ_s
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
    domain_symbol::Symbol = :Γ_d
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


# ─────────────────────────────────────────────────────────────
# Triangulations
# ─────────────────────────────────────────────────────────────

"""
    build_triangulations(domain::TankDomain2D, model) -> TankTriangulations

Label the model, extract base triangulations (`Ω`, `Γ`, `Γin`, `Γot`),
then partition the top surface `Γ` using the masks derived from the
structure domains and damping zones in `domain`.

Returns a [`TankTriangulations`](@ref) holding every sub-triangulation
needed for simulation assembly.

Labelling convention (Cartesian 2D, entity ids):
- `"surface"` → entities 3, 4, 6 (top side + top corners)
- `"bottom"`  → entities 1, 2, 5 (bottom side + bottom corners)
- `"inlet"`   → entity 7 (left wall)
- `"outlet"`  → entity 8 (right wall)
- `"water"`   → entity 9 (interior)
"""
function build_triangulations(domain::TankDomain2D, model)
    # — Label model faces ————————————————————————————————
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "surface", [3, 4, 6])
    add_tag_from_tags!(labels, "bottom",  [1, 2, 5])
    add_tag_from_tags!(labels, "inlet",   [7])
    add_tag_from_tags!(labels, "outlet",  [8])
    add_tag_from_tags!(labels, "water",   [9])

    # — Base triangulations ——————————————————————————————
    Ω    = Interior(model)
    Γ    = Boundary(model, tags="surface")
    Γbot = Boundary(model, tags="bottom")
    Γin  = Boundary(model, tags="inlet")
    Γout = Boundary(model, tags="outlet")

    # — Build boolean masks on surface cell coordinates ——
    xΓ = get_cell_coordinates(Γ)
    smasks, dmasks = surface_masks(domain)

    s_bits = [lazy_map(m, xΓ) for m in smasks]  # per-structure
    d_bits = [lazy_map(m, xΓ) for m in dmasks]  # per-damping

    # Build structure triangulations and map to their domain_symbol
    Γ_structures = Vector{Any}(undef, length(s_bits))
    structure_syms = Symbol[]
    structure_trians = Dict{Symbol, Any}()
    seen_structure_syms = Set{Symbol}()
    for (i, b) in enumerate(s_bits)
        tri = Triangulation(Γ, findall(b))
        Γ_structures[i] = tri
        sym = domain.structure_domains[i].domain_symbol
        if sym in seen_structure_syms
            error("Duplicate domain_symbol :$sym found in structure_domains. Each structure must have a unique domain_symbol.")
        end
        push!(seen_structure_syms, sym)
        structure_trians[sym] = tri
        push!(structure_syms, sym)
    end

    # Damping triangulations and map to their domain_symbol
    Γ_dampings = Vector{Any}(undef, length(d_bits))
    damping_syms = Symbol[]
    damping_trians = Dict{Symbol, Any}()
    seen_damping_syms = Set{Symbol}()
    for (i, b) in enumerate(d_bits)
        tri = Triangulation(Γ, findall(b))
        Γ_dampings[i] = tri
        sym = domain.damping_zones[i].domain_symbol
        if sym in seen_damping_syms
            error("Duplicate domain_symbol :$sym found in damping_zones. Each damping zone must have a unique domain_symbol.")
        end
        push!(seen_damping_syms, sym)
        damping_trians[sym] = tri
        push!(damping_syms, sym)
    end

    # — Composed masks ——————————————————————————————————
    n = length(xΓ)
    any_structure = _or_bits(s_bits, n)
    any_damping   = _or_bits(d_bits, n)
    any_zone      = any_structure .| any_damping

    # Free surface = surface cells in no zone at all
    Γfs = Triangulation(Γ, findall(!, any_zone))
    # κ surface = everything except structures (free surface + damping)
    Γκ  = Triangulation(Γ, findall(!, any_structure))
    # η surface = all structures
    Γη  = Triangulation(Γ, findall(any_structure))

    # Compose dictionary for TankTriangulations
    trian_dict = Dict(
        :Ω => Ω,
        :Γ => Γ,
        :Γbot => Γbot,
        :Γin => Γin,
        :Γout => Γout,
        :Γ_structures => Γ_structures,
        :Γ_dampings => Γ_dampings,
        :Γfs => Γfs,
        :Γκ => Γκ,
        :Γη => Γη,
    )
    # Add each structure triangulation under its domain_symbol, error if duplicate
    for sym in structure_syms
        if haskey(trian_dict, sym)
            error("domain_symbol :$sym for structure triangulation would overwrite an existing triangulation key in TankTriangulations. All keys must be unique.")
        end
        trian_dict[sym] = structure_trians[sym]
    end
    # Add each damping triangulation under its domain_symbol, error if duplicate
    for sym in damping_syms
        if haskey(trian_dict, sym)
            error("domain_symbol :$sym for damping triangulation would overwrite an existing triangulation key in TankTriangulations. All keys must be unique.")
        end
        trian_dict[sym] = damping_trians[sym]
    end
    return TankTriangulations(trian_dict)
end

"""
    _or_bits(bits, n) -> BitVector

Element-wise OR of a vector of boolean arrays.
Returns a `falses(n)` vector when `bits` is empty.
"""
function _or_bits(bits, n)
    result = falses(n)
    for b in bits
        result .= result .| b
    end
    return result
end


# ─────────────────────────────────────────────────────────────
# IntegrationDomains bridge
# ─────────────────────────────────────────────────────────────

"""
    get_integration_domains(tri::TankTriangulations; degree::Int=4) -> IntegrationDomains

Build an `IntegrationDomains` container of integration measures and normals
from `tri`.

# Populated keys

| Key | Source | Used by |
|-----|--------|---------|
| `:dΩ` | `tri.Ω` | PotentialFlow, Resonator |
| `:dΓfs` | `tri.Γfs` | Free surface without structures or damping zones |
| `:dΓ_s` | `tri.Γη` (all structures) | Membrane2D, Beam, PF↔struct coupling |
| `:dΓ_in` | `tri.Γin` | Inlet boundary terms |
| `:dΓ_out` | `tri.Γout` | Outlet boundary terms |
| `:dΓ_bot` | `tri.Γbot` | Bottom boundary terms |
| `:dΓ_d_i` | `tri.Γ_dampings[i]` | Damping zone terms (`:dΓ_d_1`, `:dΓ_d_2`, …) |

# Optional skeleton keys (beams only)

If the structure triangulation has interior edges,
add `:dΛ_s`, `:n_Λ_s`, and `:h_s` manually after calling this function.

# Example

```julia
dom = G.get_integration_domains(tank_trians; degree=4)
```
"""
function get_integration_domains(tri::TankTriangulations; degree::Union{Int, Dict{Symbol, Int}}=4)

    d = Dict{Symbol, Any}()

    # Helper to get degree for a key
    get_deg(key) = isa(degree, Dict) ? get(degree, key, 4) : degree

    # Fluid interior
    d[:dΩ]     = Measure(tri[:Ω], get_deg(:dΩ))

    # Free surface without structures or damping zones
    d[:dΓfs] = Measure(tri[:Γfs], get_deg(:dΓfs))
    d[:nΓfs] = get_normal_vector(tri[:Γfs])

    # Free surface including damping zones, but excluding structures
    d[:dΓκ]  = Measure(tri[:Γκ], get_deg(:dΓκ))

    # All-structure surface
    d[:dΓη]   = Measure(tri[:Γη], get_deg(:dΓη))

    # Per-structure measures
    for (i, Γs) in enumerate(tri[:Γ_structures])
        key = Symbol("dΓη_$i")
        d[key] = Measure(Γs, get_deg(key))
    end

    # Walls
    d[:dΓin]  = Measure(tri[:Γin], get_deg(:dΓin))
    d[:dΓout] = Measure(tri[:Γout], get_deg(:dΓout))
    d[:dΓbot] = Measure(tri[:Γbot], get_deg(:dΓbot))

    # Per-damping-zone measures (:dΓd_1, :dΓd_2, …)
    for (i, Γd) in enumerate(tri[:Γ_dampings])
        key = Symbol("dΓd_$i")
        d[key] = Measure(Γd, get_deg(key))
        d[Symbol("nΓd_$i")] = get_normal_vector(Γd)
    end

    return IntegrationDomains(d)
end
