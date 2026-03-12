"""
    module Geometry

Geometry, mesh construction, integration domains and field helpers.

Provides:
- `TankDomain2D`, `StructureDomain1D`, `DampingZone1D` — geometry specs
- `build_model`, `build_triangulations` — mesh construction
- `IntegrationDomains` — dict-based container for Gridap measures/normals
- `get_integration_domains` — build `IntegrationDomains` from triangulations
"""
module Geometry
using Parameters
using Gridap


"""
    TankTriangulations(; key=value, ...)
    TankTriangulations(dict::Dict{Symbol,Any})

Dict-based container for triangulations.
Allows flexible access to triangulations by symbol keys.

# Standard key conventions (not enforced)
- `:Ω`  — Interior (fluid domain)
- `:Γ`  — Full top-surface boundary
- `:Γbot` — Bottom boundary
- `:Γin` — Inlet (left wall) boundary
- `:Γout` — Outlet (right wall) boundary
- `:Γ_structures` — `Vector`: one triangulation per structure domain, ordered as in `TankDomain2D.structure_domains`
- `:Γ_dampings`   — `Vector`: one triangulation per damping zone, ordered as in `TankDomain2D.damping_zones`
- `:Γfs` — Free surface: surface cells that belong to no structure and no damping zone
- `:Γκ`  — Non-structure surface (free surface ∪ damping zones)
- `:Γη`  — All-structure surface (union of all structure triangulations)
"""
struct TankTriangulations
    data::Dict{Symbol, Any}
end

TankTriangulations(; kwargs...) =
    TankTriangulations(Dict{Symbol, Any}(k => v for (k, v) in pairs(kwargs)))

Base.getindex(t::TankTriangulations, k::Symbol)            = t.data[k]
Base.haskey(t::TankTriangulations, k::Symbol)              = haskey(t.data, k)
Base.get(t::TankTriangulations, k::Symbol, default)        = get(t.data, k, default)
Base.setindex!(t::TankTriangulations, val, k::Symbol)      = (t.data[k] = val)
Base.keys(t::TankTriangulations)                           = keys(t.data)


"""
    IntegrationDomains(; key=value, ...)
    IntegrationDomains(dict::Dict{Symbol,Any})

Dict-based container for Gridap measures, normals, DiracDeltas,
and any other domain data needed by weak form methods.

Each `weakform` dispatch accesses only the keys it needs via
`dom[:key]`.  No fixed schema — new keys can be added without
changing this type.

# Standard key conventions (not enforced)
- `:dΩ`     — fluid interior measure
- `:dΓ_fs`  — free-surface measure (outside structure)
- `:dΓ_s`   — structure surface measure
- `:dΛ_s`, `:n_Λ_s`, `:h_s` — beam skeleton measures/normals + mesh size
- `:dΛ_sb`, `:n_Λ_sb`  — structure boundary (fixed BC Neumann)
- `:dΓ_in`, `:dΓ_out`  — inlet / outlet radiation boundaries
- `:dΓ_d_1`, `:dΓ_d_2` — damping zone measures
- `:δ_p`    — vector of DiracDelta functionals (resonator points)
"""
struct IntegrationDomains
    data::Dict{Symbol, Any}
end

IntegrationDomains(; kwargs...) =
    IntegrationDomains(Dict{Symbol, Any}(k => v for (k, v) in pairs(kwargs)))

Base.getindex(d::IntegrationDomains, k::Symbol)            = d.data[k]
Base.haskey(d::IntegrationDomains, k::Symbol)               = haskey(d.data, k)
Base.get(d::IntegrationDomains, k::Symbol, default)         = get(d.data, k, default)
Base.setindex!(d::IntegrationDomains, val, k::Symbol)       = (d.data[k] = val)
Base.keys(d::IntegrationDomains)                            = keys(d.data)

include("CartesianGeometry.jl")

export TankTriangulations
export IntegrationDomains

end # module Geometry