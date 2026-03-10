"""
    module Geometry

Geometry, mesh construction, integration domains and field helpers.

Provides:
- `TankDomain2D`, `StructureDomain1D`, `DampingZone1D` — geometry specs
- `build_model`, `build_triangulations` — mesh construction
- `IntegrationDomains` — dict-based container for Gridap measures/normals
- `FieldMap` — symbol-indexed wrapper around Gridap field tuples
- `get_integration_domains` — build `IntegrationDomains` from triangulations
"""
module Geometry
using Parameters
using Gridap

include("CartesianGeometry.jl")

"""
    FieldMap{T}

Wraps a positional tuple of FE fields (from Gridap's multi-field
decomposition) and maps `Symbol` keys to positional indices. It allows
symbol-based access to fields, which is more intuitive when writing weak
form contributions that involve multiple fields.

# Usage
```julia
fmap = Dict(:ϕ => 1, :κ => 2, :η_m => 3)
x = FieldMap((ϕ, κ, η), fmap)
x[:ϕ]   # returns ϕ
x[:η_m] # returns η
```
"""
struct FieldMap{T}
    _data::T
    _map::Dict{Symbol, Int}
end

Base.getindex(fd::FieldMap, s::Symbol) = fd._data[fd._map[s]]
Base.haskey(fd::FieldMap, s::Symbol)   = haskey(fd._map, s)
Base.keys(fd::FieldMap)                = keys(fd._map)


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


export IntegrationDomains
export FieldMap

end # module Geometry