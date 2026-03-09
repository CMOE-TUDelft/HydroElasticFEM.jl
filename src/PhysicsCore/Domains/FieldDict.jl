# ─────────────────────────────────────────────────────────────
# FieldDict: symbol-indexed wrapper around Gridap field tuples
# ─────────────────────────────────────────────────────────────

"""
    FieldDict{T}

Wraps a positional tuple of FE fields (from Gridap's multi-field
decomposition) and maps `Symbol` keys to positional indices.

# Usage
```julia
fmap = Dict(:ϕ => 1, :κ => 2, :η_m => 3)
x = FieldDict((ϕ, κ, η), fmap)
x[:ϕ]   # returns ϕ
x[:η_m] # returns η
```
"""
struct FieldDict{T}
    _data::T
    _map::Dict{Symbol, Int}
end

Base.getindex(fd::FieldDict, s::Symbol) = fd._data[fd._map[s]]
Base.haskey(fd::FieldDict, s::Symbol)   = haskey(fd._map, s)
Base.keys(fd::FieldDict)                = keys(fd._map)
