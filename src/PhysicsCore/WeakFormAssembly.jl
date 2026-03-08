"""
    module WeakFormAssembly

Generic weak form assembler for HydroElasticFEM.

Provides the `FieldDict` wrapper for symbol-based field access,
and `assemble_*` functions that loop over a collection of physics
terms and sum contributions.

All concrete methods (`mass`, `damping`, `stiffness`, `rhs`,
`weakform`, `residual`, `jacobian`, `jacobian_t`, `jacobian_tt`)
are defined in PhysicalEntities (inside each entity file);
this module only provides composition logic and the field mapping.
"""
module WeakFormAssembly

using ..PhysicalEntities

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

# ─────────────────────────────────────────────────────────────
# Helper: sum a given form function over all terms
# ─────────────────────────────────────────────────────────────

function _assemble(f, terms, args...)
    val = f(first(terms), args...)
    for i in 2:length(terms)
        val += f(terms[i], args...)
    end
    return val
end

# ─────────────────────────────────────────────────────────────
# Internal: wrap raw tuples in FieldDict
# ─────────────────────────────────────────────────────────────

_wrap(t, fmap::Dict{Symbol,Int}) = FieldDict(t, fmap)

# ─────────────────────────────────────────────────────────────
# Linear form assemblers
# ─────────────────────────────────────────────────────────────

"""
    assemble_mass(terms, dom, fmap, x_tt, y)

Sum `mass(term, dom, x_tt, y)` over all `terms`,
wrapping the raw tuples with `fmap`.
"""
function assemble_mass(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int}, x_tt, y)
    xd = _wrap(x_tt, fmap)
    yd = _wrap(y, fmap)
    _assemble(mass, terms, dom, xd, yd)
end

"""
    assemble_damping(terms, dom, fmap, x_t, y)

Sum `damping(term, dom, x_t, y)` over all `terms`.
"""
function assemble_damping(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int}, x_t, y)
    xd = _wrap(x_t, fmap)
    yd = _wrap(y, fmap)
    _assemble(damping, terms, dom, xd, yd)
end

"""
    assemble_stiffness(terms, dom, fmap, x, y)

Sum `stiffness(term, dom, x, y)` over all `terms`.
"""
function assemble_stiffness(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int}, x, y)
    xd = _wrap(x, fmap)
    yd = _wrap(y, fmap)
    _assemble(stiffness, terms, dom, xd, yd)
end

"""
    assemble_rhs(terms, dom, fmap, f, y)

Sum `rhs(term, dom, f, y)` over all `terms`.
"""
function assemble_rhs(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int}, f, y)
    fd = _wrap(f, fmap)
    yd = _wrap(y, fmap)
    _assemble(rhs, terms, dom, fd, yd)
end

# ─────────────────────────────────────────────────────────────
# Frequency-domain assembler
# ─────────────────────────────────────────────────────────────

"""
    assemble_weakform(terms, dom, ω, fmap, x, y)

Sum `weakform(term, dom, ω, x, y)` over all `terms`,
wrapping the raw tuples with `fmap`.
"""
function assemble_weakform(terms, dom::WeakFormDomains, ω, fmap::Dict{Symbol,Int}, x, y)
    xd = _wrap(x, fmap)
    yd = _wrap(y, fmap)
    _assemble(weakform, terms, dom, ω, xd, yd)
end

# ─────────────────────────────────────────────────────────────
# Nonlinear form assemblers
# ─────────────────────────────────────────────────────────────

"""
    assemble_residual(terms, dom, fmap, x, x_t, x_tt, f, y)
"""
function assemble_residual(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int},
                           x, x_t, x_tt, f, y)
    xd    = _wrap(x, fmap)
    xd_t  = _wrap(x_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    fd    = _wrap(f, fmap)
    yd    = _wrap(y, fmap)
    _assemble(residual, terms, dom, xd, xd_t, xd_tt, fd, yd)
end

"""
    assemble_jacobian(terms, dom, fmap, dx, x_t, x_tt, y)
"""
function assemble_jacobian(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int},
                           dx, x_t, x_tt, y)
    dxd   = _wrap(dx, fmap)
    xd_t  = _wrap(x_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    yd    = _wrap(y, fmap)
    _assemble(jacobian, terms, dom, dxd, xd_t, xd_tt, yd)
end

"""
    assemble_jacobian_t(terms, dom, fmap, x, dx_t, x_tt, y)
"""
function assemble_jacobian_t(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int},
                             x, dx_t, x_tt, y)
    xd    = _wrap(x, fmap)
    dxd_t = _wrap(dx_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    yd    = _wrap(y, fmap)
    _assemble(jacobian_t, terms, dom, xd, dxd_t, xd_tt, yd)
end

"""
    assemble_jacobian_tt(terms, dom, fmap, x, x_t, dx_tt, y)
"""
function assemble_jacobian_tt(terms, dom::WeakFormDomains, fmap::Dict{Symbol,Int},
                              x, x_t, dx_tt, y)
    xd    = _wrap(x, fmap)
    xd_t  = _wrap(x_t, fmap)
    dxd_tt = _wrap(dx_tt, fmap)
    yd    = _wrap(y, fmap)
    _assemble(jacobian_tt, terms, dom, xd, xd_t, dxd_tt, yd)
end

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export FieldDict
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module WeakFormAssembly
