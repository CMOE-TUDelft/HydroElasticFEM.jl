"""
    module FEOperators

Generic weak form assembler for HydroElasticFEM.

Provides `assemble_*` functions that loop over a collection of physics
terms and sum contributions.  Field tuples are wrapped in `FieldMap`
for symbol-based access.

Also provides coupling-detection and multi-entity assembly helpers used
by the `simulate` orchestrator.

All concrete methods (`mass`, `damping`, `stiffness`, `rhs`,
`weakform`, `residual`, `jacobian`, `jacobian_t`, `jacobian_tt`)
are defined in Entities (inside each entity file);
this module only provides composition logic and the field mapping.
"""
module FEOperators

import ...PhysicsCore.Entities as E
import ...Geometry as G

# ─────────────────────────────────────────────────────────────
# FieldMap — symbol-indexed wrapper for FE field tuples
# ─────────────────────────────────────────────────────────────

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

# Internal: wrap raw tuples in FieldMap
_wrap(t, fmap::Dict{Symbol,Int}) = FieldMap(t, fmap)

# ─────────────────────────────────────────────────────────────
# Helpers: sum only active form contributions (single-entity)
# ─────────────────────────────────────────────────────────────

function _assemble_active(f, has_form, terms, args...)
    val = nothing
    for term in terms
        if has_form(term)
            contrib = f(term, args...)
            val = isnothing(val) ? contrib : (val + contrib)
        end
    end
    isnothing(val) && error("No active $(nameof(f)) contributions in provided terms")
    return val
end

_has_weakform(term) =
    E.has_mass_form(term) || E.has_damping_form(term) || E.has_stiffness_form(term)
_has_weakform(a, b) =
    E.has_mass_form(a, b) || E.has_damping_form(a, b) || E.has_stiffness_form(a, b)
_has_residual(term) = _has_weakform(term) || E.has_rhs_form(term)

# ─────────────────────────────────────────────────────────────
# Coupling detection
# ─────────────────────────────────────────────────────────────

"""
    detect_couplings(entities) -> Vector{Tuple}

Auto-detect active coupling pairs among `entities` by checking all
ordered pairs `(a, b)` against the trait functions `has_mass_form`,
`has_damping_form`, `has_stiffness_form`, and `has_rhs_form`.

Returns a vector of `(entity_a, entity_b)` tuples that have at
least one active coupling form.
"""
function detect_couplings(entities)
    pairs = Tuple[]
    for (i, ea) in enumerate(entities)
        for (j, eb) in enumerate(entities)
            i == j && continue
            if (E.has_mass_form(ea, eb) || E.has_damping_form(ea, eb) ||
                E.has_stiffness_form(ea, eb) || E.has_rhs_form(ea, eb))
                push!(pairs, (ea, eb))
            end
        end
    end
    return pairs
end

# ─────────────────────────────────────────────────────────────
# Multi-entity assembly helpers (single + coupling entities)
# ─────────────────────────────────────────────────────────────

_add(a, b) = isnothing(a) ? b : a + b

"""Sum single-entity + coupling contributions for a given form."""
function _assemble_form(f, has_f, entities, coupling_pairs, dom, fmap, x, y)
    xd = FieldMap(x, fmap)
    yd = FieldMap(y, fmap)
    val = nothing
    for e in entities
        if has_f(e)
            val = _add(val, f(e, dom, xd, yd))
        end
    end
    for (ea, eb) in coupling_pairs
        if has_f(ea, eb)
            val = _add(val, f(ea, eb, dom, xd, yd))
        end
    end
    isnothing(val) && error("No active $(nameof(f)) contributions found")
    return val
end

"""Assemble frequency-domain bilinear form (single + coupling entities)."""
function _assemble_bilinear(entities, coupling_pairs, dom, ω, fmap, x, y)
    xd = FieldMap(x, fmap)
    yd = FieldMap(y, fmap)
    val = nothing
    for e in entities
        if _has_weakform(e)
            val = _add(val, E.weakform(e, dom, ω, xd, yd))
        end
    end
    for (ea, eb) in coupling_pairs
        if _has_weakform(ea, eb)
            val = _add(val, E.weakform(ea, eb, dom, ω, xd, yd))
        end
    end
    isnothing(val) && error("No active weak form contributions found")
    return val
end

# ─────────────────────────────────────────────────────────────
# Linear form assemblers
# ─────────────────────────────────────────────────────────────

"""
    assemble_mass(terms, dom, fmap, x_tt, y)

Sum `mass(term, dom, x_tt, y)` over all `terms`,
wrapping the raw tuples with `fmap`.
"""
function assemble_mass(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int}, x_tt, y)
    xd = _wrap(x_tt, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(E.mass, E.has_mass_form, terms, dom, xd, yd)
end

"""
    assemble_damping(terms, dom, fmap, x_t, y)

Sum `damping(term, dom, x_t, y)` over all `terms`.
"""
function assemble_damping(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int}, x_t, y)
    xd = _wrap(x_t, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(E.damping, E.has_damping_form, terms, dom, xd, yd)
end

"""
    assemble_stiffness(terms, dom, fmap, x, y)

Sum `stiffness(term, dom, x, y)` over all `terms`.
"""
function assemble_stiffness(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int}, x, y)
    xd = _wrap(x, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(E.stiffness, E.has_stiffness_form, terms, dom, xd, yd)
end

"""
    assemble_rhs(terms, dom, fmap, f, y)

Sum `rhs(term, dom, f, y)` over all `terms`.
"""
function assemble_rhs(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int}, f, y)
    fd = _wrap(f, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(E.rhs, E.has_rhs_form, terms, dom, fd, yd)
end

# ─────────────────────────────────────────────────────────────
# Frequency-domain assembler
# ─────────────────────────────────────────────────────────────

"""
    assemble_weakform(terms, dom, ω, fmap, x, y)

Sum `weakform(term, dom, ω, x, y)` over all `terms`,
wrapping the raw tuples with `fmap`.
"""
function assemble_weakform(terms, dom::G.IntegrationDomains, ω, fmap::Dict{Symbol,Int}, x, y)
    xd = _wrap(x, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(E.weakform, _has_weakform, terms, dom, ω, xd, yd)
end

# ─────────────────────────────────────────────────────────────
# Nonlinear form assemblers
# ─────────────────────────────────────────────────────────────

"""
    assemble_residual(terms, dom, fmap, x, x_t, x_tt, f, y)
"""
function assemble_residual(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int},
                           x, x_t, x_tt, f, y)
    xd    = _wrap(x, fmap)
    xd_t  = _wrap(x_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    fd    = _wrap(f, fmap)
    yd    = _wrap(y, fmap)
    _assemble_active(residual, _has_residual, terms, dom, xd, xd_t, xd_tt, fd, yd)
end

"""
    assemble_jacobian(terms, dom, fmap, dx, x_t, x_tt, y)
"""
function assemble_jacobian(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int},
                           dx, x_t, x_tt, y)
    dxd   = _wrap(dx, fmap)
    xd_t  = _wrap(x_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    yd    = _wrap(y, fmap)
    _assemble_active(
        jacobian, has_stiffness_form,
        terms, dom, dxd, xd_t, xd_tt, yd,
    )
end

"""
    assemble_jacobian_t(terms, dom, fmap, x, dx_t, x_tt, y)
"""
function assemble_jacobian_t(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int},
                             x, dx_t, x_tt, y)
    xd    = _wrap(x, fmap)
    dxd_t = _wrap(dx_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    yd    = _wrap(y, fmap)
    _assemble_active(
        jacobian_t, has_damping_form,
        terms, dom, xd, dxd_t, xd_tt, yd,
    )
end

"""
    assemble_jacobian_tt(terms, dom, fmap, x, x_t, dx_tt, y)
"""
function assemble_jacobian_tt(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int},
                              x, x_t, dx_tt, y)
    xd    = _wrap(x, fmap)
    xd_t  = _wrap(x_t, fmap)
    dxd_tt = _wrap(dx_tt, fmap)
    yd    = _wrap(y, fmap)
    _assemble_active(
        jacobian_tt, has_mass_form,
        terms, dom, xd, xd_t, dxd_tt, yd,
    )
end

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export FieldMap
export detect_couplings
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module FEOperators
