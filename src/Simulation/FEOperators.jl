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
are defined in Physics (inside each entity file);
this module only provides composition logic and the field mapping.
"""
module FEOperators

import ...Physics as P
import ...Geometry as G
import ...AssemblyContexts as AC

using Gridap
using Gridap.ODEs

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
    P.has_mass_form(term) || P.has_damping_form(term) || P.has_stiffness_form(term)
_has_weakform(a, b) =
    P.has_mass_form(a, b) || P.has_damping_form(a, b) || P.has_stiffness_form(a, b)
_has_residual(term) = _has_weakform(term) || P.has_rhs_form(term)

_has_weakform(ctx::AC.AbstractAssemblyContext, term) = begin
    forms = P.active_forms(ctx, term)
    forms.mass || forms.damping || forms.stiffness
end

_has_weakform(ctx::AC.AbstractAssemblyContext, a, b) = begin
    forms = P.active_forms(ctx, a, b)
    forms.mass || forms.damping || forms.stiffness
end

_has_residual(ctx::AC.AbstractAssemblyContext, term) = begin
    forms = P.active_forms(ctx, term)
    forms.mass || forms.damping || forms.stiffness || forms.rhs
end

_has_residual(ctx::AC.AbstractAssemblyContext, a, b) = begin
    forms = P.active_forms(ctx, a, b)
    forms.mass || forms.damping || forms.stiffness || forms.rhs
end

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
            if (P.has_mass_form(ea, eb) || P.has_damping_form(ea, eb) ||
                P.has_stiffness_form(ea, eb) || P.has_rhs_form(ea, eb))
                push!(pairs, (ea, eb))
            end
        end
    end
    return pairs
end

function detect_couplings(entities, ctx::AC.AbstractAssemblyContext)
    pairs = Tuple[]
    for (i, ea) in enumerate(entities)
        for (j, eb) in enumerate(entities)
            i == j && continue
            if _has_residual(ctx, ea, eb)
                push!(pairs, (ea, eb))
            end
        end
    end
    return pairs
end

# ─────────────────────────────────────────────────────────────
# Multi-entity assembly helpers (single + coupling entities)
# ─────────────────────────────────────────────────────────────

"""Sum single-entity + coupling contributions for a given form."""
function _assemble_form(form_sym::Symbol, f, entities, coupling_pairs, ctx::AC.AbstractAssemblyContext, fmap, x, y)
    xd = FieldMap(x, fmap)
    yd = FieldMap(y, fmap)
    val = nothing
    for e in entities
        if getproperty(P.active_forms(ctx, e), form_sym)
            val = P._add_contribution(val, f(e, ctx, xd, yd))
        end
    end
    for (ea, eb) in coupling_pairs
        if getproperty(P.active_forms(ctx, ea, eb), form_sym)
            val = P._add_contribution(val, f(ea, eb, ctx, xd, yd))
        end
    end
    isnothing(val) && error("No active $(nameof(f)) contributions found")
    return val
end

function _has_active_form(form_sym::Symbol, entities, coupling_pairs, ctx::AC.AbstractAssemblyContext)
    single_active = any(e -> getproperty(P.active_forms(ctx, e), form_sym), entities)
    pair_active = any(((ea, eb),) -> getproperty(P.active_forms(ctx, ea, eb), form_sym), coupling_pairs)
    return single_active || pair_active
end

function _find_volume_symbol(entities, fmap)
    for e in entities
        sym = P.variable_symbol(e)
        if haskey(fmap, sym) && e.space_domain_symbol == :Ω
            return sym
        end
    end
    return nothing
end

function _zero_mass_contribution(ctx::AC.TimeAssemblyContext, fmap, x_tt, y, sym::Symbol)
    haskey(AC.domains(ctx), :dΩ) || error("Cannot build zero mass fallback: domain measure `:dΩ` is missing.")
    xd = FieldMap(x_tt, fmap)
    yd = FieldMap(y, fmap)
    return ∫(0.0 * xd[sym] * yd[sym])AC.domains(ctx)[:dΩ]
end

"""Assemble frequency-domain bilinear form (single + coupling entities)."""
function _assemble_bilinear(entities, coupling_pairs, ctx::AC.FrequencyAssemblyContext, fmap, x, y)
    xd = FieldMap(x, fmap)
    yd = FieldMap(y, fmap)
    val = nothing
    for e in entities
        if _has_weakform(ctx, e)
            val = P._add_contribution(val, P.weakform(e, ctx, xd, yd))
        end
    end
    for (ea, eb) in coupling_pairs
        if _has_weakform(ctx, ea, eb)
            val = P._add_contribution(val, P.weakform(ea, eb, ctx, xd, yd))
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
    _assemble_active(P.mass, P.has_mass_form, terms, dom, xd, yd)
end

function assemble_mass(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int}, x_tt, y)
    _assemble_form(:mass, P.mass, terms, Tuple[], ctx, fmap, x_tt, y)
end

"""
    assemble_damping(terms, dom, fmap, x_t, y)

Sum `damping(term, dom, x_t, y)` over all `terms`.
"""
function assemble_damping(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int}, x_t, y)
    xd = _wrap(x_t, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(P.damping, P.has_damping_form, terms, dom, xd, yd)
end

function assemble_damping(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int}, x_t, y)
    _assemble_form(:damping, P.damping, terms, Tuple[], ctx, fmap, x_t, y)
end

"""
    assemble_stiffness(terms, dom, fmap, x, y)

Sum `stiffness(term, dom, x, y)` over all `terms`.
"""
function assemble_stiffness(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int}, x, y)
    xd = _wrap(x, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(P.stiffness, P.has_stiffness_form, terms, dom, xd, yd)
end

function assemble_stiffness(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int}, x, y)
    _assemble_form(:stiffness, P.stiffness, terms, Tuple[], ctx, fmap, x, y)
end

"""
    assemble_rhs(terms, dom, fmap, f, y)

Sum `rhs(term, dom, f, y)` over all `terms`.
"""
function assemble_rhs(terms, dom::G.IntegrationDomains, fmap::Dict{Symbol,Int}, f, y)
    fd = _wrap(f, fmap)
    yd = _wrap(y, fmap)
    _assemble_active(P.rhs, P.has_rhs_form, terms, dom, fd, yd)
end

function assemble_rhs(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int}, f, y)
    _assemble_form(:rhs, P.rhs, terms, Tuple[], ctx, fmap, f, y)
end

function _assemble_rhs_total(entities, coupling_pairs, ctx::AC.AbstractAssemblyContext, fmap, f, y)
    fd = _wrap(f, fmap)
    yd = _wrap(y, fmap)
    val = nothing
    for e in entities
        if P.active_forms(ctx, e).rhs
            val = P._add_contribution(val, P.rhs(e, ctx, fd, yd))
        end
    end
    for (ea, eb) in coupling_pairs
        if P.active_forms(ctx, ea, eb).rhs
            val = P._add_contribution(val, P.rhs(ea, eb, ctx, fd, yd))
        end
    end
    isnothing(val) && error("No active rhs contributions found")
    return val
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
    _assemble_active(P.weakform, _has_weakform, terms, dom, ω, xd, yd)
end

function assemble_weakform(terms, ctx::AC.FrequencyAssemblyContext, fmap::Dict{Symbol,Int}, x, y)
    xd = _wrap(x, fmap)
    yd = _wrap(y, fmap)
    val = nothing
    for term in terms
        if _has_weakform(ctx, term)
            val = P._add_contribution(val, P.weakform(term, ctx, xd, yd))
        end
    end
    isnothing(val) && error("No active weak form contributions in provided terms")
    return val
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
    _assemble_active(P.residual, _has_residual, terms, dom, xd, xd_t, xd_tt, fd, yd)
end

function assemble_residual(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int},
                           x, x_t, x_tt, f, y)
    xd    = _wrap(x, fmap)
    xd_t  = _wrap(x_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    fd    = _wrap(f, fmap)
    yd    = _wrap(y, fmap)
    val = nothing
    for term in terms
        if _has_residual(ctx, term)
            val = P._add_contribution(val, P.residual(term, ctx, xd, xd_t, xd_tt, fd, yd))
        end
    end
    isnothing(val) && error("No active residual contributions in provided terms")
    return val
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
        P.jacobian, P.has_stiffness_form,
        terms, dom, dxd, xd_t, xd_tt, yd,
    )
end

function assemble_jacobian(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int},
                           dx, x_t, x_tt, y)
    dxd   = _wrap(dx, fmap)
    xd_t  = _wrap(x_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    yd    = _wrap(y, fmap)
    val = nothing
    for term in terms
        if P.active_forms(ctx, term).stiffness
            val = P._add_contribution(val, P.jacobian(term, ctx, dxd, xd_t, xd_tt, yd))
        end
    end
    isnothing(val) && error("No active jacobian contributions in provided terms")
    return val
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
        P.jacobian_t, P.has_damping_form,
        terms, dom, xd, dxd_t, xd_tt, yd,
    )
end

function assemble_jacobian_t(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int},
                             x, dx_t, x_tt, y)
    xd    = _wrap(x, fmap)
    dxd_t = _wrap(dx_t, fmap)
    xd_tt = _wrap(x_tt, fmap)
    yd    = _wrap(y, fmap)
    val = nothing
    for term in terms
        if P.active_forms(ctx, term).damping
            val = P._add_contribution(val, P.jacobian_t(term, ctx, xd, dxd_t, xd_tt, yd))
        end
    end
    isnothing(val) && error("No active jacobian_t contributions in provided terms")
    return val
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
        P.jacobian_tt, P.has_mass_form,
        terms, dom, xd, xd_t, dxd_tt, yd,
    )
end

function assemble_jacobian_tt(terms, ctx::AC.AbstractAssemblyContext, fmap::Dict{Symbol,Int},
                              x, x_t, dx_tt, y)
    xd    = _wrap(x, fmap)
    xd_t  = _wrap(x_t, fmap)
    dxd_tt = _wrap(dx_tt, fmap)
    yd    = _wrap(y, fmap)
    val = nothing
    for term in terms
        if P.active_forms(ctx, term).mass
            val = P._add_contribution(val, P.jacobian_tt(term, ctx, xd, xd_t, dxd_tt, yd))
        end
    end
    isnothing(val) && error("No active jacobian_tt contributions in provided terms")
    return val
end

function _adapt_frequency_rhs(rhs_fn)
    (ctx, y) -> (applicable(rhs_fn, ctx, y) && !applicable(rhs_fn, y)) ? rhs_fn(ctx, y) : rhs_fn(y)
end

function _adapt_time_rhs(rhs_fn)
    (ctx, y) -> applicable(rhs_fn, AC.current_time(ctx), y) ? rhs_fn(AC.current_time(ctx), y) : rhs_fn(ctx, y)
end

_zero_rhs(fmap) = (ctx, y) -> zeros(length(fmap))

# ─────────────────────────────────────────────────────────────
# FE operator construction
# ─────────────────────────────────────────────────────────────

"""
    build_fe_operator(entities, coupling_pairs, dom, ω, fmap, X, Y; rhs_fn=nothing)

Build a **frequency-domain** `AffineFEOperator`.

# Arguments
- `entities` — physics entities
- `dom::IntegrationDomains` — integration measures/normals
- `ω` — angular frequency
- `fmap::Dict{Symbol,Int}` — field-symbol to positional-index map
- `X` — trial multi-field FE space
- `Y` — test multi-field FE space
- `rhs_fn` — optional callable `rhs_fn(y::FieldMap) -> DomainContribution`;
  if `nothing`, uses zero right-hand side (requires `:dΩ` in `dom`)
"""
function build_frequency_fe_operator(entities::Vector{<:P.PhysicsParameters},
                                     ctx::AC.FrequencyAssemblyContext,
                                     fmap::Dict{Symbol,Int}, X, Y;
                                     rhs_fn=nothing)
    coupling_pairs = detect_couplings(entities, ctx)
    rhs_cb = rhs_fn === nothing ? _zero_rhs(fmap) : _adapt_frequency_rhs(rhs_fn)

    a(x, y) = _assemble_bilinear(entities, coupling_pairs, ctx, fmap, x, y)
    l(y) = _assemble_rhs_total(entities, coupling_pairs, ctx, fmap, rhs_cb(ctx, FieldMap(y, fmap)), y)

    AffineFEOperator(a, l, X, Y)
end

"""
    build_fe_operator(entities, coupling_pairs, dom, fmap, X, Y; rhs_fn=nothing)

Build a **time-domain** `TransientLinearFEOperator`.

# Arguments
- `entities::Vector{<:P.PhysicsParameters}` — physics entities
- `dom::IntegrationDomains` — integration measures/normals
- `fmap::Dict{Symbol,Int}` — field-symbol to positional-index map
- `X` — trial multi-field FE space (transient)
- `Y` — test multi-field FE space
- `rhs_fn` — optional callable `rhs_fn(t, y::FieldMap) -> DomainContribution`;
  if `nothing`, uses zero right-hand side (requires `:dΩ` in `dom`)
"""
function build_time_fe_operator(entities::Vector{<:P.PhysicsParameters},
                                base_ctx::AC.TimeAssemblyContext,
                                fmap::Dict{Symbol,Int}, X, Y;
                                rhs_fn=nothing)
    coupling_pairs = detect_couplings(entities, base_ctx)
    rhs_cb = rhs_fn === nothing ? _zero_rhs(fmap) : _adapt_time_rhs(rhs_fn)
    volume_sym = _find_volume_symbol(entities, fmap)

    a(t, x, y) = _assemble_form(:stiffness, P.stiffness, entities, coupling_pairs, AC.with_time(base_ctx, t), fmap, x, y)
    c(t, x_t, y) = _assemble_form(:damping, P.damping, entities, coupling_pairs, AC.with_time(base_ctx, t), fmap, x_t, y)
    m(t, x_tt, y) = begin
        ctx = AC.with_time(base_ctx, t)
        if _has_active_form(:mass, entities, coupling_pairs, ctx)
            _assemble_form(:mass, P.mass, entities, coupling_pairs, ctx, fmap, x_tt, y)
        else
            isnothing(volume_sym) && error("No active mass contributions found and no volume field is available for zero-mass fallback.")
            _zero_mass_contribution(ctx, fmap, x_tt, y, volume_sym)
        end
    end
    l(t, y) = begin
        ctx = AC.with_time(base_ctx, t)
        _assemble_rhs_total(entities, coupling_pairs, ctx, fmap, rhs_cb(ctx, FieldMap(y, fmap)), y)
    end

    TransientLinearFEOperator((a, c, m), l, X, Y;
        constant_forms=(true, true, true))
end

function build_fe_operator(entities,
                           dom::G.IntegrationDomains, ω,
                           fmap::Dict{Symbol,Int}, X, Y;
                           rhs_fn=nothing)
    ctx = AC.FrequencyAssemblyContext(dom, ω, nothing)
    build_frequency_fe_operator(entities, ctx, fmap, X, Y; rhs_fn=rhs_fn)
end

function build_fe_operator(entities::Vector{<:P.PhysicsParameters},
                           dom::G.IntegrationDomains,
                           fmap::Dict{Symbol,Int}, X, Y;
                           rhs_fn=nothing)
    ctx = AC.TimeAssemblyContext(dom, 0.0, nothing)
    build_time_fe_operator(entities, ctx, fmap, X, Y; rhs_fn=rhs_fn)
end

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export FieldMap
export detect_couplings
export build_fe_operator
export build_frequency_fe_operator, build_time_fe_operator
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module FEOperators
