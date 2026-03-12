"""
    module Entities

Unified type hierarchy for all physical entities in HydroElasticFEM.

Each entity file may define any subset of `mass`, `damping`,
`stiffness`, and `rhs` methods that access FE fields by symbol
(e.g. `x[:ϕ]`, `x[:η_m]`) via a `FieldMap` wrapper provided by
`FEOperators`.

Generic composed forms (`weakform`, `residual`, `jacobian`, ...)
are defined at module level and dispatch to the linear forms.
"""
module Entities

using Parameters
using Gridap
using Printf
using ...Geometry
using ...ParameterHandler

# ─────────────────────────────────────────────────────────────
# Abstract base
# ─────────────────────────────────────────────────────────────

"""
    abstract type PhysicsParameters

Abstract base type for all physics parameter structures.
"""
abstract type PhysicsParameters end

"""
    print_parameters(params::PhysicsParameters)

Print the parameters of a concrete `PhysicsParameters` subtype.
Must be implemented for every concrete subtype.
"""
function print_parameters(params::PhysicsParameters)
    error("print_parameters not implemented for $(typeof(params))")
end

"""
    variable_symbol(s::PhysicsParameters) -> Symbol

Return the field symbol for the primary variable of entity `s`.
Must be implemented for every concrete subtype.
"""
function variable_symbol(s::PhysicsParameters)
    error("variable_symbol not implemented for $(typeof(s))")
end

# ─────────────────────────────────────────────────────────────
# Abstract weak form interface
# ─────────────────────────────────────────────────────────────

"""
    mass(s, dom, x_tt, y)

Mass weak form.
"""
function mass(s::PhysicsParameters, dom, x_tt, y)
    error("mass not implemented for $(typeof(s))")
end

"""
    damping(s, dom, x_t, y)

Damping weak form.
"""
function damping(s::PhysicsParameters, dom, x_t, y)
    error("damping not implemented for $(typeof(s))")
end

"""
    stiffness(s, dom, x, y)

Stiffness weak form.
"""
function stiffness(s::PhysicsParameters, dom, x, y)
    error("stiffness not implemented for $(typeof(s))")
end

"""
    rhs(s, dom, f, y)

Right-hand side weak form.
"""
function rhs(s::PhysicsParameters, dom, f, y)
    error("rhs not implemented for $(typeof(s))")
end

# Optional form-presence traits.
# Single-entity forms default to present for backward compatibility.
has_mass_form(::PhysicsParameters) = true
has_damping_form(::PhysicsParameters) = true
has_stiffness_form(::PhysicsParameters) = true
has_rhs_form(::PhysicsParameters) = true

# Two-entity coupling forms default to absent unless enabled.
has_mass_form(a, b) = false
has_damping_form(a, b) = false
has_stiffness_form(a, b) = false
has_rhs_form(a, b) = false

function _sum_present(terms...)
    val = nothing
    for t in terms
        if !isnothing(t)
            val = isnothing(val) ? t : (val + t)
        end
    end
    return val
end

_require_nonempty(val, kind, obj) =
    isnothing(val) ? error("$kind has no active contributions for $(typeof(obj))") : val


# ─────────────────────────────────────────────────────────────
# Structural entities
# ─────────────────────────────────────────────────────────────

abstract type AbstractStructure <: PhysicsParameters end

# Entity files (struct definition + single-variable weak forms)
include("PotentialFlow.jl")
include("FreeSurface.jl")
include("Membrane2D.jl")
include("EulerBernoulliBeam.jl")
include("Resonator.jl")

# Coupling weak forms (cross-terms between pairs of entities)
include("CouplingTerms.jl")

# ─────────────────────────────────────────────────────────────
# Frequency-domain weakform: derived from linear forms
# ─────────────────────────────────────────────────────────────

"""
    weakform(s, dom, ω, x, y)

Frequency-domain bilinear form composed from the linear `mass`,
`damping`, and `stiffness` weak forms via the substitution
`∂ₜₜ → -ω²`, `∂ₜ → -iω`:

    weakform = -ω² mass + (-iω) damping + stiffness
"""
weakform(s, dom::IntegrationDomains, ω, x, y) =
    _require_nonempty(
        _sum_present(
            has_mass_form(s) ? (-ω^2) * mass(s, dom, x, y) : nothing,
            has_damping_form(s) ? (-im * ω) * damping(s, dom, x, y) : nothing,
            has_stiffness_form(s) ? stiffness(s, dom, x, y) : nothing,
        ),
        "weakform", s
    )

# Two-entity coupling variant
weakform(a, b, dom::IntegrationDomains, ω, x, y) =
    _require_nonempty(
        _sum_present(
            has_mass_form(a, b) ? (-ω^2) * mass(a, b, dom, x, y) : nothing,
            has_damping_form(a, b) ? (-im * ω) * damping(a, b, dom, x, y) : nothing,
            has_stiffness_form(a, b) ? stiffness(a, b, dom, x, y) : nothing,
        ),
        "weakform", (a, b)
    )

# ─────────────────────────────────────────────────────────────
# Nonlinear weak forms: residual, jacobian, jacobian_t, jacobian_tt
# ─────────────────────────────────────────────────────────────

"""
    residual(s, dom, x, x_t, x_tt, f, y)

Nonlinear residual: `mass(x_tt, y) + damping(x_t, y) + stiffness(x, y) - rhs(f, y)`.
"""
residual(s, dom::IntegrationDomains, x, x_t, x_tt, f, y) =
    _require_nonempty(
        _sum_present(
            has_mass_form(s) ? mass(s, dom, x_tt, y) : nothing,
            has_damping_form(s) ? damping(s, dom, x_t, y) : nothing,
            has_stiffness_form(s) ? stiffness(s, dom, x, y) : nothing,
            has_rhs_form(s) ? (-rhs(s, dom, f, y)) : nothing,
        ),
        "residual", s
    )

# Two-entity coupling variant
residual(a, b, dom::IntegrationDomains, x, x_t, x_tt, f, y) =
    _require_nonempty(
        _sum_present(
            has_mass_form(a, b) ? mass(a, b, dom, x_tt, y) : nothing,
            has_damping_form(a, b) ? damping(a, b, dom, x_t, y) : nothing,
            has_stiffness_form(a, b) ? stiffness(a, b, dom, x, y) : nothing,
            has_rhs_form(a, b) ? (-rhs(a, b, dom, f, y)) : nothing,
        ),
        "residual", (a, b)
    )

"""
    jacobian(s, dom, dx, x_t, x_tt, y)

Jacobian w.r.t. displacement: `∂R/∂x · dx = stiffness(dx, y)`.
"""
jacobian(s, dom::IntegrationDomains, dx, x_t, x_tt, y) =
    _require_nonempty(
        has_stiffness_form(s) ? stiffness(s, dom, dx, y) : nothing,
        "jacobian", s
    )

# Two-entity coupling variant
jacobian(a, b, dom::IntegrationDomains, dx, x_t, x_tt, y) =
    _require_nonempty(
        has_stiffness_form(a, b) ? stiffness(a, b, dom, dx, y) : nothing,
        "jacobian", (a, b)
    )

"""
    jacobian_t(s, dom, x, dx_t, x_tt, y)

Jacobian w.r.t. velocity: `∂R/∂x_t · dx_t = damping(dx_t, y)`.
"""
jacobian_t(s, dom::IntegrationDomains, x, dx_t, x_tt, y) =
    _require_nonempty(
        has_damping_form(s) ? damping(s, dom, dx_t, y) : nothing,
        "jacobian_t", s
    )

# Two-entity coupling variant
jacobian_t(a, b, dom::IntegrationDomains, x, dx_t, x_tt, y) =
    _require_nonempty(
        has_damping_form(a, b) ? damping(a, b, dom, dx_t, y) : nothing,
        "jacobian_t", (a, b)
    )

"""
    jacobian_tt(s, dom, x, x_t, dx_tt, y)

Jacobian w.r.t. acceleration: `∂R/∂x_tt · dx_tt = mass(dx_tt, y)`.
"""
jacobian_tt(s, dom::IntegrationDomains, x, x_t, dx_tt, y) =
    _require_nonempty(
        has_mass_form(s) ? mass(s, dom, dx_tt, y) : nothing,
        "jacobian_tt", s
    )

# Two-entity coupling variant
jacobian_tt(a, b, dom::IntegrationDomains, x, x_t, dx_tt, y) =
    _require_nonempty(
        has_mass_form(a, b) ? mass(a, b, dom, dx_tt, y) : nothing,
        "jacobian_tt", (a, b)
    )

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export PhysicsParameters, print_parameters
export mass, damping, stiffness, rhs
export weak_form, residual, jacobian, jacobian_t, jacobian_tt
export AbstractStructure, PotentialFlow, FreeSurface, Membrane2D, EulerBernoulliBeam
export ResonatorSingle, resonator_array

end # module Entities
