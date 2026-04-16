"""
     module Physics

Unified type hierarchy for all physical entities in HydroElasticFEM.

Each entity file may define any subset of `mass`, `damping`,
`stiffness`, and `rhs` methods that access FE fields by symbol
(e.g. `x[:ϕ]`, `x[:η_m]`) via a `FieldMap` wrapper provided by
`FEOperators`.

Generic composed forms (`weakform`, `residual`, `jacobian`, ...)
are defined at module level and dispatch to the linear forms.
"""
module Physics

using Parameters
using Printf
using Gridap
using WaveSpec
using ..Geometry
using ..ParameterHandler
import ..AssemblyContexts as AC

# -----────────────────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────────────────

include("Helpers.jl")

# ─────────────────────────────────────────────────────────────
# Abstract base
# ─────────────────────────────────────────────────────────────

"""
    abstract type PhysicsParameters

Abstract base type for all physics parameter structures.
"""
abstract type PhysicsParameters end

abstract type Structure <:PhysicsParameters end

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

mass(s::PhysicsParameters, ctx::AC.AbstractAssemblyContext, x_tt, y) =
    mass(s, AC.domains(ctx), x_tt, y)
damping(s::PhysicsParameters, ctx::AC.AbstractAssemblyContext, x_t, y) =
    damping(s, AC.domains(ctx), x_t, y)
stiffness(s::PhysicsParameters, ctx::AC.AbstractAssemblyContext, x, y) =
    stiffness(s, AC.domains(ctx), x, y)
rhs(s::PhysicsParameters, ctx::AC.AbstractAssemblyContext, f, y) =
    rhs(s, AC.domains(ctx), f, y)

mass(a, b, ctx::AC.AbstractAssemblyContext, x_tt, y) =
    mass(a, b, AC.domains(ctx), x_tt, y)
damping(a, b, ctx::AC.AbstractAssemblyContext, x_t, y) =
    damping(a, b, AC.domains(ctx), x_t, y)
stiffness(a, b, ctx::AC.AbstractAssemblyContext, x, y) =
    stiffness(a, b, AC.domains(ctx), x, y)
rhs(a, b, ctx::AC.AbstractAssemblyContext, f, y) =
    rhs(a, b, AC.domains(ctx), f, y)

# Optional form-presence traits.
# Single-entity forms default to present for backward compatibility.
"""
    has_mass_form(s::PhysicsParameters)
    has_mass_form(a, b)

Return whether mass-form contributions are active for the given entity or
entity pair. Single-entity forms default to `true`; coupling forms default
to `false` unless specialized.
"""
has_mass_form(::PhysicsParameters) = true
"""
    has_damping_form(s::PhysicsParameters)
    has_damping_form(a, b)

Return whether damping-form contributions are active for the given entity or
entity pair. Single-entity forms default to `true`; coupling forms default
to `false` unless specialized.
"""
has_damping_form(::PhysicsParameters) = true
"""
    has_stiffness_form(s::PhysicsParameters)
    has_stiffness_form(a, b)

Return whether stiffness-form contributions are active for the given entity or
entity pair. Single-entity forms default to `true`; coupling forms default
to `false` unless specialized.
"""
has_stiffness_form(::PhysicsParameters) = true
"""
    has_rhs_form(s::PhysicsParameters)
    has_rhs_form(a, b)

Return whether right-hand-side contributions are active for the given entity or
entity pair. Single-entity forms default to `true`; coupling forms default
to `false` unless specialized.
"""
has_rhs_form(::PhysicsParameters) = true

# Two-entity coupling forms default to absent unless enabled.
has_mass_form(a, b) = false
has_damping_form(a, b) = false
has_stiffness_form(a, b) = false
has_rhs_form(a, b) = false

"""
    active_forms(ctx::AC.AbstractAssemblyContext, s::PhysicsParameters)
    active_forms(ctx::AC.AbstractAssemblyContext, a, b)

Return a named tuple with active form flags for `mass`, `damping`,
`stiffness`, and `rhs` in the provided assembly context.
"""
active_forms(::AC.AbstractAssemblyContext, s::PhysicsParameters) = (
    mass=has_mass_form(s),
    damping=has_damping_form(s),
    stiffness=has_stiffness_form(s),
    rhs=has_rhs_form(s),
)

active_forms(::AC.AbstractAssemblyContext, a, b) = (
    mass=has_mass_form(a, b),
    damping=has_damping_form(a, b),
    stiffness=has_stiffness_form(a, b),
    rhs=has_rhs_form(a, b),
)

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
weakform(s, ctx::AC.FrequencyAssemblyContext, x, y) =
    _require_nonempty(
        _sum_present(
            active_forms(ctx, s).mass ? (-AC.frequency(ctx)^2) * mass(s, ctx, x, y) : nothing,
            active_forms(ctx, s).damping ? (-im * AC.frequency(ctx)) * damping(s, ctx, x, y) : nothing,
            active_forms(ctx, s).stiffness ? stiffness(s, ctx, x, y) : nothing,
        ),
        "weakform", s
    )

# Two-entity coupling variant
weakform(a, b, ctx::AC.FrequencyAssemblyContext, x, y) =
    _require_nonempty(
        _sum_present(
            active_forms(ctx, a, b).mass ? (-AC.frequency(ctx)^2) * mass(a, b, ctx, x, y) : nothing,
            active_forms(ctx, a, b).damping ? (-im * AC.frequency(ctx)) * damping(a, b, ctx, x, y) : nothing,
            active_forms(ctx, a, b).stiffness ? stiffness(a, b, ctx, x, y) : nothing,
        ),
        "weakform", (a, b)
    )

function _compat_frequency_context(obj, dom::IntegrationDomains, ω)
    αₕ = nothing
    if obj isa FreeSurface
        αₕ = -im * ω / obj.g * (1.0 - obj.βₕ) / obj.βₕ
    elseif obj isa Tuple && length(obj) == 2 && obj[2] isa FreeSurface
        fs = obj[2]
        αₕ = -im * ω / fs.g * (1.0 - fs.βₕ) / fs.βₕ
    end
    return AC.FrequencyAssemblyContext(dom, ω, αₕ)
end

weakform(s, dom::IntegrationDomains, ω, x, y) = weakform(s, _compat_frequency_context(s, dom, ω), x, y)
weakform(a, b, dom::IntegrationDomains, ω, x, y) = weakform(a, b, _compat_frequency_context((a, b), dom, ω), x, y)

# ─────────────────────────────────────────────────────────────
# Nonlinear weak forms: residual, jacobian, jacobian_t, jacobian_tt
# ─────────────────────────────────────────────────────────────

"""
    residual(s, dom, x, x_t, x_tt, f, y)

Nonlinear residual: `mass(x_tt, y) + damping(x_t, y) + stiffness(x, y) - rhs(f, y)`.
"""
residual(s, ctx::AC.AbstractAssemblyContext, x, x_t, x_tt, f, y) =
    _require_nonempty(
        _sum_present(
            active_forms(ctx, s).mass ? mass(s, ctx, x_tt, y) : nothing,
            active_forms(ctx, s).damping ? damping(s, ctx, x_t, y) : nothing,
            active_forms(ctx, s).stiffness ? stiffness(s, ctx, x, y) : nothing,
            active_forms(ctx, s).rhs ? (-rhs(s, ctx, f, y)) : nothing,
        ),
        "residual", s
    )

# Two-entity coupling variant
residual(a, b, ctx::AC.AbstractAssemblyContext, x, x_t, x_tt, f, y) =
    _require_nonempty(
        _sum_present(
            active_forms(ctx, a, b).mass ? mass(a, b, ctx, x_tt, y) : nothing,
            active_forms(ctx, a, b).damping ? damping(a, b, ctx, x_t, y) : nothing,
            active_forms(ctx, a, b).stiffness ? stiffness(a, b, ctx, x, y) : nothing,
            active_forms(ctx, a, b).rhs ? (-rhs(a, b, ctx, f, y)) : nothing,
        ),
        "residual", (a, b)
    )

residual(s, dom::IntegrationDomains, x, x_t, x_tt, f, y) =
    residual(s, AC.TimeAssemblyContext(dom, 0.0, nothing), x, x_t, x_tt, f, y)
residual(a, b, dom::IntegrationDomains, x, x_t, x_tt, f, y) =
    residual(a, b, AC.TimeAssemblyContext(dom, 0.0, nothing), x, x_t, x_tt, f, y)

"""
    jacobian(s, dom, dx, x_t, x_tt, y)

Jacobian w.r.t. displacement: `∂R/∂x · dx = stiffness(dx, y)`.
"""
jacobian(s, ctx::AC.AbstractAssemblyContext, dx, x_t, x_tt, y) =
    _require_nonempty(
        active_forms(ctx, s).stiffness ? stiffness(s, ctx, dx, y) : nothing,
        "jacobian", s
    )

# Two-entity coupling variant
jacobian(a, b, ctx::AC.AbstractAssemblyContext, dx, x_t, x_tt, y) =
    _require_nonempty(
        active_forms(ctx, a, b).stiffness ? stiffness(a, b, ctx, dx, y) : nothing,
        "jacobian", (a, b)
    )

jacobian(s, dom::IntegrationDomains, dx, x_t, x_tt, y) =
    jacobian(s, AC.TimeAssemblyContext(dom, 0.0, nothing), dx, x_t, x_tt, y)
jacobian(a, b, dom::IntegrationDomains, dx, x_t, x_tt, y) =
    jacobian(a, b, AC.TimeAssemblyContext(dom, 0.0, nothing), dx, x_t, x_tt, y)

"""
    jacobian_t(s, dom, x, dx_t, x_tt, y)

Jacobian w.r.t. velocity: `∂R/∂x_t · dx_t = damping(dx_t, y)`.
"""
jacobian_t(s, ctx::AC.AbstractAssemblyContext, x, dx_t, x_tt, y) =
    _require_nonempty(
        active_forms(ctx, s).damping ? damping(s, ctx, dx_t, y) : nothing,
        "jacobian_t", s
    )

# Two-entity coupling variant
jacobian_t(a, b, ctx::AC.AbstractAssemblyContext, x, dx_t, x_tt, y) =
    _require_nonempty(
        active_forms(ctx, a, b).damping ? damping(a, b, ctx, dx_t, y) : nothing,
        "jacobian_t", (a, b)
    )

jacobian_t(s, dom::IntegrationDomains, x, dx_t, x_tt, y) =
    jacobian_t(s, AC.TimeAssemblyContext(dom, 0.0, nothing), x, dx_t, x_tt, y)
jacobian_t(a, b, dom::IntegrationDomains, x, dx_t, x_tt, y) =
    jacobian_t(a, b, AC.TimeAssemblyContext(dom, 0.0, nothing), x, dx_t, x_tt, y)

"""
    jacobian_tt(s, dom, x, x_t, dx_tt, y)

Jacobian w.r.t. acceleration: `∂R/∂x_tt · dx_tt = mass(dx_tt, y)`.
"""
jacobian_tt(s, ctx::AC.AbstractAssemblyContext, x, x_t, dx_tt, y) =
    _require_nonempty(
        active_forms(ctx, s).mass ? mass(s, ctx, dx_tt, y) : nothing,
        "jacobian_tt", s
    )

# Two-entity coupling variant
jacobian_tt(a, b, ctx::AC.AbstractAssemblyContext, x, x_t, dx_tt, y) =
    _require_nonempty(
        active_forms(ctx, a, b).mass ? mass(a, b, ctx, dx_tt, y) : nothing,
        "jacobian_tt", (a, b)
    )

jacobian_tt(s, dom::IntegrationDomains, x, x_t, dx_tt, y) =
    jacobian_tt(s, AC.TimeAssemblyContext(dom, 0.0, nothing), x, x_t, dx_tt, y)
jacobian_tt(a, b, dom::IntegrationDomains, x, x_t, dx_tt, y) =
    jacobian_tt(a, b, AC.TimeAssemblyContext(dom, 0.0, nothing), x, x_t, dx_tt, y)

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export PhysicsParameters, print_parameters, variable_symbol
export mass, damping, stiffness, rhs
export has_mass_form, has_damping_form, has_stiffness_form, has_rhs_form
export active_forms
export weakform, residual, jacobian, jacobian_t, jacobian_tt
export PotentialFlow, FreeSurface, Membrane2D, EulerBernoulliBeam, Resonator
export AbstractPotentialFlowBC, RadiationBC, PrescribedInletPotentialBC, DampingZoneBC
export CouplingTerms

end # module Physics
