"""
    module PhysicalEntities

Unified type hierarchy for all physical entities in HydroElasticFEM.

Each entity file defines `mass`, `damping`, `stiffness`, and `rhs`
methods that access FE fields by symbol (e.g. `x[:ϕ]`, `x[:η_m]`)
via a `FieldDict` wrapper provided by `WeakFormAssembly`.

Generic composed forms (`weakform`, `residual`, `jacobian`, ...)
are defined at module level and dispatch to the linear forms.
"""
module PhysicalEntities

using Parameters
using Gridap
using Printf

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

# ─────────────────────────────────────────────────────────────
# Boundary conditions
# ─────────────────────────────────────────────────────────────

abstract type BoundaryCondition end
struct FreeBoundary  <: BoundaryCondition end
struct FixedBoundary <: BoundaryCondition end

# ─────────────────────────────────────────────────────────────
# Structural entities
# ─────────────────────────────────────────────────────────────

abstract type AbstractStructure <: PhysicsParameters end

# Domain container (must be included before entity files that reference it)
include("WeakFormDomains.jl")

# Entity files (struct definition + single-variable weak forms)
include("PotentialFlow.jl")
include("Membrane2D.jl")
include("Beam2D.jl")
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
weakform(s, dom::WeakFormDomains, ω, x, y) =
    -ω^2       * mass(s, dom, x, y) +
    (-im * ω)  * damping(s, dom, x, y) +
    stiffness(s, dom, x, y)

# Two-entity coupling variant
weakform(a, b, dom::WeakFormDomains, ω, x, y) =
    -ω^2       * mass(a, b, dom, x, y) +
    (-im * ω)  * damping(a, b, dom, x, y) +
    stiffness(a, b, dom, x, y)

# ─────────────────────────────────────────────────────────────
# Nonlinear weak forms: residual, jacobian, jacobian_t, jacobian_tt
# ─────────────────────────────────────────────────────────────

"""
    residual(s, dom, x, x_t, x_tt, f, y)

Nonlinear residual: `mass(x_tt, y) + damping(x_t, y) + stiffness(x, y) - rhs(f, y)`.
"""
residual(s, dom::WeakFormDomains, x, x_t, x_tt, f, y) =
    mass(s, dom, x_tt, y) +
    damping(s, dom, x_t, y) +
    stiffness(s, dom, x, y) -
    rhs(s, dom, f, y)

# Two-entity coupling variant
residual(a, b, dom::WeakFormDomains, x, x_t, x_tt, f, y) =
    mass(a, b, dom, x_tt, y) +
    damping(a, b, dom, x_t, y) +
    stiffness(a, b, dom, x, y) -
    rhs(a, b, dom, f, y)

"""
    jacobian(s, dom, dx, x_t, x_tt, y)

Jacobian w.r.t. displacement: `∂R/∂x · dx = stiffness(dx, y)`.
"""
jacobian(s, dom::WeakFormDomains, dx, x_t, x_tt, y) =
    stiffness(s, dom, dx, y)

# Two-entity coupling variant
jacobian(a, b, dom::WeakFormDomains, dx, x_t, x_tt, y) =
    stiffness(a, b, dom, dx, y)

"""
    jacobian_t(s, dom, x, dx_t, x_tt, y)

Jacobian w.r.t. velocity: `∂R/∂x_t · dx_t = damping(dx_t, y)`.
"""
jacobian_t(s, dom::WeakFormDomains, x, dx_t, x_tt, y) =
    damping(s, dom, dx_t, y)

# Two-entity coupling variant
jacobian_t(a, b, dom::WeakFormDomains, x, dx_t, x_tt, y) =
    damping(a, b, dom, dx_t, y)

"""
    jacobian_tt(s, dom, x, x_t, dx_tt, y)

Jacobian w.r.t. acceleration: `∂R/∂x_tt · dx_tt = mass(dx_tt, y)`.
"""
jacobian_tt(s, dom::WeakFormDomains, x, x_t, dx_tt, y) =
    mass(s, dom, dx_tt, y)

# Two-entity coupling variant
jacobian_tt(a, b, dom::WeakFormDomains, x, x_t, dx_tt, y) =
    mass(a, b, dom, dx_tt, y)

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export PhysicsParameters, print_parameters
export BoundaryCondition, FreeBoundary, FixedBoundary
export AbstractStructure, PotentialFlow, Membrane2D, Beam2D
export ResonatorSingle, resonator_array
export WeakFormDomains
export η_symbol
export weakform, mass, damping, stiffness, rhs
export residual, jacobian, jacobian_t, jacobian_tt

end # module PhysicalEntities
