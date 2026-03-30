"""
    AbstractPotentialFlowBC

Abstract type for potential flow boundary conditions.
"""
abstract type AbstractPotentialFlowBC end

"""
    RadiationBC <: AbstractPotentialFlowBC

Radiation boundary condition for potential flow, based on the linearized free-surface 
condition and the dispersion relation. Can be applied to any boundary, but typically 
used on open boundaries (e.g., `:dΓout`) to allow outgoing waves to radiate without 
reflection.
"""
@with_kw struct RadiationBC <: AbstractPotentialFlowBC
    domain::Symbol = :dΓout
    enabled::Bool = true
end

"""
    PrescribedInletPotentialBC <: AbstractPotentialFlowBC

Prescribed potential boundary condition for potential flow, typically used on inlet
boundaries (e.g., `:dΓin`) to specify incoming wave conditions. The `forcing` field can
be a constant value or a function of space (and time, if needed) that defines the
potential on the specified boundary. The `quantity` field indicates whether the forcing
represents a potential value, a normal gradient (Neumann condition), or a traction
condition, which affects how the contribution is added to the weak form.
"""
@with_kw struct PrescribedInletPotentialBC <: AbstractPotentialFlowBC
    domain::Symbol = :dΓin
    forcing::Any
    quantity::Symbol = :potential
end

"""
    DampingZoneBC <: AbstractPotentialFlowBC
Damping zone boundary condition for potential flow, used to model wave absorption in
specific regions (e.g., near boundaries) to minimize reflections. The `σ` field can be a 
constant damping coefficient or a function that defines spatially varying damping. 
The `rhs_forcing` field allows for additional forcing terms to be included in the weak 
form, which can be useful for modeling active wave absorption or other effects within 
the damping zone.
"""
@with_kw struct DampingZoneBC <: AbstractPotentialFlowBC
    domain::Symbol
    σ::Any = 0.0
    rhs_forcing::Any = nothing
end

# ─────────────────────────────────────────────────────────────
# PotentialFlow physics parameters and weak form implementations
# ─────────────────────────────────────────────────────────────

"""
    PotentialFlow <: PhysicsParameters

Parameters for the 2D fluid potential (Laplace equation).

# Fields
- `ρw::Float64` — Density of water
- `g::Float64`  — Gravitational acceleration
"""
@with_kw struct PotentialFlow <: PhysicsParameters
    ρw::Float64 = 1025.0
    g::Float64  = 9.81
    fe::FESpaceConfig = FESpaceConfig()
    symbol::Symbol = :ϕ
    space_domain_symbol::Symbol = :Ω
    sea_state::Union{WaveSpec.AiryWaves.AiryState, Nothing} = nothing
    boundary_conditions::Vector{AbstractPotentialFlowBC} = AbstractPotentialFlowBC[]
end

function print_parameters(f::PotentialFlow)
    @printf("\n[MSG] Fluid Properties:\n")
    @printf("[VAL] ρw = %.2f kg/m3\n", f.ρw)
    @printf("[VAL] g  = %.4f m/s2\n", f.g)
    println()
end

variable_symbol(s::PotentialFlow) = s.symbol


# Potential-flow dynamics are elliptic: no standalone mass/damping forms.
has_mass_form(::PotentialFlow) = false
has_damping_form(::PotentialFlow) = false

# ── Single-variable weak forms ─────────────────────────────
#    Field access via variable_symbol (velocity potential)

function stiffness(pf::PotentialFlow, dom::IntegrationDomains, x, y)
    sym = variable_symbol(pf)
    ϕ = x[sym]
    w = y[sym]
    val = ∫(∇(w) ⋅ ∇(ϕ))dom[:dΩ]
    bc_val = _stiffness_bc_contributions(pf, dom, ϕ, w)
    return _add_contribution(val, bc_val)
end

function rhs(pf::PotentialFlow, dom::IntegrationDomains, f, y)
    sym = variable_symbol(pf)
    w = y[sym]
    val = ∫(w * f[sym])dom[:dΓin]
    bc_val = _rhs_bc_contributions(pf, dom, w)
    return _add_contribution(val, bc_val)
end

# ----────────────────────────────────────────────────────────────
# Stiffness contributions from boundary conditions (e.g., radiation BC)
# -----────────────────────────────────────────────────────────────

"""
    _stiffness_bc_contributions(pf::PotentialFlow, dom::IntegrationDomains, ϕ, w)

Compute the total stiffness contributions from all boundary conditions defined in the `PotentialFlow` physics
instance `pf`. This function iterates over all boundary conditions, checks their types, and sums their 
contributions to the stiffness form.
"""
function _stiffness_bc_contributions(pf::PotentialFlow, dom::IntegrationDomains, ϕ, w)
    val = nothing
    for bc in pf.boundary_conditions
        val = _add_contribution(val, _stiffness_bc_contribution(pf, bc, dom, ϕ, w))
    end
    return val
end


"""
    _stiffness_bc_contribution(pf::PotentialFlow, bc::AbstractPotentialFlowBC, dom::IntegrationDomains, ϕ, w)

Compute the stiffness contribution from a given boundary condition `bc` for the potential flow physics `pf`. 
This function dispatches to specific implementations based on the type of boundary condition, allowing for 
flexible handling of different BC types (e.g., radiation, damping zones). The contributions are integrated 
over the relevant domains specified in the BC and added to the overall stiffness form.
"""
_stiffness_bc_contribution(::PotentialFlow, ::AbstractPotentialFlowBC, ::IntegrationDomains, ϕ, w) = nothing

function _stiffness_bc_contribution(pf::PotentialFlow, bc::RadiationBC, dom::IntegrationDomains, ϕ, w)
    bc.enabled || return nothing
    _ = _radiation_frequency(pf)
    k = _radiation_wavenumber(pf)
    dΓ = dom[bc.domain]
    return ∫(-im * k * w * ϕ)dΓ
end

function _stiffness_bc_contribution(::PotentialFlow, bc::DampingZoneBC, dom::IntegrationDomains, ϕ, w)
    σ = _as_space_function(bc.σ)
    dΓ = dom[bc.domain]
    return ∫(-σ * w * ϕ)dΓ
end

# -----────────────────────────────────────────────────────────────
# RHS contributions from boundary conditions (e.g., prescribed inlet potential)
# -----────────────────────────────────────────────────────────────


"""
    _rhs_bc_contributions(pf::PotentialFlow, dom::IntegrationDomains, w)

Compute the total right-hand side contributions from all boundary conditions defined in the `PotentialFlow` physics
instance `pf`. This function iterates over all boundary conditions, checks their types, and sums their contributions
to the rhs form.
"""
function _rhs_bc_contributions(pf::PotentialFlow, dom::IntegrationDomains, w)
    val = nothing
    for bc in pf.boundary_conditions
        val = _add_contribution(val, _rhs_bc_contribution(pf, bc, dom, w))
    end
    return val
end

"""
    _rhs_bc_contribution(pf::PotentialFlow, bc::AbstractPotentialFlowBC, dom::IntegrationDomains, w)

Compute the right-hand side contribution from a given boundary condition `bc` for the potential flow physics `pf`.
This function dispatches to specific implementations based on the type of boundary condition, allowing for
flexible handling of different BC types (e.g., prescribed inlet potential, damping zones). The contributions are integrated
over the relevant domains specified in the BC and added to the overall rhs form.
"""
_rhs_bc_contribution(::PotentialFlow, ::AbstractPotentialFlowBC, ::IntegrationDomains, w) = nothing

function _rhs_bc_contribution(pf::PotentialFlow, bc::PrescribedInletPotentialBC, dom::IntegrationDomains, w)
    forcing = _as_space_function(bc.forcing)
    return _prescribed_rhs_contribution(Val(bc.quantity), pf, forcing, dom[bc.domain], w)
end

function _rhs_bc_contribution(::PotentialFlow, bc::DampingZoneBC, dom::IntegrationDomains, w)
    isnothing(bc.rhs_forcing) && return nothing
    forcing = _as_space_function(bc.rhs_forcing)
    return ∫(w * forcing)dom[bc.domain]
end

_prescribed_rhs_contribution(::Val{:traction}, ::PotentialFlow, forcing, dΓ, w) = ∫(w * forcing)dΓ
_prescribed_rhs_contribution(::Val{:normal_gradient}, ::PotentialFlow, forcing, dΓ, w) = ∫(w * forcing)dΓ
function _prescribed_rhs_contribution(::Val{:potential}, pf::PotentialFlow, forcing, dΓ, w)
    k = _radiation_wavenumber(pf)
    return ∫(-im * k * w * forcing)dΓ
end
function _prescribed_rhs_contribution(::Val{Q}, ::PotentialFlow, forcing, dΓ, w) where {Q}
    error("Unsupported PrescribedInletPotentialBC quantity `$(Q)`. Expected one of :potential, :normal_gradient, or :traction.")
end

# -----────────────────────────────────────────────────────────────
# Helper functions for BC contributions
# -----────────────────────────────────────────────────────────────

function _single_frequency_wave(pf::PotentialFlow)
    isnothing(pf.sea_state) && error("PotentialFlow radiation requires `sea_state` to be defined.")
    length(pf.sea_state.ω) == 1 || error("PotentialFlow radiation requires a single-frequency `sea_state` (length(sea_state.ω) == 1).")
    length(pf.sea_state.k) == 1 || error("PotentialFlow radiation requires a single wavenumber entry (length(sea_state.k) == 1).")
    return pf.sea_state
end

_radiation_wavenumber(pf::PotentialFlow) = _single_frequency_wave(pf).k[1]
_radiation_frequency(pf::PotentialFlow) = _single_frequency_wave(pf).ω[1]
