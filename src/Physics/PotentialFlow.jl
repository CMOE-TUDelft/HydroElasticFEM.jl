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

Wave-damping boundary condition on a top-surface damping zone.

Users provide the incident free-surface elevation `η_in` and vertical velocity
`vz_in`, while the package derives the weak-form RHS terms internally via

- `ηd = μ₂ * η_in`
- `∇ₙϕd = μ₁ * vz_in`
"""
@with_kw struct DampingZoneBC <: AbstractPotentialFlowBC
    domain::Symbol
    μ₁::Any = 0.0
    μ₂::Any = 0.0
    η_in::Any = 0.0
    vz_in::Any = 0.0
    enabled::Bool = true
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

function stiffness(pf::PotentialFlow, ctx::AC.AbstractAssemblyContext, x, y)
    dom = AC.domains(ctx)
    sym = variable_symbol(pf)
    ϕ = x[sym]
    w = y[sym]
    val = ∫(∇(w) ⋅ ∇(ϕ))dom[:dΩ]
    bc_val = _stiffness_bc_contributions(pf, ctx, ϕ, w)
    return _add_contribution(val, bc_val)
end

function rhs(pf::PotentialFlow, dom::IntegrationDomains, f, y)
    sym = variable_symbol(pf)
    w = y[sym]
    val = ∫(w * f[sym])dom[:dΓin]
    bc_val = _rhs_bc_contributions(pf, dom, w)
    return _add_contribution(val, bc_val)
end

function rhs(pf::PotentialFlow, ctx::AC.AbstractAssemblyContext, f, y)
    dom = AC.domains(ctx)
    sym = variable_symbol(pf)
    w = y[sym]
    val = ∫(w * f[sym])dom[:dΓin]
    bc_val = _rhs_bc_contributions(pf, ctx, w)
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

function _stiffness_bc_contributions(pf::PotentialFlow, ctx::AC.AbstractAssemblyContext, ϕ, w)
    val = nothing
    for bc in pf.boundary_conditions
        val = _add_contribution(val, _stiffness_bc_contribution(pf, bc, ctx, ϕ, w))
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
_stiffness_bc_contribution(::PotentialFlow, ::AbstractPotentialFlowBC, ::AC.AbstractAssemblyContext, ϕ, w) = nothing

function _stiffness_bc_contribution(pf::PotentialFlow, bc::RadiationBC, dom::IntegrationDomains, ϕ, w)
    bc.enabled || return nothing
    _ = _radiation_frequency(pf)
    k = _radiation_wavenumber(pf)
    dΓ = dom[bc.domain]
    return ∫(-im * k * w * ϕ)dΓ
end

function _stiffness_bc_contribution(::PotentialFlow, bc::DampingZoneBC, dom::IntegrationDomains, ϕ, w)
    bc.enabled || return nothing
    dΓ = dom[bc.domain]
    nΓ = _boundary_normal(dom, bc.domain)
    μ₁ = _as_space_function(bc.μ₁)
    αₕ = _stabilization_parameter(dom)
    ∇ₙϕ = ∇(ϕ) ⋅ nΓ
    μ₁αₕ(x) = μ₁(x) * αₕ
    return ∫(μ₁αₕ * w * ∇ₙϕ)dΓ
end

function _stiffness_bc_contribution(::PotentialFlow, bc::DampingZoneBC, ctx::AC.AbstractAssemblyContext, ϕ, w)
    bc.enabled || return nothing
    dom = AC.domains(ctx)
    dΓ = dom[bc.domain]
    nΓ = _boundary_normal(dom, bc.domain)
    μ₁ = _as_space_function(bc.μ₁)
    αₕ = _stabilization_parameter(ctx)
    ∇ₙϕ = ∇(ϕ) ⋅ nΓ
    μ₁αₕ(x) = μ₁(x) * αₕ
    return ∫(μ₁αₕ * w * ∇ₙϕ)dΓ
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

function _rhs_bc_contributions(pf::PotentialFlow, ctx::AC.AbstractAssemblyContext, w)
    val = nothing
    for bc in pf.boundary_conditions
        val = _add_contribution(val, _rhs_bc_contribution(pf, bc, ctx, w))
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
_rhs_bc_contribution(::PotentialFlow, ::AbstractPotentialFlowBC, ::AC.AbstractAssemblyContext, w) = nothing

function _rhs_bc_contribution(pf::PotentialFlow, bc::PrescribedInletPotentialBC, dom::IntegrationDomains, w)
    forcing = _resolve_space_function(bc.forcing, dom)
    return _prescribed_rhs_contribution(Val(bc.quantity), pf, forcing, dom[bc.domain], w)
end

function _rhs_bc_contribution(pf::PotentialFlow, bc::PrescribedInletPotentialBC, ctx::AC.AbstractAssemblyContext, w)
    dom = AC.domains(ctx)
    forcing = _resolve_space_function(bc.forcing, ctx)
    return _prescribed_rhs_contribution(Val(bc.quantity), pf, forcing, dom[bc.domain], w)
end

function _rhs_bc_contribution(::PotentialFlow, bc::DampingZoneBC, dom::IntegrationDomains, w)
    bc.enabled || return nothing
    dΓ = dom[bc.domain]
    αₕ = _stabilization_parameter(dom)
    ηd = _ηd(bc, dom)
    ∇ₙϕd = _normal_damped(bc, dom)
    f(x) = -ηd(x) + αₕ * ∇ₙϕd(x)
    return ∫(f * w)dΓ
end

function _rhs_bc_contribution(::PotentialFlow, bc::DampingZoneBC, ctx::AC.AbstractAssemblyContext, w)
    bc.enabled || return nothing
    dom = AC.domains(ctx)
    dΓ = dom[bc.domain]
    αₕ = _stabilization_parameter(ctx)
    ηd = _ηd(bc, ctx)
    ∇ₙϕd = _normal_damped(bc, ctx)
    f(x) = -ηd(x) + αₕ * ∇ₙϕd(x)
    return ∫(f * w)dΓ
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

_active_damping_zone_bcs(pf::PotentialFlow) = [bc for bc in pf.boundary_conditions if bc isa DampingZoneBC && bc.enabled]

function _damping_zone_enabled(pf::PotentialFlow)
    !isempty(_active_damping_zone_bcs(pf))
end

function _boundary_normal(dom::IntegrationDomains, domain::Symbol)
    key = Symbol(replace(String(domain), "dΓ" => "nΓ", count=1))
    haskey(dom, key) || error("Normal key `$key` was not found in IntegrationDomains for damping zone `$domain`.")
    return dom[key]
end

function _stabilization_parameter(ctx::AC.AbstractAssemblyContext)
    AC.has_stabilization(ctx) || error("Damping-zone stabilization requires `αₕ` in the assembly context.")
    return AC.stabilization_parameter(ctx)
end

_stabilization_parameter(::IntegrationDomains) =
    error("Damping-zone stabilization now requires an immutable assembly context.")

function _ηd(bc::DampingZoneBC, dom)
    μ₂ = _as_space_function(bc.μ₂)
    η_in = _resolve_space_function(bc.η_in, dom)
    x -> μ₂(x) * η_in(x)
end

function _normal_damped(bc::DampingZoneBC, dom)
    μ₁ = _as_space_function(bc.μ₁)
    vz_in = _resolve_space_function(bc.vz_in, dom)
    x -> μ₁(x) * vz_in(x)
end

function _ηd(bc::DampingZoneBC, ctx::AC.AbstractAssemblyContext)
    μ₂ = _as_space_function(bc.μ₂)
    η_in = _resolve_space_function(bc.η_in, ctx)
    x -> μ₂(x) * η_in(x)
end

function _normal_damped(bc::DampingZoneBC, ctx::AC.AbstractAssemblyContext)
    μ₁ = _as_space_function(bc.μ₁)
    vz_in = _resolve_space_function(bc.vz_in, ctx)
    x -> μ₁(x) * vz_in(x)
end
