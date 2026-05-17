"""
    FreeSurface <: PhysicsParameters

Linearised free-surface boundary condition parameters.

Introduces an auxiliary surface-elevation variable `κ` [m] that lives on the
free-surface boundary `Γfs` (the portion of the top boundary outside any
structure).  The standalone stiffness form captures the gravity restoring
term `βₕ·g·u·κ`; coupling with the velocity potential `ϕ` (from
[`PotentialFlow`](@ref)) is handled in `CouplingTerms.jl`.

The time-stepping stabilisation coefficient `αₕ` is derived from `βₕ` and
the angular frequency `ω` in `build_frequency_context` / `build_time_context`.

# Fields
- `ρw::Float64` — Fluid density [kg/m³]; default 1025.0
- `g::Float64`  — Gravitational acceleration [m/s²]; default 9.81
- `βₕ::Float64` — Free-surface stabilisation parameter ∈ (0, 1]; default 0.5.
  Setting `βₕ = 1` recovers the fully explicit free-surface condition;
  `βₕ = 0.5` gives a balanced split between explicit and implicit parts.
- `fe::FESpaceConfig` — FE discretisation parameters
- `dim::Int`    — Ambient spatial dimension; default 2
- `symbol::Symbol` — Field unknown symbol; default `:κ`
- `space_domain_symbol::Symbol` — Triangulation key for FE spaces; default `:Γκ`

See also: [`PotentialFlow`](@ref)

# Reference
- [C23] Colomes, O., Verdugo, F., & Akkerman, I. (2023). A monolithic
    finite element formulation for the hydroelastic analysis of very large
    floating structures. *Int. J. Numer. Methods Eng.*, 124(3), 714-751.
    DOI: https://doi.org/10.1002/nme.7140
"""
@with_kw struct FreeSurface <: PhysicsParameters
    ρw::Float64 = 1025.0
    g::Float64  = 9.81
    βₕ::Float64 = 0.5
    fe::FESpaceConfig = FESpaceConfig()
    dim::Int    = 2
    symbol::Symbol = :κ
    space_domain_symbol::Symbol = :Γκ
end

function print_parameters(fs::FreeSurface)
    @printf("\n[MSG] Free Surface Properties:\n")
    @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", fs.ρw)
    @printf("[VAL] Gravitational acceleration, g = %.4f m/s2\n", fs.g)
    @printf("[VAL] Stabilisation parameter, βₕ = %.4f\n", fs.βₕ)
    println()
end

variable_symbol(s::FreeSurface) = s.symbol
ambient_dimension(fs::FreeSurface) = fs.dim

function stabilization_parameter(fs::FreeSurface, ctx::AC.FrequencyAssemblyContext)
    AC.has_stabilization(ctx) || error("Frequency-domain free-surface stabilization requires `αₕ` in the assembly context.")
    return AC.stabilization_parameter(ctx)
end

function stabilization_parameter(fs::FreeSurface, ctx::AC.TimeAssemblyContext)
    AC.has_stabilization(ctx) || error("Time-domain damping-zone stabilization requires `TimeConfig.αₕ`.")
    return AC.stabilization_parameter(ctx)
end

# κ has no standalone mass, damping or rhs; stiffness captures the
# gravity restoring term βₕ·g·u·κ on Γκ.
has_mass_form(::FreeSurface) = false
has_damping_form(::FreeSurface) = false
has_stiffness_form(::FreeSurface) = true
has_rhs_form(::FreeSurface) = false

function stiffness(fs::FreeSurface, dom::IntegrationDomains, x, y)
    κ_sym = variable_symbol(fs)
    κ = x[κ_sym]
    u = y[κ_sym]
    βₕ = fs.βₕ
    g  = fs.g
    # Free-surface gravity restoring term.
    # Standalone weak form: ∫_Γκ βₕ·g·u·κ dΓ.
    # Full stabilized PF-FS coupling is in CouplingTerms.jl.
    # Reference: [C23] Section 2.2, Eq. (7)-(10).
    ∫(βₕ * g * u * κ)dom[:dΓκ]
end
