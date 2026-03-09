"""
    FreeSurface <: PhysicsParameters

Linearised free-surface boundary condition parameters.

Introduces an auxiliary surface-elevation variable `κ` that lives on the
free-surface boundary `Γ_fs` (the portion of the top boundary outside any
structure).  The standalone stiffness form captures the gravity restoring
term `βₕ·g·u·κ`; coupling with the velocity potential `ϕ` (from
`PotentialFlow`) is handled in CouplingTerms.jl.

# Fields
- `ρw::Float64` — Density of water
- `g::Float64`  — Gravitational acceleration
- `βₕ::Float64` — Stabilisation parameter for the free-surface formulation
"""
@with_kw struct FreeSurface <: PhysicsParameters
    ρw::Float64 = 1025.0
    g::Float64  = 9.81
    βₕ::Float64 = 0.5
    fe::FESpaceConfig = FESpaceConfig()
end

function print_parameters(fs::FreeSurface)
    @printf("\n[MSG] Free Surface Properties:\n")
    @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", fs.ρw)
    @printf("[VAL] Gravitational acceleration, g = %.4f m/s2\n", fs.g)
    @printf("[VAL] Stabilisation parameter, βₕ = %.4f\n", fs.βₕ)
    println()
end

variable_symbol(::FreeSurface) = :κ

# κ has no standalone mass, damping or rhs; stiffness captures the
# gravity restoring term βₕ·g·u·κ on Γ_fs.
has_mass_form(::FreeSurface) = false
has_damping_form(::FreeSurface) = false
has_stiffness_form(::FreeSurface) = true
has_rhs_form(::FreeSurface) = false

function stiffness(fs::FreeSurface, dom::WeakFormDomains, x, y)
    κ_sym = variable_symbol(fs)
    κ = x[κ_sym]
    u = y[κ_sym]
    βₕ = fs.βₕ
    g  = fs.g
    ∫(βₕ * g * u * κ)dom[:dΓ_fs]
end
