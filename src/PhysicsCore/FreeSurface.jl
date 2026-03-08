"""
    FreeSurface <: PhysicsParameters

Linearised free-surface boundary condition parameters.

Introduces an auxiliary surface-elevation variable `κ` that lives on the
free-surface boundary `Γ_fs` (the portion of the top boundary outside any
structure).  The coupling with the velocity potential `ϕ` (from
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
end

function print_parameters(fs::FreeSurface)
    @printf("\n[MSG] Free Surface Properties:\n")
    @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", fs.ρw)
    @printf("[VAL] Gravitational acceleration, g = %.4f m/s2\n", fs.g)
    @printf("[VAL] Stabilisation parameter, βₕ = %.4f\n", fs.βₕ)
    println()
end

variable_symbol(::FreeSurface) = :κ

# ── Single-variable weak forms ──
#    κ has no standalone dynamics; all physics enters via
#    coupling with PotentialFlow (see CouplingTerms.jl).

function mass(fs::FreeSurface, dom::WeakFormDomains, x_tt, y)
    κ_sym = variable_symbol(fs)
    u = y[κ_sym]
    ∫(0.0 * u)dom[:dΓ_fs]
end

function damping(fs::FreeSurface, dom::WeakFormDomains, x_t, y)
    κ_sym = variable_symbol(fs)
    u = y[κ_sym]
    ∫(0.0 * u)dom[:dΓ_fs]
end

function stiffness(fs::FreeSurface, dom::WeakFormDomains, x, y)
    κ_sym = variable_symbol(fs)
    u = y[κ_sym]
    ∫(0.0 * u)dom[:dΓ_fs]
end

function rhs(fs::FreeSurface, dom::WeakFormDomains, f, y)
    κ_sym = variable_symbol(fs)
    u = y[κ_sym]
    ∫(0.0 * u)dom[:dΓ_fs]
end
