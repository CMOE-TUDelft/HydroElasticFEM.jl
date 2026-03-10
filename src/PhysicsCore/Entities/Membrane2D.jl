"""
    Membrane2D <: AbstractStructure

Parameters for a 2D membrane model, normalised by fluid density ρw.

# Fields
- `L::Float64`   — Length of membrane
- `mᵨ::Float64`  — Mass per unit length per unit width / ρw
- `Tᵨ::Float64`  — Pre-Tension / ρw
- `τ::Float64`   — Proportional Structural Damping coefficient
- `g::Float64`   — Gravitational acceleration
- `bndType::BoundaryCondition` — Boundary Type
- `ωn1::Float64`  — Dry Analytical Natural frequency (derived: `(π/L) * √(Tᵨ/mᵨ)`)
"""
@with_kw struct Membrane2D <: AbstractStructure
    L::Float64
    mᵨ::Float64
    Tᵨ::Float64
    τ::Float64     = 0.0
    g::Float64     = 9.81
    symbol::Symbol = :η_m

    # Derived quantities
    ωn1::Float64    = (π / L) * sqrt(Tᵨ / mᵨ)
    fe::FESpaceConfig = FESpaceConfig()
end

function print_parameters(memb::Membrane2D)
    @printf("\n[MSG] Membrane Properties:\n")
    @printf("[VAL] L = %.4f m\n", memb.L)
    @printf("[VAL] mᵨ = %.4f m\n", memb.mᵨ)
    @printf("[VAL] Tᵨ = %.4f m3/s2\n", memb.Tᵨ)
    @printf("[VAL] τ = %.4f \n", memb.τ)
    @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s \n", memb.ωn1)
    println()
end

variable_symbol(s::Membrane2D) = s.symbol

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only η_m terms — no coupling to ϕ or other fields

function mass(s::Membrane2D, dom::IntegrationDomains, x_tt, y)
    sym = variable_symbol(s)
    ηₜₜ = x_tt[sym]
    v   = y[sym]
    ∫(s.mᵨ * v * ηₜₜ)dom[:dΓ_s]
end

function damping(s::Membrane2D, dom::IntegrationDomains, x_t, y)
    sym = variable_symbol(s)
    ηₜ = x_t[sym]
    v  = y[sym]
    Tᵨ = s.Tᵨ
    τ  = s.τ
    val = ∫(Tᵨ * τ * ∇(v) ⋅ ∇(ηₜ))dom[:dΓ_s]
    return val
end

function stiffness(s::Membrane2D, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    Tᵨ = s.Tᵨ
    val = ∫(v * (s.g * η) + Tᵨ * ∇(v) ⋅ ∇(η))dom[:dΓ_s]
    return val
end

function rhs(s::Membrane2D, dom::IntegrationDomains, f, y)
    sym = variable_symbol(s)
    v = y[sym]
    ∫(v * f[sym])dom[:dΓ_s]
end
