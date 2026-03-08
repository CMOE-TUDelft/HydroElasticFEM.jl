"""
    Membrane2D <: AbstractStructure

Parameters for a 2D membrane model.

# Fields
- `L::Float64` — Length of membrane
- `m::Float64` — Mass per unit length per unit width
- `T::Float64` — Pre-Tension
- `τ::Float64` — Proportional Structural Damping coefficient
- `ρw::Float64` — Density of water
- `g::Float64` — Gravitational acceleration
- `bndType::BoundaryCondition` — Boundary Type
- `MTotal::Float64` — Total Mass per unit width (derived: `m * L`)
- `ωn1::Float64` — Dry Analytical Natural frequency (derived: `(π/L) * √(T/m)`)
"""
@with_kw struct Membrane2D <: AbstractStructure
    L::Float64
    m::Float64
    T::Float64
    τ::Float64     = 0.0
    ρw::Float64    = 1025.0
    g::Float64     = 9.81
    bndType::BoundaryCondition = FreeBoundary()

    # Derived quantities
    MTotal::Float64 = m * L
    ωn1::Float64    = (π / L) * sqrt(T / m)
end

function print_parameters(memb::Membrane2D)
    mᵨ = memb.m / memb.ρw
    Tᵨ = memb.T / memb.ρw
    @printf("\n[MSG] Membrane Properties:\n")
    @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", memb.ρw)
    @printf("[VAL] Lm = %.4f m\n", memb.L)
    @printf("[VAL] m, mᵨ = %.4f kg/m2, %.4f m\n", memb.m, mᵨ)
    @printf("[VAL] T, Tᵨ = %.4f N/m, %.4f m3/s2\n", memb.T, Tᵨ)
    @printf("[VAL] τ = %.4f \n", memb.τ)
    @printf("[VAL] memBndType = %s \n", string(memb.bndType))
    @printf("[VAL] MTotal = %.4f kg/m \n", memb.MTotal)
    @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s \n", memb.ωn1)
    println()
end

η_symbol(::Membrane2D) = :η_m

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only η_m terms — no coupling to ϕ or other fields

function mass(s::Membrane2D, dom::WeakFormDomains, x_tt, y)
    ηₜₜ = x_tt[:η_m]
    v   = y[:η_m]
    ∫((s.m / s.ρw) * v * ηₜₜ)dom[:dΓ_s]
end

function damping(s::Membrane2D, dom::WeakFormDomains, x_t, y)
    ηₜ = x_t[:η_m]
    v  = y[:η_m]
    Tᵨ = s.T / s.ρw
    τ  = s.τ
    val = ∫(Tᵨ * τ * ∇(v) ⋅ ∇(ηₜ))dom[:dΓ_s]
    if s.bndType isa FixedBoundary
        val += ∫(-Tᵨ * τ * v * ∇(ηₜ) ⋅ dom[:n_Λ_sb])dom[:dΛ_sb]
    end
    return val
end

function stiffness(s::Membrane2D, dom::WeakFormDomains, x, y)
    η = x[:η_m]
    v = y[:η_m]
    Tᵨ = s.T / s.ρw
    val = ∫(v * (s.g * η) + Tᵨ * ∇(v) ⋅ ∇(η))dom[:dΓ_s]
    if s.bndType isa FixedBoundary
        val += ∫(-Tᵨ * v * ∇(η) ⋅ dom[:n_Λ_sb])dom[:dΛ_sb]
    end
    return val
end

function rhs(s::Membrane2D, dom::WeakFormDomains, f, y)
    v = y[:η_m]
    ∫(v * f[:η_m])dom[:dΓ_s]
end
