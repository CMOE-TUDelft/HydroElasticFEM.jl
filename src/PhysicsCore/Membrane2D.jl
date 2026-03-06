"""
    Membrane2D <: AbstractStructure

Parameters for a 2D membrane model.

# Fields
- `L::Float64` — Length of membrane
- `m::Float64` — Mass per unit length per unit width
- `T::Float64` — Pre-Tension
- `τ::Float64` — Proportional Structural Damping coefficient
- `bndType::BoundaryCondition` — Boundary Type
- `MTotal::Float64` — Total Mass per unit width (derived: `m * L`)
- `ωn1::Float64` — Dry Analytical Natural frequency (derived: `(π/L) * √(T/m)`)
"""
@with_kw struct Membrane2D <: AbstractStructure
    L::Float64
    m::Float64
    T::Float64
    τ::Float64     = 0.0
    bndType::BoundaryCondition = FreeBoundary()

    # Derived quantities
    MTotal::Float64 = m * L
    ωn1::Float64    = (π / L) * sqrt(T / m)
end

function print_parameters(memb::Membrane2D, ρw::Real=1025)
    mᵨ = memb.m / ρw
    Tᵨ = memb.T / ρw
    @printf("\n[MSG] Membrane Properties:\n")
    @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", ρw)
    @printf("[VAL] Lm = %.4f m\n", memb.L)
    @printf("[VAL] m, mᵨ = %.4f kg/m2, %.4f m\n", memb.m, mᵨ)
    @printf("[VAL] T, Tᵨ = %.4f N/m, %.4f m3/s2\n", memb.T, Tᵨ)
    @printf("[VAL] τ = %.4f \n", memb.τ)
    @printf("[VAL] memBndType = %s \n", string(memb.bndType))
    @printf("[VAL] MTotal = %.4f kg/m \n", memb.MTotal)
    @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s \n", memb.ωn1)
    println()
end
