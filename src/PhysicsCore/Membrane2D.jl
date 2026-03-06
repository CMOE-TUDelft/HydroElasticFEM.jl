"""
    Membrane2D <: AbstractStructure

Parameters for a 2D membrane model.

# Fields
- `L::Real` — Length of membrane
- `m::Real` — Mass per unit length per unit width
- `T::Real` — Pre-Tension
- `τ::Real` — Proportional Structural Damping coefficient
- `bndType::BoundaryCondition` — Boundary Type
- `MTotal::Real` — Total Mass per unit width (derived)
- `ωn1::Real` — Dry Analytical Natural frequency (derived)
"""
struct Membrane2D <: AbstractStructure
    L::Real
    m::Real
    T::Real
    τ::Real
    bndType::BoundaryCondition

    # Derived quantities
    MTotal::Real
    ωn1::Real

    function Membrane2D(L, m, T, τ, bndType::BoundaryCondition)
        MTotal = m * L
        ωn1 = (π / L) * sqrt(T / m)
        new(L, m, T, τ, bndType, MTotal, ωn1)
    end
end

function Membrane2D(bndType::BoundaryCondition=FreeBoundary())
    Membrane2D(0.0, 0.0, 0.0, 0.0, bndType)
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
