"""
    Beam2D <: AbstractStructure

Parameters for a 2D Euler-Bernoulli beam model (no joints).

# Fields
- `L::Real` — Length of beam
- `m::Real` — Mass per unit length per unit width
- `E::Real` — Young's Modulus
- `I::Real` — Second Moment of Area
- `τ::Real` — Stiffness Proportional Structural Damping coefficient
- `bndType::BoundaryCondition` — Boundary Type
- `EI::Real` — Flexural Rigidity (derived)
- `τEI::Real` — Damping Rigidity (derived)
- `MTotal::Real` — Total Mass per unit width (derived)
- `ωn1::Real` — Dry Analytical Natural frequency (derived)
"""
struct Beam2D <: AbstractStructure
    L::Real
    m::Real
    E::Real
    I::Real
    τ::Real
    bndType::BoundaryCondition

    # Derived quantities
    EI::Real
    τEI::Real
    MTotal::Real
    ωn1::Real

    function Beam2D(L, m, E, I, τ, bndType::BoundaryCondition)
        EI = E * I
        τEI = τ * EI
        MTotal = m * L
        ωn1 = 22.3733 * sqrt(EI / (m * L^4))
        new(L, m, E, I, τ, bndType, EI, τEI, MTotal, ωn1)
    end
end

function Beam2D(bndType::BoundaryCondition=FreeBoundary())
    Beam2D(0.0, 0.0, 0.0, 0.0, 0.0, bndType)
end

function print_parameters(beam::Beam2D, ρw::Real=1025)
    mᵨ = beam.m / ρw
    EIᵨ = beam.EI / ρw
    @printf("\n[MSG] Beam Properties:\n")
    @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", ρw)
    @printf("[VAL] L = %.4f m\n", beam.L)
    @printf("[VAL] m, mᵨ = %.4f kg/m2, %.4f m\n", beam.m, mᵨ)
    @printf("[VAL] E = %.4f Pa\n", beam.E)
    @printf("[VAL] I = %.4f m4/m\n", beam.I)
    @printf("[VAL] τ = %.4f \n", beam.τ)
    @printf("[VAL] EI, EIᵨ = %.4f Nm2/m, %.4f m5/s2\n", beam.EI, EIᵨ)
    @printf("[VAL] τEI = %.4f \n", beam.τEI)
    @printf("[VAL] beamBndType = %s \n", string(beam.bndType))
    @printf("[VAL] MTotal = %.4f kg \n", beam.MTotal)
    @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s\n", beam.ωn1)
    @printf("[MSG] See free-free vibration frequency formula involving 22.3733 from Wiki.\n")
    println()
end
