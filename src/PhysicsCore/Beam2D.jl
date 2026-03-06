"""
    Beam2D <: AbstractStructure

Parameters for a 2D Euler-Bernoulli beam model (no joints).

# Fields
- `L::Float64` — Length of beam
- `m::Float64` — Mass per unit length per unit width
- `E::Float64` — Young's Modulus
- `I::Float64` — Second Moment of Area
- `τ::Float64` — Stiffness Proportional Structural Damping coefficient
- `bndType::BoundaryCondition` — Boundary Type
- `EI::Float64` — Flexural Rigidity (derived: `E * I`)
- `τEI::Float64` — Damping Rigidity (derived: `τ * EI`)
- `MTotal::Float64` — Total Mass per unit width (derived: `m * L`)
- `ωn1::Float64` — Dry Analytical Natural frequency (derived)
"""
@with_kw struct Beam2D <: AbstractStructure
    L::Float64
    m::Float64
    E::Float64
    I::Float64
    τ::Float64     = 0.0
    bndType::BoundaryCondition = FreeBoundary()

    # Derived quantities
    EI::Float64    = E * I
    τEI::Float64   = τ * EI
    MTotal::Float64 = m * L
    ωn1::Float64   = 22.3733 * sqrt(EI / (m * L^4))
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
