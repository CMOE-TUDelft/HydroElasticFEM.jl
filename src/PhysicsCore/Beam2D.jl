"""
    Beam2D <: AbstractStructure

Parameters for a 2D Euler-Bernoulli beam model (no joints).

# Fields
- `L::Float64` — Length of beam
- `m::Float64` — Mass per unit length per unit width
- `E::Float64` — Young's Modulus
- `I::Float64` — Second Moment of Area
- `τ::Float64` — Stiffness Proportional Structural Damping coefficient
- `ρw::Float64` — Density of water
- `g::Float64` — Gravitational acceleration
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
    ρw::Float64    = 1025.0
    g::Float64     = 9.81
    bndType::BoundaryCondition = FreeBoundary()

    # Derived quantities
    EI::Float64    = E * I
    τEI::Float64   = τ * EI
    MTotal::Float64 = m * L
    ωn1::Float64   = 22.3733 * sqrt(EI / (m * L^4))
end

function print_parameters(beam::Beam2D)
    mᵨ = beam.m / beam.ρw
    EIᵨ = beam.EI / beam.ρw
    @printf("\n[MSG] Beam Properties:\n")
    @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", beam.ρw)
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

variable_symbol(::Beam2D) = :η_b

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only η_b terms — no coupling to ϕ or other fields

function mass(s::Beam2D, dom::WeakFormDomains, x_tt, y)
    sym = variable_symbol(s)
    ηₜₜ = x_tt[sym]
    v   = y[sym]
    ∫((s.m / s.ρw) * v * ηₜₜ)dom[:dΓ_s]
end

function damping(s::Beam2D, dom::WeakFormDomains, x_t, y)
    sym = variable_symbol(s)
    ηₜ = x_t[sym]
    v  = y[sym]
    EIᵨ = s.EI / s.ρw
    τ   = s.τ
    γ   = dom[:γ_s]
    h   = dom[:h_s]
    n_Λ = dom[:n_Λ_s]

    val = ∫(EIᵨ * τ * Δ(v) * Δ(ηₜ))dom[:dΓ_s] +
          ∫(EIᵨ * τ * (
              -jump(∇(v) ⋅ n_Λ) * mean(Δ(ηₜ))
              - mean(Δ(v)) * jump(∇(ηₜ) ⋅ n_Λ)
              + (γ / h) * jump(∇(v) ⋅ n_Λ) * jump(∇(ηₜ) ⋅ n_Λ)))dom[:dΛ_s]
    if s.bndType isa FixedBoundary
        val += ∫(-EIᵨ * τ * v * Δ(ηₜ) * dom[:n_Λ_sb])dom[:dΛ_sb]
    end
    return val
end

function stiffness(s::Beam2D, dom::WeakFormDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    EIᵨ = s.EI / s.ρw
    γ   = dom[:γ_s]
    h   = dom[:h_s]
    n_Λ = dom[:n_Λ_s]

    val = ∫(v * (s.g * η) + EIᵨ * Δ(v) * Δ(η))dom[:dΓ_s] +
          ∫(EIᵨ * (
              -jump(∇(v) ⋅ n_Λ) * mean(Δ(η))
              - mean(Δ(v)) * jump(∇(η) ⋅ n_Λ)
              + (γ / h) * jump(∇(v) ⋅ n_Λ) * jump(∇(η) ⋅ n_Λ)))dom[:dΛ_s]
    if s.bndType isa FixedBoundary
        val += ∫(-EIᵨ * v * Δ(η) * dom[:n_Λ_sb])dom[:dΛ_sb]
    end
    return val
end

function rhs(s::Beam2D, dom::WeakFormDomains, f, y)
    sym = variable_symbol(s)
    v = y[sym]
    ∫(v * f[sym])dom[:dΓ_s]
end
