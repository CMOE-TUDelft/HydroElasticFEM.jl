"""
    EulerBernoulliBeam <: AbstractStructure

Parameters for a 2D Euler-Bernoulli beam model (no joints),
normalised by fluid density ρw.

# Fields
- `L::Float64`    — Length of beam
- `mᵨ::Float64`   — Mass per unit length per unit width / ρw
- `EIᵨ::Float64`  — Flexural Rigidity / ρw
- `τ::Float64`    — Stiffness Proportional Structural Damping coefficient
- `g::Float64`    — Gravitational acceleration
- `bndType::BoundaryCondition` — Boundary Type
- `ωn1::Float64`   — Dry Analytical Natural frequency (derived)
"""
@with_kw struct EulerBernoulliBeam <: AbstractStructure
    L::Float64
    mᵨ::Float64
    EIᵨ::Float64
    τ::Float64     = 0.0
    g::Float64     = 9.81
    bndType::BoundaryCondition = FreeBoundary()

    # Derived quantities
    ωn1::Float64   = 22.3733 * sqrt(EIᵨ / (mᵨ * L^4))
    fe::FESpaceConfig = FESpaceConfig()
end

function print_parameters(beam::EulerBernoulliBeam)
    @printf("\n[MSG] Beam Properties:\n")
    @printf("[VAL] L = %.4f m\n", beam.L)
    @printf("[VAL] mᵨ = %.4f m\n", beam.mᵨ)
    @printf("[VAL] EIᵨ = %.4f m5/s2\n", beam.EIᵨ)
    @printf("[VAL] τ = %.4f \n", beam.τ)
    @printf("[VAL] beamBndType = %s \n", string(beam.bndType))
    @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s\n", beam.ωn1)
    @printf("[MSG] See free-free vibration frequency formula involving 22.3733 from Wiki.\n")
    println()
end

variable_symbol(::EulerBernoulliBeam) = :η_b

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only η_b terms — no coupling to ϕ or other fields

function mass(s::EulerBernoulliBeam, dom::WeakFormDomains, x_tt, y)
    sym = variable_symbol(s)
    ηₜₜ = x_tt[sym]
    v   = y[sym]
    ∫(s.mᵨ * v * ηₜₜ)dom[:dΓ_s]
end

function damping(s::EulerBernoulliBeam, dom::WeakFormDomains, x_t, y)
    sym = variable_symbol(s)
    ηₜ = x_t[sym]
    v  = y[sym]
    EIᵨ = s.EIᵨ
    τ   = s.τ
    γ   = s.fe.γ
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

function stiffness(s::EulerBernoulliBeam, dom::WeakFormDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    EIᵨ = s.EIᵨ
    γ   = s.fe.γ
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

function rhs(s::EulerBernoulliBeam, dom::WeakFormDomains, f, y)
    sym = variable_symbol(s)
    v = y[sym]
    ∫(v * f[sym])dom[:dΓ_s]
end
