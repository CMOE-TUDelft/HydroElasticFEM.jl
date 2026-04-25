"""
    JointRotationalSpring

Rotational spring stiffness contribution at an interior joint of an
`EulerBernoulliBeam`.  Each joint adds the term

```math
\\int_{\\Lambda_j} k_r \\,
[\\![ \\nabla(v) \\cdot n_{\\Lambda_j} ]\\!] \\,
[\\![ \\nabla(\\eta) \\cdot n_{\\Lambda_j} ]\\!]
\\, \\mathrm{d}\\Lambda_j
```

to the beam stiffness form, where ``[\\![\\cdot]\\!]`` denotes the jump across
the skeleton facet ``\\Lambda_j`` and ``n_{\\Lambda_j}`` is its outward normal.

The `domain_symbol` and `normal_symbol` must match the keys registered in
`IntegrationDomains` — this is done automatically by `get_integration_domains`
when the corresponding `JointDomain1D` is declared in `TankDomain2D`.

# Fields
- `domain_symbol::Symbol` — Key for the joint skeleton `Measure` in
  `IntegrationDomains` (e.g. `:dΛj_1`).
- `normal_symbol::Symbol` — Key for the joint outward-normal field in
  `IntegrationDomains` (e.g. `:n_Λ_j_1`).
- `kᵣ::Float64` — Rotational spring stiffness ``[\\mathrm{N}\\cdot\\mathrm{m}^2/\\rho_w]``.

# See also
[`JointDomain1D`](@ref HydroElasticFEM.Geometry.JointDomain1D),
[`EulerBernoulliBeam`](@ref)
"""
struct JointRotationalSpring
    domain_symbol::Symbol
    normal_symbol::Symbol
    kᵣ::Float64
end

"""
    EulerBernoulliBeam <: Structure

Parameters for a 2D Euler-Bernoulli beam model,
normalised by fluid density ρw.

# Fields
- `L::Float64`    — Length of beam
- `mᵨ::Float64`   — Mass per unit length per unit width / ρw
- `EIᵨ::Float64`  — Flexural Rigidity / ρw
- `τ::Float64`    — Stiffness Proportional Structural Damping coefficient
- `g::Float64`    — Gravitational acceleration
- `joints::Vector{JointRotationalSpring}` — Optional rotational spring joints
- `bndType::BoundaryCondition` — Boundary Type
- `ωn1::Float64`   — Dry Analytical Natural frequency (derived)
"""
@with_kw struct EulerBernoulliBeam <: Structure
    L::Float64
    mᵨ::Float64
    EIᵨ::Float64
    τ::Float64     = 0.0
    g::Float64     = 9.81
    joints::Vector{JointRotationalSpring} = JointRotationalSpring[]
    symbol::Symbol = :η_b
    space_domain_symbol::Symbol = :Γη
    fe::FESpaceConfig = FESpaceConfig()

    # Derived quantities
    ωn1::Float64   = 22.3733 * sqrt(EIᵨ / (mᵨ * L^4))
end

function print_parameters(beam::EulerBernoulliBeam)
    @printf("\n[MSG] Beam Properties:\n")
    @printf("[VAL] L = %.4f m\n", beam.L)
    @printf("[VAL] mᵨ = %.4f m\n", beam.mᵨ)
    @printf("[VAL] EIᵨ = %.4f m5/s2\n", beam.EIᵨ)
    @printf("[VAL] τ = %.4f \n", beam.τ)
    @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s\n", beam.ωn1)
    @printf("[MSG] See free-free vibration frequency formula involving 22.3733 from Wiki.\n")
    println()
end

variable_symbol(s::EulerBernoulliBeam) = s.symbol

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only η_b terms — no coupling to ϕ or other fields

function mass(s::EulerBernoulliBeam, dom::IntegrationDomains, x_tt, y)
    sym = variable_symbol(s)
    ηₜₜ = x_tt[sym]
    v   = y[sym]
    ∫(s.mᵨ * v * ηₜₜ)dom[:dΓη]
end

function damping(s::EulerBernoulliBeam, dom::IntegrationDomains, x_t, y)
    sym = variable_symbol(s)
    ηₜ = x_t[sym]
    v  = y[sym]
    EIᵨ = s.EIᵨ
    τ   = s.τ
    γ   = s.fe.γ
    h   = dom[:h_η]
    n_Λ = dom[:n_Λ_η]

    val = ∫(EIᵨ * τ * Δ(v) * Δ(ηₜ))dom[:dΓη] +
          ∫(EIᵨ * τ * (
              -jump(∇(v) ⋅ n_Λ) * mean(Δ(ηₜ))
              - mean(Δ(v)) * jump(∇(ηₜ) ⋅ n_Λ)
              + (γ / h) * jump(∇(v) ⋅ n_Λ) * jump(∇(ηₜ) ⋅ n_Λ)))dom[:dΛη]
    return val
end

function stiffness(s::EulerBernoulliBeam, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    EIᵨ = s.EIᵨ
    γ   = s.fe.γ
    h   = dom[:h_η]
    n_Λ = dom[:n_Λ_η]

    val = ∫(v * (s.g * η) + EIᵨ * Δ(v) * Δ(η))dom[:dΓη] +
          ∫(EIᵨ * (
              -jump(∇(v) ⋅ n_Λ) * mean(Δ(η))
              - mean(Δ(v)) * jump(∇(η) ⋅ n_Λ)
              + (γ / h) * jump(∇(v) ⋅ n_Λ) * jump(∇(η) ⋅ n_Λ)))dom[:dΛη]

    for joint in s.joints
        dΛj  = dom[joint.domain_symbol]
        n_Λj = dom[joint.normal_symbol]
        val  = _add_contribution(val, ∫(joint.kᵣ * jump(∇(v) ⋅ n_Λj) * jump(∇(η) ⋅ n_Λj))dΛj)
    end

    return val
end

function rhs(s::EulerBernoulliBeam, dom::IntegrationDomains, f, y)
    sym = variable_symbol(s)
    v = y[sym]
    ∫(v * f[sym])dom[:dΓη]
end
