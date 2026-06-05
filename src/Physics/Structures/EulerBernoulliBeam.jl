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
when the corresponding `JointDomain` is declared in a 2D `TankDomain`.

# Fields
- `domain_symbol::Symbol` — Key for the joint skeleton `Measure` in
  `IntegrationDomains` (e.g. `:dΛj_1`).
- `normal_symbol::Symbol` — Key for the joint outward-normal field in
  `IntegrationDomains` (e.g. `:n_Λ_j_1`).
- `kᵣ::Float64` — Rotational spring stiffness ``[\\mathrm{N}\\cdot\\mathrm{m}^2/\\rho_w]``.

# See also
[`JointDomain`](@ref HydroElasticFEM.Geometry.JointDomain),
[`EulerBernoulliBeam`](@ref)
"""
struct JointRotationalSpring
    domain_symbol::Symbol
    normal_symbol::Symbol
    kᵣ::Float64
end

"""
    EulerBernoulliBeam <: Structure

Parameters for a 2D Euler-Bernoulli beam model, normalised by the ambient
fluid density `ρw`.

The interior-penalty C/DG formulation uses a Symmetric Interior Penalty (SIP)
consistency + penalty scheme for the fourth-order bending operator.

# Fields
- `L::Float64`    — Beam span [m]
- `mᵨ::Float64`   — Mass per unit length / ρw [m²], i.e. `ρb·hb / ρw`
- `EIᵨ::Union{Float64,Function}` — Flexural rigidity / ρw [m⁵/s²];
  may be a scalar `Float64` or a univariate `Function(x) -> Float64`
  for spatially varying stiffness.
- `τ::Float64`    — Stiffness-proportional structural damping coefficient; default 0
- `g::Float64`    — Gravitational acceleration [m/s²]; default 9.81
- `joints::Vector{JointRotationalSpring}` — Rotational spring connections at
  interior skeleton facets; leave empty for a continuous beam.
- `symbol::Symbol` — Field unknown symbol; default `:η_b`
- `space_domain_symbol::Symbol` — Triangulation key used for FE spaces; default `:Γη`
- `fe::FESpaceConfig` — FE discretisation parameters
- `ωn1::Union{Float64,Nothing}` — First dry analytical natural frequency [rad/s], derived
  from `EIᵨ` and `mᵨ` for clamped-free; `nothing` when `EIᵨ` is a `Function`

# Example
```julia
beam = EulerBernoulliBeam(L=10.0, mᵨ=0.5, EIᵨ=1.0e4)
```

See also: [`JointRotationalSpring`](@ref), [`JointDomain`](@ref)

# Reference
- [C23] Colomes, O., Verdugo, F., & Akkerman, I. (2023). A monolithic
    finite element formulation for the hydroelastic analysis of very large
    floating structures. *Int. J. Numer. Methods Eng.*, 124(3), 714-751.
    DOI: https://doi.org/10.1002/nme.7140
"""
@with_kw struct EulerBernoulliBeam <: Structure
    L::Float64
    mᵨ::Float64
    EIᵨ::Union{Float64, Function}
    τ::Float64     = 0.0
    g::Float64     = 9.81
    joints::Vector{JointRotationalSpring} = JointRotationalSpring[]
    symbol::Symbol = :η_b
    space_domain_symbol::Symbol = :Γη
    fe::FESpaceConfig = FESpaceConfig()

    # Derived quantities
    ωn1::Union{Float64, Nothing} = EIᵨ isa Float64 ? 22.3733 * sqrt(EIᵨ / (mᵨ * L^4)) : nothing
end

function print_parameters(beam::EulerBernoulliBeam)
    @printf("\n[MSG] Beam Properties:\n")
    @printf("[VAL] L = %.4f m\n", beam.L)
    @printf("[VAL] mᵨ = %.4f m\n", beam.mᵨ)
    beam.EIᵨ isa Float64 ? @printf("[VAL] EIᵨ = %.4f m5/s2\n", beam.EIᵨ) : @printf("[VAL] EIᵨ = <function>\n")
    @printf("[VAL] τ = %.4f \n", beam.τ)
    beam.ωn1 !== nothing ? @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s\n", beam.ωn1) : @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = <not available for variable EIᵨ>\n")
    @printf("[MSG] See free-free vibration frequency formula involving 22.3733 from Wiki.\n")
    println()
end

variable_symbol(s::EulerBernoulliBeam) = s.symbol

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only η_b terms — no coupling to ϕ or other fields

"""
    mass(s::EulerBernoulliBeam, dom::IntegrationDomains, x_tt, y)

Euler-Bernoulli beam inertia (mass) bilinear form.

Assembles:
```math
\\int_{\\Gamma_\\eta} m_\\varrho \\, v \\, \\partial_{tt}\\eta \\, \\mathrm{d}\\Gamma_\\eta
```

# Arguments
- `s::EulerBernoulliBeam`: beam parameters (provides `mρ`)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`)
- `x_tt`: second time-derivative trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`

# Reference
[C23] Colomés et al. (2023), Section 3.1, Eq. (16).
"""
function mass(s::EulerBernoulliBeam, dom::IntegrationDomains, x_tt, y)
    sym = variable_symbol(s)
    ηₜₜ = x_tt[sym]
    v   = y[sym]
    dΩ = _space_measure(dom, s)
    ∫(s.mᵨ * v * ηₜₜ)dΩ
end

"""
    damping(s::EulerBernoulliBeam, dom::IntegrationDomains, x_t, y)

Euler-Bernoulli beam stiffness-proportional Rayleigh damping bilinear form.

Uses the symmetric interior-penalty C/DG formulation with penalty parameter `γ`:

```math
\\int_{\\Gamma_\\eta} EI\\tau \\Delta v \\Delta\\partial_t\\eta \\, \\mathrm{d}\\Gamma_\\eta
- \\int_{\\Lambda_\\eta} \\text{(consistency + symmetry + penalty skeleton terms)}
```

# Arguments
- `s::EulerBernoulliBeam`: beam parameters (provides `EIρ`, `τ`, `fe.gamma`)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`, `:dΛη`, `:h_η`, `:n_Λ_η`)
- `x_t`: first time-derivative trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`

# Reference
[C23] Colomés et al. (2023), Section 3.1, Eq. (16)-(20).
"""
function damping(s::EulerBernoulliBeam, dom::IntegrationDomains, x_t, y)
    sym = variable_symbol(s)
    ηₜ = x_t[sym]
    v  = y[sym]
    τ   = s.τ
    trian = get_triangulation(v)
    EIτ_param = s.EIᵨ isa Float64 ? s.EIᵨ * τ : (x -> s.EIᵨ(x) * τ)
    EIτ = materialize(EIτ_param, trian)
    γ   = s.fe.γ
    h   = dom[:h_η]
    n_Λ = dom[:n_Λ_η]
    dΩ  = _space_measure(dom, s)

    val = ∫(EIτ * Δ(v) * Δ(ηₜ))dΩ +
          ∫(
              -jump(∇(v) ⋅ n_Λ) * mean(EIτ * Δ(ηₜ))
              - mean(EIτ * Δ(v)) * jump(∇(ηₜ) ⋅ n_Λ)
              + (γ / h) * mean(EIτ) * jump(∇(v) ⋅ n_Λ) * jump(∇(ηₜ) ⋅ n_Λ))dom[:dΛη]
    return val
end

"""
    stiffness(s::EulerBernoulliBeam, dom::IntegrationDomains, x, y)

Euler-Bernoulli beam stiffness bilinear form (gravity + C/DG bending + joints).

Assembles the bulk bending term, symmetric interior-penalty skeleton terms,
and optional rotational-spring contributions at declared joints.

# Arguments
- `s::EulerBernoulliBeam`: beam parameters (`EIᵨ`, `g`, `fe.γ`, optional `joints`)
- `dom::IntegrationDomains`: integration domains (`:dΓη`, `:dΛη`, `:h_η`, `:n_Λ_η`, and joint keys)
- `x`: trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`

# Reference
[C23] Colomés et al. (2023), Section 3.1, Eq. (16)-(20).
"""
function stiffness(s::EulerBernoulliBeam, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    trian = get_triangulation(v)
    EI  = materialize(s.EIᵨ, trian)
    γ   = s.fe.γ
    h   = dom[:h_η]
    n_Λ = dom[:n_Λ_η]
    dΩ  = _space_measure(dom, s)

    # Euler-Bernoulli C/DG bending formulation on Γb and Skeleton(Γb).
    # Bulk: ∫_Γb a1·Δη·Δv dΓ, with a1 = EI/ρ.
    # Skeleton: consistency + symmetry + penalty terms.
    # Reference: [C23] Section 3.1, Eq. (16)-(20).
    val = ∫(v * (s.g * η) + EI * Δ(v) * Δ(η))dΩ +
          ∫(
              -jump(∇(v) ⋅ n_Λ) * mean(EI * Δ(η))
              - mean(EI * Δ(v)) * jump(∇(η) ⋅ n_Λ)
              + (γ / h) * mean(EI) * jump(∇(v) ⋅ n_Λ) * jump(∇(η) ⋅ n_Λ))dom[:dΛη]

    for joint in s.joints
        dΛj  = dom[joint.domain_symbol]
        n_Λj = dom[joint.normal_symbol]
        val  = _add_contribution(val, ∫(joint.kᵣ * jump(∇(v) ⋅ n_Λj) * jump(∇(η) ⋅ n_Λj))dΛj)
    end

    return val
end

"""
    rhs(s::EulerBernoulliBeam, dom::IntegrationDomains, f, y)

Euler-Bernoulli beam right-hand side (applied load) linear form.

Assembles the distributed load contribution:
```math
\\int_{\\Gamma_\\eta} v \\, f_\\eta \\, \\mathrm{d}\\Gamma_\\eta
```

# Arguments
- `s::EulerBernoulliBeam`: beam parameters (provides `symbol` for field lookup)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`)
- `f`: forcing `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`
"""
function rhs(s::EulerBernoulliBeam, dom::IntegrationDomains, f, y)
    sym = variable_symbol(s)
    v = y[sym]
    dΩ = _space_measure(dom, s)
    ∫(v * f[sym])dΩ
end
