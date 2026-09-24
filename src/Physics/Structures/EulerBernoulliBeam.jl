"""
    JointRotationalSpring

Rotational spring stiffness contribution at an interior joint of an
`EulerBernoulliBeam` (or a [`TensionedEulerBernoulliBeam`](@ref)).  Each
joint adds the term

```math
\\int_{\\Lambda_j} k_r \\,
[\\![ \\nabla(v) \\cdot n_{\\Lambda_j} ]\\!] \\,
[\\![ \\nabla(\\eta) \\cdot n_{\\Lambda_j} ]\\!]
\\, \\mathrm{d}\\Lambda_j
```

to the beam stiffness form, where ``[\\![\\cdot]\\!]`` denotes the jump across
the skeleton facet ``\\Lambda_j`` and ``n_{\\Lambda_j}`` is its outward normal.
This contribution is *not* subject to stiffness-proportional Rayleigh
damping (see [`extra_stiffness_form`](@ref)): joints are purely elastic
connections.

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
[`EulerBernoulliBeam`](@ref), [`TensionedEulerBernoulliBeam`](@ref)
"""
struct JointRotationalSpring
    domain_symbol::Symbol
    normal_symbol::Symbol
    kᵣ::Float64
end

"""
    _joint_stiffness_form(joints, η, v, dom)

Shared rotational-spring joint penalty contribution, reused by both
[`EulerBernoulliBeam`](@ref) and [`TensionedEulerBernoulliBeam`](@ref) (both
discretised with the same C/DG skeleton and `JointRotationalSpring`
connections). Returns `nothing` when `joints` is empty, matching the
`_add_contribution` no-op convention.
"""
function _joint_stiffness_form(joints::Vector{JointRotationalSpring}, η, v, dom::IntegrationDomains)
    val = nothing
    for joint in joints
        dΛj  = dom[joint.domain_symbol]
        n_Λj = dom[joint.normal_symbol]
        val  = _add_contribution(val, ∫(joint.kᵣ * jump(∇(v) ⋅ n_Λj) * jump(∇(η) ⋅ n_Λj))dΛj)
    end
    return val
end

"""
    _eb_bending_stiffness_operator(EIᵨ, γ, η, v, dom, dΩ)

Shared Euler-Bernoulli C/DG bending stiffness operator (bulk + interior-penalty
skeleton), reused by both [`EulerBernoulliBeam`](@ref) and
[`TensionedEulerBernoulliBeam`](@ref):

```math
\\int_{\\Gamma_\\eta} EI\\,\\Delta v\\,\\Delta \\eta\\,\\mathrm{d}\\Gamma
- \\int_{\\Lambda_\\eta} \\text{(consistency + symmetry + penalty skeleton terms)}
```

`EIᵨ` may be a `Float64` or a univariate `Function(x) -> Float64`; it is
materialised on `get_triangulation(v)`. `γ` is the Nitsche/SIP penalty
parameter (`fe.γ`).

# Reference
[C23] Colomés et al. (2023), Section 3.1, Eq. (16)-(20).
"""
function _eb_bending_stiffness_operator(EIᵨ, γ::Float64, η, v, dom::IntegrationDomains, dΩ)
    trian = get_triangulation(v)
    EI  = materialize(EIᵨ, trian)
    h   = dom[:h_η]
    n_Λ = dom[:n_Λ_η]

    ∫(EI * Δ(v) * Δ(η))dΩ +
    ∫(
        -jump(∇(v) ⋅ n_Λ) * mean(EI * Δ(η))
        - mean(EI * Δ(v)) * jump(∇(η) ⋅ n_Λ)
        + (γ / h) * mean(EI) * jump(∇(v) ⋅ n_Λ) * jump(∇(η) ⋅ n_Λ))dom[:dΛη]
end

"""
    EulerBernoulliBeam <: AbstractHydroelasticStructure

Parameters for a 2D Euler-Bernoulli beam model, normalised by the ambient
fluid density `ρw`.

The interior-penalty C/DG formulation uses a Symmetric Interior Penalty (SIP)
consistency + penalty scheme for the fourth-order bending operator.

`mass`, `damping`, `stiffness`, and `rhs` are inherited from
[`AbstractHydroelasticStructure`](@ref); this file supplies the beam's
elastic operator via [`stiffness_operator`](@ref) (bending, via
`_eb_bending_stiffness_operator`) and the joint contribution via
[`extra_stiffness_form`](@ref) (which is *not* subject to Rayleigh damping).

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

See also: [`JointRotationalSpring`](@ref), [`JointDomain`](@ref),
[`TensionedEulerBernoulliBeam`](@ref)

# Reference
- [C23] Colomes, O., Verdugo, F., & Akkerman, I. (2023). A monolithic
    finite element formulation for the hydroelastic analysis of very large
    floating structures. *Int. J. Numer. Methods Eng.*, 124(3), 714-751.
    DOI: https://doi.org/10.1002/nme.7140
"""
@with_kw struct EulerBernoulliBeam <: AbstractHydroelasticStructure
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

mass_density(s::EulerBernoulliBeam) = s.mᵨ
damping_parameter(s::EulerBernoulliBeam) = s.τ

"""
    stiffness_operator(s::EulerBernoulliBeam, dom::IntegrationDomains, x, y)

Euler-Bernoulli C/DG bending stiffness operator (bulk + interior-penalty
skeleton), via `_eb_bending_stiffness_operator`.

Combined with the shared hydrostatic term (in [`stiffness`](@ref)) and the
joint contribution (in [`extra_stiffness_form`](@ref)) this reproduces the
original beam bilinear form; scaled by `τ` (in [`damping`](@ref)) it
reproduces the beam's stiffness-proportional Rayleigh damping form exactly,
since `τ` commutes linearly through `mean`/`jump`/`∫`.

# Reference
[C23] Colomés et al. (2023), Section 3.1, Eq. (16)-(20).
"""
function stiffness_operator(s::EulerBernoulliBeam, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    dΩ = _space_measure(dom, s)
    _eb_bending_stiffness_operator(s.EIᵨ, s.fe.γ, η, v, dom, dΩ)
end

"""
    extra_stiffness_form(s::EulerBernoulliBeam, dom::IntegrationDomains, x, y)

Rotational-spring joint contributions at `s.joints`, via
`_joint_stiffness_form`. Not subject to Rayleigh damping (joints are
purely elastic connections).
"""
function extra_stiffness_form(s::EulerBernoulliBeam, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    _joint_stiffness_form(s.joints, x[sym], y[sym], dom)
end
