# =============================================================================
# TensionedEulerBernoulliBeam.jl
#
# Combines membrane pre-tension stiffness with Euler-Bernoulli bending
# stiffness in a single AbstractHydroelasticStructure. Reduces exactly to
# Membrane when EIᵨ = 0, and to EulerBernoulliBeam when Tᵨ = 0 — see
# test/examples/TensionedEulerBernoulliBeamWeakFormTests.jl and
# docs/src/guide/theory.md.
#
# Reuses `_eb_bending_stiffness_operator` and `_joint_stiffness_form` from
# EulerBernoulliBeam.jl, so this file only combines them with the membrane
# tension term — no bending/skeleton logic is duplicated.
# =============================================================================

"""
    TensionedEulerBernoulliBeam <: AbstractHydroelasticStructure

Parameters for a 2D beam combining Euler-Bernoulli bending stiffness with
membrane-like axial pre-tension, normalised by the ambient fluid density
`ρw`. The governing equation is

```math
m\\,\\eta_{tt} + EI\\,\\Delta^2\\eta - \\nabla\\cdot(T\\nabla\\eta) = p .
```

This model reduces exactly to [`Membrane`](@ref) when `EIᵨ = 0`, and to
[`EulerBernoulliBeam`](@ref) when `Tᵨ = 0` — the assembled system matrices
agree with the corresponding pure `Membrane`/`EulerBernoulliBeam` matrices
to machine precision, since `stiffness_operator` is literally the sum of a
[`Membrane`](@ref)-style tension term and an
[`EulerBernoulliBeam`](@ref)-style bending term (via the same
`_eb_bending_stiffness_operator` helper). See
`docs/src/guide/theory.md` for the governing equation, non-dimensionalisation,
and asymptotic (tension-dominated / bending-dominated) regimes.

`mass`, `damping`, `stiffness`, and `rhs` are inherited from
[`AbstractHydroelasticStructure`](@ref); this file supplies the combined
elastic operator via [`stiffness_operator`](@ref) and the joint contribution
via [`extra_stiffness_form`](@ref) (not subject to Rayleigh damping).

# Fields
- `L::Float64`    — Beam span [m]
- `mᵨ::Float64`   — Mass per unit length / ρw [m²]
- `EIᵨ::Union{Float64,Function}` — Flexural rigidity / ρw [m⁵/s²]; may be a
  scalar `Float64` or a univariate `Function(x) -> Float64` for spatially
  varying stiffness.
- `Tᵨ::Float64`   — Pre-tension / ρw [m³/s²]
- `τ::Float64`    — Stiffness-proportional structural damping coefficient
  (applied to both the bending and tension parts of `stiffness_operator`);
  default 0
- `g::Float64`    — Gravitational acceleration [m/s²]; default 9.81
- `joints::Vector{JointRotationalSpring}` — Rotational spring connections at
  interior skeleton facets; leave empty for a continuous beam.
- `symbol::Symbol` — Field unknown symbol; default `:η_tb`
- `space_domain_symbol::Symbol` — Triangulation key used for FE spaces; default `:Γη`
- `fe::FESpaceConfig` — FE discretisation parameters
- `ωn1::Union{Float64,Nothing}` — First dry analytical natural frequency
  [rad/s] for the canonical simply-supported 1D case, `√(EIᵨ/mᵨ·(π/L)⁴ +
  Tᵨ/mᵨ·(π/L)²)` (the same simply-supported convention as [`Membrane`](@ref)'s
  `ωn1`); `nothing` when `EIᵨ` is a `Function`.

# Example
```julia
ρ = 1025.0; g = 9.81
beam = TensionedEulerBernoulliBeam(L=10.0, mᵨ=922.5/ρ, EIᵨ=1e5/ρ, Tᵨ=10*ρ*g/ρ)
```

See also: [`Membrane`](@ref), [`EulerBernoulliBeam`](@ref),
[`JointRotationalSpring`](@ref)

# References
- [C23] Colomes, O., Verdugo, F., & Akkerman, I. (2023). A monolithic
  finite element formulation for the hydroelastic analysis of very large
  floating structures. *Int. J. Numer. Methods Eng.*, 124(3), 714-751.
  DOI: https://doi.org/10.1002/nme.7140
- [A24] Agarwal, S., Colomes, O., & Metrikine, A. V. (2024).
  Dynamic analysis of viscoelastic floating membranes using monolithic
  finite element method. *Journal of Fluids and Structures*, 129, 104167.
  DOI: https://doi.org/10.1016/j.jfluidstructs.2024.104167
"""
@with_kw struct TensionedEulerBernoulliBeam <: AbstractHydroelasticStructure
    L::Float64
    mᵨ::Float64
    EIᵨ::Union{Float64, Function}
    Tᵨ::Float64
    τ::Float64     = 0.0
    g::Float64     = 9.81
    joints::Vector{JointRotationalSpring} = JointRotationalSpring[]
    symbol::Symbol = :η_tb
    space_domain_symbol::Symbol = :Γη
    fe::FESpaceConfig = FESpaceConfig()

    # Derived quantity: dry natural frequency for the canonical simply-supported
    # 1D case, ω₁² = (EI/m)(π/L)⁴ + (T/m)(π/L)² — same convention as Membrane's ωn1.
    ωn1::Union{Float64, Nothing} = EIᵨ isa Float64 ?
        sqrt(EIᵨ / mᵨ * (π / L)^4 + Tᵨ / mᵨ * (π / L)^2) : nothing
end

function print_parameters(beam::TensionedEulerBernoulliBeam)
    @printf("\n[MSG] Tensioned Euler-Bernoulli Beam Properties:\n")
    @printf("[VAL] L = %.4f m\n", beam.L)
    @printf("[VAL] mᵨ = %.4f m\n", beam.mᵨ)
    beam.EIᵨ isa Float64 ? @printf("[VAL] EIᵨ = %.4f m5/s2\n", beam.EIᵨ) : @printf("[VAL] EIᵨ = <function>\n")
    @printf("[VAL] Tᵨ = %.4f m3/s2\n", beam.Tᵨ)
    @printf("[VAL] τ = %.4f \n", beam.τ)
    beam.ωn1 !== nothing ? @printf("[VAL] 1st Dry Analytical Natural Freq (simply-supported), ωn1 = %.4f rad/s\n", beam.ωn1) : @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = <not available for variable EIᵨ>\n")
    println()
end

variable_symbol(s::TensionedEulerBernoulliBeam) = s.symbol

mass_density(s::TensionedEulerBernoulliBeam) = s.mᵨ
damping_parameter(s::TensionedEulerBernoulliBeam) = s.τ

"""
    stiffness_operator(s::TensionedEulerBernoulliBeam, dom::IntegrationDomains, x, y)

Combined elastic stiffness operator: Euler-Bernoulli C/DG bending (via
`_eb_bending_stiffness_operator`, identical to
[`EulerBernoulliBeam`](@ref)'s) plus membrane pre-tension (identical in form
to [`Membrane`](@ref)'s `stiffness_operator`):

```math
\\int_{\\Gamma_\\eta} \\Big[ EI\\,\\Delta v\\,\\Delta\\eta \\;+\\; T\\,\\nabla v \\cdot \\nabla \\eta \\Big]\\,\\mathrm{d}\\Gamma
- \\int_{\\Lambda_\\eta} \\text{(bending consistency + symmetry + penalty terms)}
```

Setting `EIᵨ = 0` makes this identical to [`Membrane`](@ref)'s
`stiffness_operator`; setting `Tᵨ = 0` makes this identical to
[`EulerBernoulliBeam`](@ref)'s `stiffness_operator`. Since [`damping`](@ref)
scales this whole operator by `τ` (see [`AbstractHydroelasticStructure`](@ref)),
Rayleigh damping is applied proportionally to both the bending and tension
parts of the stiffness.
"""
function stiffness_operator(s::TensionedEulerBernoulliBeam, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    dΩ = _space_measure(dom, s)
    _eb_bending_stiffness_operator(s.EIᵨ, s.fe.γ, η, v, dom, dΩ) + ∫(s.Tᵨ * ∇(v) ⋅ ∇(η))dΩ
end

"""
    extra_stiffness_form(s::TensionedEulerBernoulliBeam, dom::IntegrationDomains, x, y)

Rotational-spring joint contributions at `s.joints`, via
`_joint_stiffness_form` (identical mechanism to
[`EulerBernoulliBeam`](@ref)'s). Not subject to Rayleigh damping.
"""
function extra_stiffness_form(s::TensionedEulerBernoulliBeam, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    _joint_stiffness_form(s.joints, x[sym], y[sym], dom)
end
