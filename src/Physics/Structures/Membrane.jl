"""
    Membrane <: AbstractHydroelasticStructure

Parameters for an nD membrane model, normalised by fluid density `ρw`.

The membrane manifold can be 1D or 2D, embedded in a 2D or 3D fluid.
Weak forms are written in a dimension-agnostic way and evaluated on
`space_domain_symbol`.  The structural damping form uses stiffness-proportional
Rayleigh damping with coefficient `τ`.

`mass`, `damping`, `stiffness`, and `rhs` are inherited from
[`AbstractHydroelasticStructure`](@ref); this file only supplies the
membrane's elastic operator, [`stiffness_operator`](@ref) (pre-tension), via
[`mass_density`](@ref) and [`damping_parameter`](@ref).

# Fields
- `L::Float64`         — Characteristic membrane length [m]
- `mᵨ::Float64`        — Mass per unit manifold measure / ρw [dimensionless]
- `Tᵨ::Float64`        — Pre-tension / ρw [m³/s²]
- `τ::Float64`         — Stiffness-proportional structural damping coefficient; default 0
- `g::Float64`         — Gravitational acceleration [m/s²]; default 9.81
- `ambient_dim::Int`   — Fluid ambient dimension: 2 or 3; default 2
- `manifold_dim::Int`  — Membrane manifold dimension: 1 or 2; default 1
- `symbol::Symbol`     — Field unknown symbol; default `:η_m`
- `space_domain_symbol::Symbol` — Triangulation key for FE spaces; default `:Γη`
- `fe::FESpaceConfig`  — FE discretisation parameters
- `ωn1::Float64`       — First dry analytical natural frequency [rad/s], derived as
  `(π/L) * √(Tᵨ/mᵨ)` for the 1D canonical case

# References
- [A24] Agarwal, S., Colomes, O., & Metrikine, A. V. (2024).
  Dynamic analysis of viscoelastic floating membranes using monolithic
  finite element method. *Journal of Fluids and Structures*, 129, 104167.
  DOI: https://doi.org/10.1016/j.jfluidstructs.2024.104167
- [C23] Colomes, O., Verdugo, F., & Akkerman, I. (2023). A monolithic
  finite element formulation for the hydroelastic analysis of very large
  floating structures. *Int. J. Numer. Methods Eng.*, 124(3), 714-751.
  DOI: https://doi.org/10.1002/nme.7140
"""
@with_kw struct Membrane <: AbstractHydroelasticStructure
  L::Float64
  mᵨ::Float64
  Tᵨ::Float64
  τ::Float64 = 0.0
  g::Float64 = 9.81
  ambient_dim::Int = 2
  manifold_dim::Int = 1
  symbol::Symbol = :η_m
  space_domain_symbol::Symbol = :Γη
  fe::FESpaceConfig = FESpaceConfig()

  # 1st dry analytical natural frequency for the canonical 1D case.
  ωn1::Float64 = (π / L) * sqrt(Tᵨ / mᵨ)
end

function print_parameters(memb::Membrane)
  @printf("\n[MSG] Membrane Properties:\n")
  @printf("[VAL] L = %.4f m\n", memb.L)
  @printf("[VAL] mᵨ = %.4f m\n", memb.mᵨ)
  @printf("[VAL] Tᵨ = %.4f m3/s2\n", memb.Tᵨ)
  @printf("[VAL] τ = %.4f \n", memb.τ)
  @printf("[VAL] ambient_dim = %d\n", memb.ambient_dim)
  @printf("[VAL] manifold_dim = %d\n", memb.manifold_dim)
  @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s \n", memb.ωn1)
  println()
end

variable_symbol(s::Membrane) = s.symbol
ambient_dimension(s::Membrane) = s.ambient_dim
manifold_dimension(s::Membrane) = s.manifold_dim

mass_density(s::Membrane) = s.mᵨ
damping_parameter(s::Membrane) = s.τ

"""
    stiffness_operator(m::Membrane, dom::IntegrationDomains, x, y)

Membrane pre-tension elastic stiffness operator.

Assembles:
```math
\\int_{\\Gamma_\\eta} T_\\varrho \\, \\nabla v \\cdot \\nabla \\eta \\, \\mathrm{d}\\Gamma_\\eta
```

Combined with the shared hydrostatic term (in [`stiffness`](@ref)) this
reproduces the classical membrane bilinear form
`∫ (g·v·η + Tρ·∇v·∇η) dΓ`; scaled by `τ` (in [`damping`](@ref)) it
reproduces the Rayleigh damping form `∫ Tρ·τ·∇v·∇ηₜ dΓ`.

# Arguments
- `m::Membrane`: membrane parameters (provides `Tρ`)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`)
- `x`: trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`

# Reference
[A24] Agarwal et al. (2024), J. Fluids Struct., 129, 104167.
"""
function stiffness_operator(m::Membrane, dom::IntegrationDomains, x, y)
  sym = variable_symbol(m)
  η = x[sym]
  v = y[sym]
  dΩ = _space_measure(dom, m)
  ∫(m.Tᵨ * ∇(v) ⋅ ∇(η))dΩ
end
