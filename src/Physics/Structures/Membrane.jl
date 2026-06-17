"""
    Membrane <: Structure

Parameters for an nD membrane model, normalised by fluid density `ρw`.

The membrane manifold can be 1D or 2D, embedded in a 2D or 3D fluid.
Weak forms are written in a dimension-agnostic way and evaluated on
`space_domain_symbol`.  The structural damping form uses stiffness-proportional
Rayleigh damping with coefficient `τ`.

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
@with_kw struct Membrane <: Structure
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

"""
    mass(s::Membrane, dom::IntegrationDomains, x_tt, y)

Membrane inertia (mass) bilinear form.

Assembles:
```math
\\int_{\\Gamma_\\eta} m_\\varrho \\, v \\, \\partial_{tt}\\eta \\, \\mathrm{d}\\Gamma_\\eta
```

# Arguments
- `s::Membrane`: membrane parameters (provides `mρ`)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`)
- `x_tt`: second time-derivative trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`

# Reference
[A24] Agarwal et al. (2024), J. Fluids Struct., 129, 104167.
"""
function mass(s::Membrane, dom::IntegrationDomains, x_tt, y)
  sym = variable_symbol(s)
  ηₜₜ = x_tt[sym]
  v = y[sym]
  dΩ = _space_measure(dom, s)
  ∫(s.mᵨ * v * ηₜₜ)dΩ
end

"""
    damping(s::Membrane, dom::IntegrationDomains, x_t, y)

Membrane stiffness-proportional Rayleigh damping bilinear form.

Assembles:
```math
\\int_{\\Gamma_\\eta} T_\\varrho \\tau \\, \\nabla v \\cdot \\nabla(\\partial_t\\eta) \\, \\mathrm{d}\\Gamma_\\eta
```

The coefficient `τ` is the stiffness-proportional Rayleigh damping parameter.
Set `τ = 0` (default) to disable structural damping.

# Arguments
- `s::Membrane`: membrane parameters (provides `Tρ`, `τ`)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`)
- `x_t`: first time-derivative trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`

# Reference
[A24] Agarwal et al. (2024), J. Fluids Struct., 129, 104167.
"""
function damping(s::Membrane, dom::IntegrationDomains, x_t, y)
  sym = variable_symbol(s)
  ηₜ = x_t[sym]
  v = y[sym]
  dΩ = _space_measure(dom, s)
  ∫(s.Tᵨ * s.τ * ∇(v) ⋅ ∇(ηₜ))dΩ
end

"""
    stiffness(s::Membrane, dom::IntegrationDomains, x, y)

Membrane structural stiffness bilinear form (hydrostatic restoring + pre-tension).

Assembles:
```math
\\int_{\\Gamma_\\eta} \\bigl( g \\, v \\, \\eta + T_\\varrho \\, \\nabla v \\cdot \\nabla \\eta \\bigr) \\, \\mathrm{d}\\Gamma_\\eta
```

# Arguments
- `s::Membrane`: membrane parameters (provides `g`, `Tρ`)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`)
- `x`: trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`

# Reference
[A24] Agarwal et al. (2024), J. Fluids Struct., 129, 104167.
"""
function stiffness(s::Membrane, dom::IntegrationDomains, x, y)
  sym = variable_symbol(s)
  η = x[sym]
  v = y[sym]
  dΩ = _space_measure(dom, s)
  # Membrane structural bilinear form:
  # ∫_Γη (g·v·η + Tρ·∇v·∇η) dΓ.
  # Reference: [A24] (viscoelastic floating membrane formulation).
  ∫(v * (s.g * η) + s.Tᵨ * ∇(v) ⋅ ∇(η))dΩ
end

"""
    rhs(s::Membrane, dom::IntegrationDomains, f, y)

Membrane right-hand side (applied load) linear form.

Assembles the body-force or pressure load contribution:
```math
\\int_{\\Gamma_\\eta} v \\, f_\\eta \\, \\mathrm{d}\\Gamma_\\eta
```

# Arguments
- `s::Membrane`: membrane parameters (provides `symbol` for field lookup)
- `dom::IntegrationDomains`: integration measures (requires `:dΓη`)
- `f`: forcing `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`
"""
function rhs(s::Membrane, dom::IntegrationDomains, f, y)
  sym = variable_symbol(s)
  v = y[sym]
  dΩ = _space_measure(dom, s)
  ∫(v * f[sym])dΩ
end
