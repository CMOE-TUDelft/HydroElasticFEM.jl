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

function mass(s::Membrane, dom::IntegrationDomains, x_tt, y)
  sym = variable_symbol(s)
  ηₜₜ = x_tt[sym]
  v = y[sym]
  ∫(s.mᵨ * v * ηₜₜ)dom[:dΓη]
end

function damping(s::Membrane, dom::IntegrationDomains, x_t, y)
  sym = variable_symbol(s)
  ηₜ = x_t[sym]
  v = y[sym]
  ∫(s.Tᵨ * s.τ * ∇(v) ⋅ ∇(ηₜ))dom[:dΓη]
end

function stiffness(s::Membrane, dom::IntegrationDomains, x, y)
  sym = variable_symbol(s)
  η = x[sym]
  v = y[sym]
  ∫(v * (s.g * η) + s.Tᵨ * ∇(v) ⋅ ∇(η))dom[:dΓη]
end

function rhs(s::Membrane, dom::IntegrationDomains, f, y)
  sym = variable_symbol(s)
  v = y[sym]
  ∫(v * f[sym])dom[:dΓη]
end
