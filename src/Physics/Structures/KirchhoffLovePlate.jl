# =============================================================================
# KirchhoffLovePlate — implementation status (last updated 2026-05-08)
#
# A. KirchhoffLovePlate struct EXISTS here (fields E, ν, hb, ρ, ρb, g,
#    ambient_dim, manifold_dim, symbol, space_domain_symbol, fe, C).
#
# B. build_kl_tensor / build_KL_tensor EXISTS here.
#    build_kl_tensor(ambient_dim, manifold_dim, E, ν, hb, ρ) builds the
#    SymFourthOrderTensorValue{ambient_dim} constitutive tensor C = E·I/(ρ·...)
#    The scalar C[1,1,1,1] = D/ρ where D = E·h³/(12(1-ν²)).
#
# C. REGISTERED in Physics.jl via
#    include("Structures/KirchhoffLovePlate.jl") and exported in
#    HydroElasticFEM.jl.  Inherits Structure <: PhysicsParameters so the
#    PotentialFlow↔Structure coupling damping in CouplingTerms.jl applies
#    automatically — no additional coupling code is needed.
#
# D. mass / damping / stiffness / rhs were NOT implemented before this edit.
#    They are added below.
#
# E. BeamPlateConsistencyTests.jl has three tests (1D beam-plate consistency
#    on the tensor scalar C[1,1,1,1]).  No weak-form assembly tests existed.
#    New weak-form validation tests live in
#    test/Physics/KirchhoffLovePlateTests.jl.
#
# Reference: [C23] Colomes, Verdugo, Akkerman (2023), NME.
#   Section 3.2, eqs. (24)–(25).
#   DOI: 10.1002/nme.7140
# =============================================================================

using Gridap.TensorValues

"""
    build_kl_tensor(ambient_dim, manifold_dim, E, ν, hb, ρ)

Build an nD Kirchhoff-Love constitutive tensor per fluid density.

- `ambient_dim` is the embedding space dimension (2 or 3)
- `manifold_dim` is the plate manifold dimension (1 or 2)

For `manifold_dim == 1`, the model reduces to Euler-Bernoulli curvature
rigidity with `C[1,1,1,1] = E*I/ρ`.
"""
function build_kl_tensor(
  ambient_dim::Int,
  manifold_dim::Int,
  E::Float64,
  ν::Float64,
  hb::Float64,
  ρ::Float64,
)
  ambient_dim in (2, 3) ||
    error("ambient_dim must be 2 or 3, got $ambient_dim")
  manifold_dim in (1, 2) ||
    error("manifold_dim must be 1 or 2, got $manifold_dim")
  manifold_dim <= ambient_dim ||
    error("manifold_dim ($manifold_dim) must be <= ambient_dim ($ambient_dim)")

  I = hb^3 / 12
  δ(x, y) = ==(x, y)
  C_type = SymFourthOrderTensorValue{ambient_dim, Float64}

  max_index = 0
  for i in 1:ambient_dim, j in 1:ambient_dim, k in 1:ambient_dim, l in 1:ambient_dim
    max_index = max(max_index, data_index(C_type, i, j, k, l))
  end
  Cvals = zeros(Float64, max_index)

  if manifold_dim == 1
    Cvals[data_index(C_type, 1, 1, 1, 1)] = E * I / ρ
  else
    μ = E / (2 * (1 + ν))
    λ = ν * E / (1 - ν^2)
    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
      Cvals[data_index(C_type, i, j, k, l)] =
        I / ρ * (
          μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) +
          λ * (δ(i, j) * δ(k, l))
        )
    end
  end

  SymFourthOrderTensorValue(Cvals...)
end

"""
    build_KL_tensor(E, ν, hb, ρ)

Backward-compatible constructor for the classical 2D manifold in 3D ambient
Kirchhoff-Love tensor.
"""
build_KL_tensor(E, ν, hb, ρ) =
  build_kl_tensor(3, 2, Float64(E), Float64(ν), Float64(hb), Float64(ρ))

"""
    check_major_symmetry(C; atol=1e-10, dim=3)

Verify major symmetry `C[i,j,k,l] ≈ C[k,l,i,j]` on the first `dim` indices.
"""
function check_major_symmetry(C; atol=1e-10, dim=3)
  for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
    isapprox(C[i, j, k, l], C[k, l, i, j]; atol=atol) ||
      error("C not major-symmetric at ($i,$j,$k,$l)")
  end
  true
end

"""
    KirchhoffLovePlate <: Structure

Generic Kirchhoff-Love plate parameters for 1D or 2D structural manifolds
embedded in 2D or 3D fluids.

The interior-penalty C/DG discretisation follows Colomes, Verdugo &
Akkerman (2023), NME, §3.2, eqs. (24)–(25).
DOI: 10.1002/nme.7140

All weak forms are normalised by the ambient fluid density `ρ` so that
the assembled system matrices are dimensionally consistent with the
`PotentialFlow` and `FreeSurface` entities.

# Fields
- `E::Float64`            — Young's modulus [Pa]
- `ν::Float64`            — Poisson's ratio [dimensionless]
- `hb::Float64`           — Plate thickness [m]
- `ρ::Float64`            — Ambient fluid density [kg/m³]; default 1025.0
- `ρb::Float64`           — Plate material density [kg/m³]; default 256.25
- `g::Float64`            — Gravitational acceleration [m/s²]; default 9.81
- `ambient_dim::Int`      — Embedding-space dimension: 2 or 3; default 3
- `manifold_dim::Int`     — Plate manifold dimension: 1 (beam-like) or 2; default 2
- `symbol::Symbol`        — Field unknown symbol; default `:η`
- `space_domain_symbol::Symbol` — Triangulation key used for FE spaces; default `:Γη`
- `fe::FESpaceConfig`     — FE discretisation parameters
- `C`                     — Constitutive tensor `SymFourthOrderTensorValue{ambient_dim}`,
                            computed automatically via [`build_kl_tensor`](@ref).
                            The scalar `C[1,1,1,1] = D/ρ` where `D = E·h³/(12(1-ν²))`.

# Notes
- For `manifold_dim == 1` the model reduces to Euler-Bernoulli curvature
  rigidity; set `ambient_dim = 2` for 2D problems.
- `ambient_dim = 3, manifold_dim = 2` is the canonical 3D floating plate.
- The rotational penalty coefficient `γ` in `fe.γ` should be set to
  `O(p²)` (the default `10 * p^2` in `FESpaceConfig` is recommended).

# Example
```julia
plate = KirchhoffLovePlate(
    E  = 1.19e10,   # Pa  (e.g. ice)
    ν  = 0.13,
    hb = 2.0,        # m
    ρb = 256.25,     # kg/m³
)
```

See also: [`build_kl_tensor`](@ref), [`equivalent_beam_rigidity`](@ref)

# References
- [C23] Colomes, O., Verdugo, F., & Akkerman, I. (2023). A monolithic
  finite element formulation for the hydroelastic analysis of very large
  floating structures. *Int. J. Numer. Methods Eng.*, 124(3), 714-751.
  DOI: https://doi.org/10.1002/nme.7140
- `build_kl_tensor`: [C23] Section 3.2, Eq. (22)-(23)
"""
@with_kw struct KirchhoffLovePlate <: Structure
  E::Float64
  ν::Float64
  hb::Float64
  ρ::Float64 = 1025.0
  ρb::Float64 = 256.25
  g::Float64 = 9.81
  ambient_dim::Int = 3
  manifold_dim::Int = 2
  symbol::Symbol = :η
  space_domain_symbol::Symbol = :Γη
  fe::FESpaceConfig = FESpaceConfig()
  C = build_kl_tensor(ambient_dim, manifold_dim, E, ν, hb, ρ)
end

function print_parameters(plate::KirchhoffLovePlate)
  @printf("\n[MSG] Kirchhoff-Love Plate:\n")
  @printf("[VAL] E  = %.6e Pa\n", plate.E)
  @printf("[VAL] ν  = %.4f\n", plate.ν)
  @printf("[VAL] hb = %.4f m\n", plate.hb)
  @printf("[VAL] ρ  = %.2f kg/m3\n", plate.ρ)
  @printf("[VAL] ρb = %.2f kg/m3\n", plate.ρb)
  @printf("[VAL] g  = %.4f m/s2\n", plate.g)
  @printf("[VAL] ambient_dim = %d\n", plate.ambient_dim)
  @printf("[VAL] manifold_dim = %d\n", plate.manifold_dim)
  println()
end

variable_symbol(plate::KirchhoffLovePlate) = plate.symbol
ambient_dimension(plate::KirchhoffLovePlate) = plate.ambient_dim
manifold_dimension(plate::KirchhoffLovePlate) = plate.manifold_dim

"""
    equivalent_beam_rigidity(plate::KirchhoffLovePlate; width=1.0)

Equivalent Euler-Bernoulli bending rigidity per fluid density, `EI/ρ`.
"""
function equivalent_beam_rigidity(plate::KirchhoffLovePlate; width=1.0)
  width * plate.C[1, 1, 1, 1]
end


# ── Single-variable weak forms: mass, stiffness, rhs ──────────────────────
#
# The Kirchhoff-Love plate is discretised with an interior-penalty C/DG
# approach (Engel et al. 2002; Colomes et al. 2023, §3.2, eqs. 24-25).
# No intrinsic structural damping is implemented (has_damping_form = false).
# Fluid-structure coupling damping is inherited via Structure <: PhysicsParameters
# through the PotentialFlow↔Structure pair in CouplingTerms.jl.
#
# The Hessian of the scalar deflection η is computed as ε(∇(η)), where
# ε denotes the symmetric gradient operator. For a scalar field, ε(∇(η)) is
# the symmetric Hessian and therefore coincides with ∇∇(η):
#   ε(∇(η))ᵢⱼ = ∂²η/∂xᵢ∂xⱼ
# This is fully compatible with C::SymFourthOrderTensorValue{D} via ⊙.
#
# Bilinear form as assembled on the plate midsurface dom[:dΓη] with interior
# skeleton dom[:dΛη]:
#
#   a(η,v) = ∫ ( ∇∇(v) ⊙ (C ⊙ ∇∇(η)) + g·v·η ) dΓη
#            + ∫ ( - jump(∇(v)) ⋅ mean((C ⊙ ∇∇(η)) ⋅ n)
#                 - mean((C ⊙ ∇∇(v)) ⋅ n) ⋅ jump(∇(η))
#                 + (γ/hₑ) D_ρ jump(∇(v)) ⋅ jump(∇(η)) ) dΛη
#
# where D_ρ = C[1,1,1,1] = D/ρ_fluid = E·h³/(12(1-ν²)·ρ),
#       γ   = plate.fe.γ (= 10·p² from FESpaceConfig),
#       hₑ  = dom[:h_η] (representative element size on Γη).
#
# Simply-supported plate BC on Γη: η=0 enforced as Dirichlet; M_n=0 is natural.
# ==========================================================================

"""
    has_damping_form(::KirchhoffLovePlate)

The plate has no structural damping term; returns `false`.
FSI coupling damping is provided separately via `CouplingTerms.jl`.
"""
has_damping_form(::KirchhoffLovePlate) = false

"""
    mass(plate, dom, x_tt, y)

Mass bilinear form: ∫ (ρb·hb/ρ) v·ηtt dΓ.
"""
function mass(s::KirchhoffLovePlate, dom::IntegrationDomains, x_tt, y)
  sym  = variable_symbol(s)
  ηₜₜ = x_tt[sym]
  v    = y[sym]
  m_ρ  = s.ρb * s.hb / s.ρ   # mass per unit area / ρ_fluid  [dimensionless]
  ∫(m_ρ * v * ηₜₜ)dom[:dΓη]
end

"""
    stiffness(plate, dom, x, y)

  C/DG bilinear form for the Kirchhoff-Love plate (eqs. 24-25 of [C23]).

Uses interior-penalty stabilisation with penalty coefficient
`(γ / h) * C[1,1,1,1]` where γ comes from `plate.fe.γ`.
"""
function stiffness(s::KirchhoffLovePlate, dom::IntegrationDomains, x, y)
  sym  = variable_symbol(s)
  η    = x[sym]
  v    = y[sym]
  D_ρ  = s.C[1, 1, 1, 1]   # bending stiffness / ρ_fluid = D/ρ [m⁴/s²]
  γ    = s.fe.γ
  h    = dom[:h_η]
  n_Λ  = dom[:n_Λ_η]

  # Kirchhoff-Love plate C/DG formulation.
  # Bulk term: ∫_Γb ∇∇v ⊙ (C ⊙ ∇∇η) dΓ + hydrostatic restoring term.
  # Skeleton terms: consistency + symmetry + penalty.
  # Reference: [C23] Section 3.2, Eq. (24)-(25).
  bulk = ∫(( ∇∇(v) ⊙ (s.C ⊙ ∇∇(η))) + s.g * v * η)dom[:dΓη]

  skeleton = ∫(
    -jump(∇(v)) ⊙ (mean(s.C ⊙ ∇∇(η))⋅n_Λ.⁺) 
    - (mean((s.C ⊙ ∇∇(v)))⋅n_Λ.⁺) ⊙ jump(∇(η)) 
    + D_ρ*γ/h*jump(∇(v))⊙jump(∇(η)) )dom[:dΛη]

  return bulk + skeleton
end

"""
    rhs(plate, dom, f, y)

Right-hand side linear form: ∫ v · f[sym] dΓ.
"""
function rhs(s::KirchhoffLovePlate, dom::IntegrationDomains, f, y)
  sym = variable_symbol(s)
  v   = y[sym]
  ∫(v * f[sym])dom[:dΓη]
end
