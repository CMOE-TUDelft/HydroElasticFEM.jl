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

Generic Kirchhoff-Love plate parameters supporting 1D/2D manifolds embedded
in 2D/3D fluids.
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
  space_domain_symbol::Symbol = :Γb
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

const KirchhoffLovePlate3D = KirchhoffLovePlate
