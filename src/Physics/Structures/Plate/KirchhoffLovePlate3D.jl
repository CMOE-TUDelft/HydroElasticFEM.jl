using Gridap.TensorValues

"""
    build_KL_tensor(E, ν, hb, ρ)

Build the 3D Kirchhoff-Love constitutive tensor used in the Yago benchmark.

Only in-plane indices `i,j,k,l ∈ {1,2}` are populated (plane-stress
Kirchhoff-Love assumption). Components involving index 3 remain zero.

# Arguments
- `E::Float64`  : Young's modulus [Pa]
- `ν::Float64`  : Poisson ratio [-]
- `hb::Float64` : plate thickness [m]
- `ρ::Float64`  : water density [kg/m^3]
"""
function build_KL_tensor(E, ν, hb, ρ)
  I = hb^3 / 12
  μ = E / (2 * (1 + ν))
  λ = ν * E / (1 - ν^2)
  δ(x, y) = ==(x, y)

  C_type = SymFourthOrderTensorValue{3, Float64}
  Cvals = zero(Array{Float64}(undef, 36))

  for i in 1:2, j in 1:2, k in 1:2, l in 1:2
    Cvals[data_index(C_type, i, j, k, l)] =
      I / ρ * (
        μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k)) +
        λ * (δ(i, j) * δ(k, l))
      )
  end

  return SymFourthOrderTensorValue(Cvals...)
end

"""
    check_major_symmetry(C; atol=1e-10)

Verify major symmetry `C[i,j,k,l] ≈ C[k,l,i,j]` for a 3D fourth-order
constitutive tensor.
"""
function check_major_symmetry(C; atol=1e-10)
  for i in 1:3, j in 1:3, k in 1:3, l in 1:3
    isapprox(C[i, j, k, l], C[k, l, i, j]; atol=atol) ||
      error("C not major-symmetric at ($i,$j,$k,$l)")
  end
  return true
end

"""
    KirchhoffLovePlate3D <: Structure

Parameter container for a 3D Kirchhoff-Love plate entity.

This type stores constitutive and numerical parameters used by the 3D Yago
benchmark. The full monolithic weak form is assembled in the benchmark
example module.
"""
@with_kw struct KirchhoffLovePlate3D <: Structure
  E::Float64
  ν::Float64
  hb::Float64
  ρ::Float64 = 1025.0
  ρb::Float64 = 256.25
  g::Float64 = 9.81
  symbol::Symbol = :η
  space_domain_symbol::Symbol = :Γb
  fe::FESpaceConfig = FESpaceConfig()
  C = build_KL_tensor(E, ν, hb, ρ)
end

function print_parameters(plate::KirchhoffLovePlate3D)
  @printf("\n[MSG] Kirchhoff-Love Plate (3D):\n")
  @printf("[VAL] E  = %.6e Pa\n", plate.E)
  @printf("[VAL] ν  = %.4f\n", plate.ν)
  @printf("[VAL] hb = %.4f m\n", plate.hb)
  @printf("[VAL] ρ  = %.2f kg/m3\n", plate.ρ)
  @printf("[VAL] ρb = %.2f kg/m3\n", plate.ρb)
  @printf("[VAL] g  = %.4f m/s2\n", plate.g)
  println()
end

variable_symbol(plate::KirchhoffLovePlate3D) = plate.symbol
