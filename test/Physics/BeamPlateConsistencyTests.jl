using Test
using Gridap

import HydroElasticFEM.Physics as P

@testset "Beam-plate 1D consistency" begin
  ρw = 1025.0
  E = 11.9e9
  ν = 0.13
  hb = 2.0
  L = 20.0

  EIρ = E * hb^3 / (12 * ρw)

  beam = P.EulerBernoulliBeam(
    L = L,
    mᵨ = 1.0,
    EIᵨ = EIρ,
  )
  plate = P.KirchhoffLovePlate(
    E = E,
    ν = ν,
    hb = hb,
    ρ = ρw,
    ambient_dim = 2,
    manifold_dim = 1,
  )

  @test P.equivalent_beam_rigidity(plate; width = 1.0) ≈ beam.EIᵨ atol = 1e-10
  @test plate.C[1, 1, 1, 1] ≈ beam.EIᵨ atol = 1e-10
  @test plate.C[1, 1, 2, 2] ≈ 0.0 atol = 1e-12
end

@testset "KirchhoffLovePlate 2D tensor — major symmetry" begin
  # 2D plane-stress tensor: ambient_dim=2, manifold_dim=2
  E  = 11.9e9
  ν  = 0.13
  hb = 2.0
  ρw = 1025.0
  C2 = P.build_kl_tensor(2, 2, E, ν, hb, ρw)

  @test P.check_major_symmetry(C2; atol = 1e-10, dim = 2)

  # Minor symmetry
  for i in 1:2, j in 1:2, k in 1:2, l in 1:2
    @test C2[i, j, k, l] ≈ C2[j, i, k, l]  atol = 1e-12
    @test C2[i, j, k, l] ≈ C2[i, j, l, k]  atol = 1e-12
  end

  # Diagonal bending rigidity C[1,1,1,1] = D/ρ
  D_rho = E * hb^3 / (12 * (1 - ν^2) * ρw)
  @test C2[1, 1, 1, 1] ≈ D_rho  atol = 1e-8
  @test C2[2, 2, 2, 2] ≈ D_rho  atol = 1e-8

  # Cross-term: C[1,1,2,2] = ν·D/ρ
  @test C2[1, 1, 2, 2] ≈ ν * D_rho  atol = 1e-8
end
