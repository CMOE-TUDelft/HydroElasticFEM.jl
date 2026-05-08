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
