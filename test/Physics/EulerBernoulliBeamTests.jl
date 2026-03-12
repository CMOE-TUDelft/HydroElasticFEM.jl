using Test
import HydroElasticFEM.Physics as P

@testset "EulerBernoulliBeam struct" begin
  ρw = 1025.0
  EIᵨ = 500e6 * 6.667e-4 / ρw
  beam = P.EulerBernoulliBeam(
    L=20.0, mᵨ=192.956/ρw, EIᵨ=EIᵨ, τ=0.0)
  @test beam.EIᵨ ≈ EIᵨ
  @test beam.ωn1 ≈ 22.3733 * sqrt(EIᵨ / ((192.956/ρw) * 20.0^4))
  @test beam isa P.PhysicsParameters

  # Defaults: τ and bndType default, derived fields auto-computed
  beam_def = P.EulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=EIᵨ)
  @test beam_def.τ == 0.0
  @test beam_def.EIᵨ ≈ EIᵨ
  @test beam_def.ωn1 ≈ 22.3733 * sqrt(beam_def.EIᵨ / ((192.956/ρw) * 20.0^4))

  # Nonzero τ
  beam_d = P.EulerBernoulliBeam(L=10.0, mᵨ=100.0/ρw, EIᵨ=1e9*1e-3/ρw, τ=0.05)
  @test beam_d.τ == 0.05
end
