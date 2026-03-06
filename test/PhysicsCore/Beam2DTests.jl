@testset "Beam2D" begin
  beam = HydroElasticFEM.Beam2D(
    L=20.0, m=192.956, E=500e6, I=6.667e-4, τ=0.0, bndType=FreeBoundary())
  @test beam.EI ≈ 500e6 * 6.667e-4
  @test beam.τEI ≈ 0.0
  @test beam.ωn1 ≈ 22.3733 * sqrt(beam.EI / (192.956 * 20.0^4))
  @test beam isa HydroElasticFEM.AbstractStructure

  # Defaults: τ and bndType default, derived fields auto-computed
  beam_def = Beam2D(L=20.0, m=192.956, E=500e6, I=6.667e-4)
  @test beam_def.τ == 0.0
  @test beam_def.bndType isa FreeBoundary
  @test beam_def.EI ≈ 500e6 * 6.667e-4
  @test beam_def.τEI ≈ 0.0
  @test beam_def.MTotal ≈ 192.956 * 20.0
  @test beam_def.ωn1 ≈ 22.3733 * sqrt(beam_def.EI / (192.956 * 20.0^4))

  # Nonzero τ
  beam_d = Beam2D(L=10.0, m=100.0, E=1e9, I=1e-3, τ=0.05)
  @test beam_d.τEI ≈ 0.05 * 1e9 * 1e-3
end
