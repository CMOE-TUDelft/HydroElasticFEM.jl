@testset "Membrane2D" begin
  ρw = 1025.0
  mem = HydroElasticFEM.Membrane2D(
    L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1, τ=0.0, bndType=FreeBoundary())
  @test mem.L == 20.0
  @test mem.mᵨ ≈ 922.5 / ρw
  @test mem.Tᵨ == 98.1
  @test mem.τ == 0.0
  @test mem.ωn1 ≈ (π / 20.0) * sqrt(98.1 / (922.5 / ρw))
  @test mem isa HydroElasticFEM.AbstractStructure

  # Fixed boundary
  mem_fix = HydroElasticFEM.Membrane2D(
    L=10.0, mᵨ=500.0/ρw, Tᵨ=1000.0/ρw, τ=0.1, bndType=FixedBoundary())
  @test mem_fix.bndType isa FixedBoundary

  # Defaults: τ and bndType default, derived fields auto-computed
  mem_def = Membrane2D(L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1)
  @test mem_def.τ == 0.0
  @test mem_def.bndType isa FreeBoundary
  @test mem_def.ωn1 ≈ (π / 20.0) * sqrt(98.1 / (922.5 / ρw))
end
