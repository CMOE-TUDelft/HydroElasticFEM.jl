@testset "Membrane2D" begin
  mem = HydroElasticFEM.Membrane2D(
    L=20.0, m=922.5, T=98.1 * 1025.0, τ=0.0, bndType=FreeBoundary())
  @test mem.L == 20.0
  @test mem.m == 922.5
  @test mem.T == 98.1 * 1025.0
  @test mem.τ == 0.0
  @test mem.MTotal ≈ 922.5 * 20.0
  @test mem.ωn1 ≈ (π / 20.0) * sqrt(98.1 * 1025.0 / 922.5)
  @test mem isa HydroElasticFEM.AbstractStructure

  # Fixed boundary
  mem_fix = HydroElasticFEM.Membrane2D(
    L=10.0, m=500.0, T=1000.0, τ=0.1, bndType=FixedBoundary())
  @test mem_fix.bndType isa FixedBoundary

  # Defaults: τ and bndType default, derived fields auto-computed
  mem_def = Membrane2D(L=20.0, m=922.5, T=98.1 * 1025.0)
  @test mem_def.τ == 0.0
  @test mem_def.bndType isa FreeBoundary
  @test mem_def.MTotal ≈ 922.5 * 20.0
  @test mem_def.ωn1 ≈ (π / 20.0) * sqrt(98.1 * 1025.0 / 922.5)
end
