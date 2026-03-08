@testset "PhysicalEntities - module interface" begin
  PE = HydroElasticFEM.PhysicalEntities

  # Abstract base type exists
  @test PE.PhysicsParameters isa Type

  # Default print_parameters throws for unknown subtypes
  struct _TestParams <: PE.PhysicsParameters end
  @test_throws ErrorException PE.print_parameters(_TestParams())

  # print_parameters works for concrete types
  ρw = 1025.0
  mem = PE.Membrane2D(L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1, τ=0.0, bndType=PE.FreeBoundary())
  @test_nowarn PE.print_parameters(mem)
  beam = PE.EulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=500e6*6.667e-4/ρw, τ=0.0, bndType=PE.FreeBoundary())
  @test_nowarn PE.print_parameters(beam)
  resn = PE.ResonatorSingle(M=1e3, K=5.9e3, C=0.0, XZ=VectorValue(10.0, 0.0))
  @test_nowarn PE.print_parameters(resn)
end
