using Test
using Gridap
import HydroElasticFEM.PhysicsCore.Entities as E

@testset "Entities - module interface" begin

  # Abstract base type exists
  @test E.PhysicsParameters isa Type

  # Default print_parameters throws for unknown subtypes
  struct _TestParams <: E.PhysicsParameters end
  @test_throws ErrorException E.print_parameters(_TestParams())

  # print_parameters works for concrete types
  ρw = 1025.0
  mem = E.Membrane2D(L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1, τ=0.0)
  @test_nowarn E.print_parameters(mem)
  beam = E.EulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=500e6*6.667e-4/ρw, τ=0.0)
  @test_nowarn E.print_parameters(beam)
  resn = E.ResonatorSingle(M=1e3, K=5.9e3, C=0.0, XZ=VectorValue(10.0, 0.0))
  @test_nowarn E.print_parameters(resn)
end

@testset "Entities - Physics modules" begin

  @testset "PotentialFlow" include("PotentialFlowTests.jl")
  @testset "Membrane2D" include("Membrane2DTests.jl")
  @testset "EulerBernoulliBeam" include("EulerBernoulliBeamTests.jl")
  @testset "Resonator" include("ResonatorTests.jl")
  @testset "FreeSurface" include("FreeSurfaceTests.jl")

end
