using Test
using Gridap
import HydroElasticFEM.Physics as P

@testset "Physics - module interface" begin

  # Abstract base type exists
  @test P.PhysicsParameters isa Type

  # Default print_parameters throws for unknown subtypes
  struct _TestParams <: P.PhysicsParameters end
  @test_throws ErrorException P.print_parameters(_TestParams())

  # print_parameters works for concrete types
  ρw = 1025.0
  mem = P.Membrane2D(L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1, τ=0.0)
  @test_nowarn P.print_parameters(mem)
  beam = P.EulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=500e6*6.667e-4/ρw, τ=0.0)
  @test_nowarn P.print_parameters(beam)
  resn = P.ResonatorSingle(M=1e3, K=5.9e3, C=0.0, XZ=VectorValue(10.0, 0.0))
  @test_nowarn P.print_parameters(resn)
end

@testset "Physics - Physics modules" begin

  @testset "PotentialFlow" include("PotentialFlowTests.jl")
  @testset "Membrane2D" include("Membrane2DTests.jl")
  @testset "EulerBernoulliBeam" include("EulerBernoulliBeamTests.jl")
  @testset "Resonator" include("ResonatorTests.jl")
  @testset "FreeSurface" include("FreeSurfaceTests.jl")

end
