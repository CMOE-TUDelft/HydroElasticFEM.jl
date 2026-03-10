using Test

@testset "PhysicsCore" begin

  @testset "Entities module tests" include("Entities/EntitiesTests.jl")

  include("FEOperatorsTests.jl")
  include("FESpaceAssemblyTests.jl")
end
