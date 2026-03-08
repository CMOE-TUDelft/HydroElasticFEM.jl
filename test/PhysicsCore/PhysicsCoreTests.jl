@testset "PhysicsCore" begin
  include("PhysicalEntitiesTests.jl")
  include("FreeSurfaceTests.jl")
  include("Membrane2DTests.jl")
  include("EulerBernoulliBeamTests.jl")
  include("ResonatorTests.jl")
  include("WeakFormAssemblyTests.jl")
  include("FESpaceAssemblyTests.jl")
end
