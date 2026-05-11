using Test

@testset "HydroElasticFEM" begin

  # ==========================================================================
  # Geometry
  # ==========================================================================
  include("Geometry/GeometryTests.jl")

  # ==========================================================================
  # Physics
  # ==========================================================================

  include("Physics/PhysicsTests.jl")

  # ==========================================================================
  # Simulation
  # ==========================================================================

  include("Simulation/SimulationTests.jl")

  # ==========================================================================
  # Utilities
  # ==========================================================================

  include("UtilitiesTests.jl")

  # ==========================================================================
  # Time-domain solvers
  # ==========================================================================

  include("Membrane/TimeDomain/MemTimeLrmmTests.jl")
  include("Membrane/TimeDomain/MemTimeDampingZonePackageTests.jl")

  # ==========================================================================
  # Examples (integration tests against analytical solutions)
  # ==========================================================================

  include("examples/EulerBernoulliBeamWeakFormTests.jl")
  include("examples/KhabakhpashevaBeamJointTests.jl")
  include("examples/test_liu_benchmark.jl")
  include("examples/test_yago_3d_freq.jl")

end
