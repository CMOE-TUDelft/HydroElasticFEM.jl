using Test

@testset "Geometry" begin

  # Run GmshDomain tests first: gmsh's C++ allocator is sensitive to heap
  # state after extensive Gridap use on macOS ARM64.
  include("test_gmsh_domain.jl")

  include("CartesianGeometryTests.jl")
  include("test_abstract_domain.jl")

end