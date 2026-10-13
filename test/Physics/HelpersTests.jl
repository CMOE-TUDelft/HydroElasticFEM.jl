using Test
using Gridap
import HydroElasticFEM.Physics as P

@testset "_resolve_space_function input shapes" begin
  x = Point(1.0, 2.0)
  t = 0.5

  # constants
  @test P._resolve_space_function(3.0, t)(x) == 3.0
  @test P._resolve_space_function(3.0, nothing)(x) == 3.0

  # space function x -> ...
  f = x -> x[1] + 2x[2]
  @test P._resolve_space_function(f, t)(x) == 5.0
  @test P._resolve_space_function(f, nothing)(x) == 5.0

  # space function whose call with a scalar time does not error (x[1] of a Float64)
  g = x -> 7 * x[1]
  @test P._resolve_space_function(g, t)(x) == 7.0

  # time-indexed space function t -> (x -> ...)
  h = t -> (x -> x[2] * t)
  @test P._resolve_space_function(h, t)(x) == 1.0

  # space-time function (x, t) -> ...
  k = (x, t) -> x[1] + t
  @test P._resolve_space_function(k, t)(x) == 1.5

  # generic function with both (x) and (x, t) methods: (x, t) wins, as before
  both(x) = -1.0
  both(x, t) = x[1] * t
  @test P._resolve_space_function(both, t)(x) == 0.5

  # errors inside the user function propagate at evaluation time
  bad = (x, t) -> error("boom")
  @test_throws ErrorException P._resolve_space_function(bad, t)(x)
end
