using Test
using Gridap
import HydroElasticFEM.Geometry as G

@testset "CartesianDomain generic 2D" begin
  d2 = G.CartesianDomain(L = 10.0, H = 2.0, nx = 20, ny = 4)

  @test G.ambient_dimension(d2) == 2
  @test G.manifold_dimension(d2) == 2

  Ω = G.triangulation(d2)
  @test num_cells(Ω) == 20 * 4

  Γfs = G.get_boundary(d2, "free_surface")
  Γin = G.get_boundary(d2, "inlet")
  Γout = G.get_boundary(d2, "outlet")

  @test num_cells(Γfs) == 20
  @test num_cells(Γin) == 4
  @test num_cells(Γout) == 4
end

@testset "CartesianDomain generic 3D" begin
  d3 = G.CartesianDomain(L = 8.0, W = 3.0, H = 2.0, nx = 8, ny = 3, nz = 2)

  @test G.ambient_dimension(d3) == 3
  @test G.manifold_dimension(d3) == 3

  Ω = G.triangulation(d3)
  @test num_cells(Ω) == 8 * 3 * 2

  Γfs = G.get_boundary(d3, "free_surface")
  Γlat = G.get_boundary(d3, "lateral_walls")

  @test num_cells(Γfs) == 8 * 3
  @test num_cells(Γlat) == 2 * 8 * 2
end
