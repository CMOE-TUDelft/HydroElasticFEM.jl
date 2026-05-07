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
  @test_throws ErrorException G.get_boundary(d2, "surface")

  trians = G.build_triangulations(d2, G.build_model(d2))
  @test num_cells(trians[:Ω]) == 20 * 4
  @test num_cells(trians[:Γfs]) == 20
  @test num_cells(trians[:Γin]) == 4
  @test num_cells(trians[:Γout]) == 4
  @test num_cells(trians[:Γη]) == 0
  @test !haskey(trians, :Γlateral)
end

@testset "TankDomain generic 2D constructor" begin
  tank = G.TankDomain(
    L = 10.0,
    H = 2.0,
    nx = 20,
    ny = 4,
    structure_domains = [G.StructureDomain(L = 1.0, x₀ = [2.0, 2.0])],
  )

  @test tank isa G.TankDomain{2}
  @test tank.L == 10.0
  @test tank.H == 2.0
  @test tank.nx == 20
  @test tank.ny == 4
  @test length(tank.structure_domains) == 1

  Γη = G.get_boundary(tank, "structure")
  @test num_cells(Γη) > 0
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

  trians = G.build_triangulations(d3, G.build_model(d3))
  @test num_cells(trians[:Ω]) == 8 * 3 * 2
  @test num_cells(trians[:Γfs]) == 8 * 3
  @test num_cells(trians[:Γin]) == 3 * 2
  @test num_cells(trians[:Γout]) == 3 * 2
  @test num_cells(trians[:Γbot]) == 8 * 3
  @test num_cells(trians[:Γlateral]) == 2 * 8 * 2
  @test num_cells(trians[:Γη]) == 0
end

@testset "TankDomain generic 3D constructor" begin
  tank = G.TankDomain(L = 8.0, W = 3.0, H = 2.0, nx = 8, ny = 3, nz = 2)

  @test tank isa G.TankDomain{3}
  @test tank.L == 8.0
  @test tank.W == 3.0
  @test tank.H == 2.0
  @test tank.nx == 8
  @test tank.ny == 3
  @test tank.nz == 2
  @test isempty(tank.structure_domains)
  @test isempty(tank.damping_zones)
  @test isempty(tank.joint_domains)
end

@testset "TankDomain{3} delegates to generic Cartesian path" begin
  tank = G.TankDomain(L = 8.0, W = 3.0, H = 2.0, nx = 8, ny = 3, nz = 2)

  @test G.ambient_dimension(tank) == 3
  @test G.manifold_dimension(tank) == 3
  @test G.boundary_tags(tank) == G.boundary_tags(
    G.CartesianDomain(L = 8.0, W = 3.0, H = 2.0, nx = 8, ny = 3, nz = 2),
  )

  Γfs = G.get_boundary(tank, "free_surface")
  Γlat = G.get_boundary(tank, "lateral_walls")
  @test num_cells(Γfs) == 8 * 3
  @test num_cells(Γlat) == 2 * 8 * 2

  trians = G.build_triangulations(tank, G.build_model(tank))
  @test num_cells(trians[:Ω]) == 8 * 3 * 2
  @test num_cells(trians[:Γfs]) == 8 * 3
  @test num_cells(trians[:Γlateral]) == 2 * 8 * 2
  @test num_cells(trians[:Γη]) == 0
end

@testset "TankDomain{3} rejects unsupported structured subdomains" begin
  @test_throws ErrorException G.TankDomain(
    structure_domains = [G.StructureDomain(L = 1.0, x₀ = [1.0, 1.0])],
    W = 2.0,
    nz = 2,
  )
  @test_throws ErrorException G.TankDomain(
    L = 8.0,
    W = 3.0,
    H = 2.0,
    nx = 8,
    ny = 3,
    nz = 2,
    damping_zones = [G.DampingZone(L = 0.5, x₀ = [0.0, 2.0])],
  )
end

@testset "CartesianDomain3D benchmark wrapper uses structured boundaries" begin
  d = G.CartesianDomain3D(LΩ = 8.0, BΩ = 3.0, H = 2.0, nx_total = 8, ny_total = 3, nz = 2)

  @test G.ambient_dimension(d) == 3
  @test G.manifold_dimension(d) == 3
  @test haskey(G.boundary_tags(d), "lateral_walls")

  @test num_cells(G.get_boundary(d, "surface")) == 8 * 3
  @test num_cells(G.get_boundary(d, "free_surface")) == 8 * 3
  @test num_cells(G.get_boundary(d, "lateral_walls")) == 2 * 8 * 2

  trians = G.build_triangulations(d, G.build_model(d))
  @test num_cells(trians[:Ω]) == 8 * 3 * 2
  @test num_cells(trians[:Γfs]) == 8 * 3
  @test num_cells(trians[:Γlateral]) == 2 * 8 * 2
end
