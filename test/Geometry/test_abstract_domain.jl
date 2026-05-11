using Test
using Gridap
import HydroElasticFEM.Geometry as G

# ─────────────────────────────────────────────────────────────
# AbstractDomain interface — TankDomain{2}
# ─────────────────────────────────────────────────────────────

@testset "AbstractDomain — STANDARD_TAGS" begin
  @test length(G.STANDARD_TAGS) == 6
  for tag in ["fluid", "free_surface", "seabed", "inlet", "outlet", "structure"]
    @test tag in G.STANDARD_TAGS
  end
end

@testset "AbstractDomain — TankDomain{2} subtype" begin
  @test G.TankDomain{2} <: G.AbstractDomain
  tank = G.TankDomain(L=4.0, H=1.0, nx=20, ny=2)
  @test tank isa G.AbstractDomain
end

@testset "ambient_dimension / manifold_dimension — TankDomain{2}" begin
  tank = G.TankDomain()
  @test G.ambient_dimension(tank) == 2
  @test G.manifold_dimension(tank) == 2
end

@testset "boundary_tags — TankDomain{2}" begin
  tank  = G.TankDomain()
  btags = G.boundary_tags(tank)
  @test btags isa Dict{String, String}
  for tag in G.STANDARD_TAGS
    @test haskey(btags, tag)
  end
  # Cartesian label mapping sanity checks
  @test btags["free_surface"] == "surface"
  @test btags["seabed"]       == "bottom"
  @test btags["inlet"]        == "inlet"
  @test btags["outlet"]       == "outlet"
  @test btags["fluid"]        == "water"
end

@testset "get_boundary — TankDomain{2} non-structure tags" begin
  tank = G.TankDomain(L=4.0, H=1.0, nx=20, ny=4)
  for name in ["free_surface", "seabed", "inlet", "outlet"]
    bnd = G.get_boundary(tank, name)
    @test bnd isa Gridap.Geometry.BoundaryTriangulation
    @test num_cells(bnd) > 0
  end
end

@testset "get_boundary — TankDomain{2} fluid (interior)" begin
  tank = G.TankDomain(L=4.0, H=1.0, nx=20, ny=4)
  Ω = G.get_boundary(tank, "fluid")
  # Returns interior — should have cells equal to nx*ny
  @test num_cells(Ω) == 20 * 4
end

@testset "get_boundary — TankDomain{2} structure (empty when no structures)" begin
  tank  = G.TankDomain(L=4.0, H=1.0, nx=20, ny=4)
  Γη    = G.get_boundary(tank, "structure")
  @test num_cells(Γη) == 0
end

@testset "get_boundary — TankDomain{2} structure (with structures)" begin
  s1   = G.StructureDomain(L=1.0, x₀=[1.5, 1.0])
  tank = G.TankDomain(L=4.0, H=1.0, nx=40, ny=4, structure_domains=[s1])
  Γη   = G.get_boundary(tank, "structure")
  # 1.0 m structure / 0.1 m element width = 10 cells
  @test num_cells(Γη) == 10
end

@testset "get_boundary — TankDomain{2} invalid tag raises error" begin
  tank = G.TankDomain()
  @test_throws ErrorException G.get_boundary(tank, "nonexistent_tag")
end

@testset "triangulation — TankDomain{2}" begin
  tank = G.TankDomain(L=4.0, H=1.0, nx=10, ny=2)
  Ω    = G.triangulation(tank)
  @test Ω isa Gridap.Geometry.BodyFittedTriangulation
  @test num_cells(Ω) == 10 * 2
end

@testset "AbstractDomain — interface error on unimplemented subtype" begin
  # Create a dummy subtype that does not implement any method
  struct _DummyDomain <: G.AbstractDomain end
  d = _DummyDomain()
  @test_throws ErrorException G.triangulation(d)
  @test_throws ErrorException G.boundary_tags(d)
  @test_throws ErrorException G.ambient_dimension(d)
  @test_throws ErrorException G.manifold_dimension(d)
  @test_throws ErrorException G.get_boundary(d, "fluid")
end
