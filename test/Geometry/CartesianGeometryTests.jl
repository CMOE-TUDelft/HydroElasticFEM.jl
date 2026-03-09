using Test
using Gridap
import HydroElasticFEM.Geometry as G

@testset "TankDomain2D" begin

  structure_a = G.StructureDomain1D(L=1.0, x₀=[0.5,1.0])
  structure_b = G.StructureDomain1D(L=2.0, x₀=[0.5,1.0])
  damping_a = G.DampingZone1D(L=0.5, x₀=[0.0,1.0])
  damping_b = G.DampingZone1D(L=0.5, x₀=[3.5,1.0])
  tank = G.TankDomain2D(L=10.0, H=1.0, nx=20, ny=2, structure_domains=[structure_a, structure_b], damping_zones=[damping_a, damping_b])

  @test tank.L == 10.0
  @test tank.H == 1.0
  @test tank.nx == 20
  @test tank.ny == 2
  @test length(tank.structure_domains) == 2
  @test length(tank.damping_zones) == 2
  @test tank.structure_domains[1].L == 1.0
  @test tank.structure_domains[2].x₀ == [0.5, 1.0]
  @test tank.damping_zones[1].σ == 1.0
  @test tank.damping_zones[2].x₀ == [3.5, 1.0]

end

@testset "build_model" begin
  tank = G.TankDomain2D()
  model = G.build_model(tank)
  @test model isa Gridap.Geometry.CartesianDiscreteModel
end

@testset "surface_mask - StructureDomain1D" begin
  # Structure at x ∈ [1.0, 3.0], y = 1.0
  sd = G.StructureDomain1D(L=2.0, x₀=[1.0, 1.0])
  mask = G.surface_mask(sd)

  # Centroid inside → true
  @test mask([VectorValue(1.5, 1.0), VectorValue(2.5, 1.0)]) == true
  # Centroid outside (left of zone) → false
  @test mask([VectorValue(0.0, 1.0), VectorValue(0.5, 1.0)]) == false
  # Centroid outside (wrong y) → false
  @test mask([VectorValue(1.5, 0.5), VectorValue(2.5, 0.5)]) == false
  # Right at boundary → true
  @test mask([VectorValue(1.0, 1.0), VectorValue(1.0, 1.0)]) == true
  @test mask([VectorValue(3.0, 1.0), VectorValue(3.0, 1.0)]) == true
end

@testset "surface_mask - DampingZone1D" begin
  # Damping zone at x ∈ [0.0, 0.5], y = 1.0
  dz = G.DampingZone1D(L=0.5, x₀=[0.0, 1.0])
  mask = G.surface_mask(dz)

  @test mask([VectorValue(0.1, 1.0), VectorValue(0.3, 1.0)]) == true
  @test mask([VectorValue(0.6, 1.0), VectorValue(0.8, 1.0)]) == false
end

@testset "surface_masks - full domain" begin
  s1 = G.StructureDomain1D(L=1.0, x₀=[2.0, 1.0])
  s2 = G.StructureDomain1D(L=1.5, x₀=[5.0, 1.0])
  d1 = G.DampingZone1D(L=0.5, x₀=[0.0, 1.0])
  d2 = G.DampingZone1D(L=0.5, x₀=[9.5, 1.0])
  tank = G.TankDomain2D(L=10.0, H=1.0, nx=20, ny=2,
    structure_domains=[s1, s2], damping_zones=[d1, d2])

  smasks, dmasks = G.surface_masks(tank)

  @test length(smasks) == 2
  @test length(dmasks) == 2

  # Structure mask 1 hits its zone
  @test smasks[1]([VectorValue(2.5, 1.0)]) == true
  @test smasks[1]([VectorValue(5.5, 1.0)]) == false

  # Structure mask 2 hits its zone
  @test smasks[2]([VectorValue(5.5, 1.0)]) == true
  @test smasks[2]([VectorValue(2.5, 1.0)]) == false

  # Damping masks
  @test dmasks[1]([VectorValue(0.2, 1.0)]) == true
  @test dmasks[2]([VectorValue(9.7, 1.0)]) == true
  @test dmasks[1]([VectorValue(9.7, 1.0)]) == false
end

@testset "build_triangulations — single structure + two damping zones" begin
  # Domain: [0, 4] × [0, 1],  40 × 4 elements
  # Structure at x ∈ [1.5, 2.5] on y = 1.0  (L = 1.0)
  # Damping inlet  x ∈ [0.0, 0.5]  on y = 1.0
  # Damping outlet x ∈ [3.5, 4.0]  on y = 1.0
  s1 = G.StructureDomain1D(L=1.0, x₀=[1.5, 1.0])
  d1 = G.DampingZone1D(L=0.5, x₀=[0.0, 1.0])
  d2 = G.DampingZone1D(L=0.5, x₀=[3.5, 1.0])
  tank = G.TankDomain2D(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains=[s1], damping_zones=[d1, d2])

  model = G.build_model(tank)
  tri   = G.build_triangulations(tank, model)

  @test tri isa G.TankTriangulations

  # Number of per-zone triangulations matches input
  @test length(tri.Γ_structures) == 1
  @test length(tri.Γ_dampings)   == 2

  # Surface partition: every surface cell is in exactly one of
  # Γ_structures[i], Γ_dampings[j], or Γfs.
  # Total surface cells = nx = 40
  n_struct = sum(num_cells(t) for t in tri.Γ_structures)
  n_damp   = sum(num_cells(t) for t in tri.Γ_dampings)
  n_fs     = num_cells(tri.Γfs)
  @test n_struct + n_damp + n_fs == num_cells(tri.Γ)

  # Γκ = non-structure = damping + free surface
  @test num_cells(tri.Γκ) == n_damp + n_fs
  # Γη = all structures
  @test num_cells(tri.Γη) == n_struct

  # Sanity: structure zone gets 10 cells (1.0 / 0.1 element size)
  @test num_cells(tri.Γ_structures[1]) == 10
  # Each damping zone gets 5 cells (0.5 / 0.1)
  @test num_cells(tri.Γ_dampings[1]) == 5
  @test num_cells(tri.Γ_dampings[2]) == 5
end

@testset "build_triangulations — no structures or damping" begin
  tank = G.TankDomain2D(L=4.0, H=1.0, nx=20, ny=2)
  model = G.build_model(tank)
  tri   = G.build_triangulations(tank, model)

  @test length(tri.Γ_structures) == 0
  @test length(tri.Γ_dampings)   == 0
  # Entire surface is free surface
  @test num_cells(tri.Γfs) == num_cells(tri.Γ)
  @test num_cells(tri.Γκ)  == num_cells(tri.Γ)
  @test num_cells(tri.Γη)  == 0
end

@testset "get_weak_form_domains" begin
  s1 = G.StructureDomain1D(L=1.0, x₀=[1.5, 1.0])
  d1 = G.DampingZone1D(L=0.5, x₀=[0.0, 1.0])
  d2 = G.DampingZone1D(L=0.5, x₀=[3.5, 1.0])
  tank = G.TankDomain2D(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains=[s1], damping_zones=[d1, d2])
  model = G.build_model(tank)
  trians   = G.build_triangulations(tank, model)

  d = G.get_weak_form_domains(trians; degree=4)

  @test d isa Dict{Symbol, Any}

  # Core keys present
  @test haskey(d, :dΩ)
  @test haskey(d, :dΓ_fs)
  @test haskey(d, :dΓ_s)
  @test haskey(d, :dΓ_in)
  @test haskey(d, :dΓ_out)
  @test haskey(d, :dΓ_bot)

  # Per-damping-zone keys
  @test haskey(d, :dΓ_d_1)
  @test haskey(d, :dΓ_d_2)

  # All values are Gridap Measures
  for v in values(d)
    @test v isa Gridap.CellData.Measure
  end
end