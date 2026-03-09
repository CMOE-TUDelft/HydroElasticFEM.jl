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