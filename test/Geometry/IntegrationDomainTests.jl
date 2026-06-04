using Test
using Gridap
import HydroElasticFEM.Geometry as G

@testset "get_integration_domains" begin
  s1 = G.StructureDomain(L=1.0, x₀=[1.5, 1.0])
  d1 = G.DampingZone(L=0.5, x₀=[0.0, 1.0], domain_symbol=:Γ_d_in)
  d2 = G.DampingZone(L=0.5, x₀=[3.5, 1.0], domain_symbol=:Γ_d_out)
  tank = G.TankDomain(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains=[s1], damping_zones=[d1, d2])
  model = G.build_model(tank)
  trians = G.build_triangulations(tank, model)

  d = G.get_integration_domains(trians; degree=4)

  @test d isa G.IntegrationDomains

  # Core keys present
  @test haskey(d, :dΩ)
  @test haskey(d, :dΓfs)
  @test haskey(d, :dΓκ)
  @test haskey(d, :dΓη)
  @test haskey(d, :dΓin)
  @test haskey(d, :dΓout)
  @test haskey(d, :dΓbot)

  # Per-damping-zone and structure-skeleton keys
  @test haskey(d, :dΓd_1)
  @test haskey(d, :dΓd_2)
  @test haskey(d, :nΓd_1)
  @test haskey(d, :nΓd_2)
  @test haskey(d, :dΛη)
  @test haskey(d, :n_Λ_η)
  @test haskey(d, :h_η)

  # Measures, normals, and scalar mesh sizes are all stored in the integration-domain dictionary.
  for (k, v) in pairs(d.data)
    if startswith(String(k), "nΓ")
      @test v isa Gridap.CellData.GenericCellField
    elseif startswith(String(k), "n_Λ")
      @test (v isa Gridap.CellData.GenericCellField) || (v isa Gridap.Geometry.SkeletonPair)
    elseif startswith(String(k), "h_")
      @test v isa Real
    else
      @test v isa Gridap.CellData.Measure
    end
  end
end

@testset "get_integration_domains — automatic joint keys" begin
  s1 = G.StructureDomain(L=1.0, x₀=[1.5, 1.0])
  j1 = G.JointDomain(location=[2.0, 1.0], domain_symbol=:dΛj_1, normal_symbol=:n_Λ_j_1)
  tank = G.TankDomain(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains=[s1],
    joint_domains=[j1])

  model = G.build_model(tank)
  trians = G.build_triangulations(tank, model)
  d = G.get_integration_domains(trians; degree=4)

  @test haskey(d, :dΛj_1)
  @test haskey(d, :n_Λ_j_1)
  @test d[:dΛj_1] isa Gridap.CellData.Measure
  @test !isnothing(d[:n_Λ_j_1])
end

@testset "get_integration_domains — user-defined space domains" begin
  structure = G.StructureDomain(L=1.0, x₀=[1.5, 1.0], domain_symbol=:Γ_s_custom)
  damping = G.DampingZone(L=0.5, x₀=[0.0, 1.0], domain_symbol=:Γ_d_custom)
  tank = G.TankDomain(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains=[structure], damping_zones=[damping])

  model = G.build_model(tank)
  trians = G.build_triangulations(tank, model)

  d = G.get_integration_domains(trians; degree=Dict(:Γ_s_custom => 6, :Γ_d_custom => 8))

  @test haskey(trians, :Γ_s_custom)
  @test haskey(trians, :Γ_d_custom)
  @test haskey(d, :dΓ_s_custom)
  @test haskey(d, :dΓ_d_custom)
  @test d[:dΓ_s_custom] isa Gridap.CellData.Measure
  @test d[:dΓ_d_custom] isa Gridap.CellData.Measure
end

@testset "get_integration_domains — user-defined measure-key degrees" begin
  structure = G.StructureDomain(L=1.0, x₀=[1.5, 1.0], domain_symbol=:Γ_s_alt)
  tank = G.TankDomain(L=4.0, H=1.0, nx=40, ny=4, structure_domains=[structure])

  model = G.build_model(tank)
  trians = G.build_triangulations(tank, model)

  d = G.get_integration_domains(trians; degree=Dict(:dΓ_s_alt => 6))

  @test haskey(trians, :Γ_s_alt)
  @test haskey(d, :dΓ_s_alt)
  @test d[:dΓ_s_alt] isa Gridap.CellData.Measure
end

@testset "get_integration_domains — metadata degree keys are ignored" begin
  tank = G.TankDomain(L=4.0, H=1.0, nx=20, ny=2)
  model = G.build_model(tank)
  trians = G.build_triangulations(tank, model)

  d = G.get_integration_domains(trians; degree=Dict(:Γ_structures => 6, :joint_domains => 6))

  @test !haskey(d, :dΓ_structures)
  @test !haskey(d, :djoint_domains)
end
