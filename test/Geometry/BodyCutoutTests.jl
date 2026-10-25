using Test
using Gridap
import HydroElasticFEM.Geometry as G

# =========================================================================
# BodyCutout: surface-piercing box bodies cut out of the tank mesh
#
# 2D tank: [0,8] × [−2,0], 32 × 8 cells (Δx = Δz = 0.25), free surface z = 0
# Box:     x ∈ [3,5], draft 0.5           → wetted length 2 + 2·0.5 = 3
# 3D tank: [0,8] × [0,4] × [−2,0], 16 × 8 × 4 cells (Δ = 0.5)
# Box:     x ∈ [3,5], y ∈ [1,3], draft 0.5 → wetted area 4 + 2·0.5·4 = 8
# =========================================================================

_map2(x) = VectorValue(x[1], x[2] - 2.0)
_map3(x) = VectorValue(x[1], x[2], x[3] - 2.0)
_bc_tank2(; kw...) = G.TankDomain(; L = 8.0, H = 2.0, nx = 32, ny = 8, map = _map2, kw...)
_bc_tank3(; kw...) = G.TankDomain(; L = 8.0, W = 4.0, H = 2.0, nx = 16, ny = 8, nz = 4, map = _map3, kw...)
_box2(; kw...) = G.BodyCutout(; L = 2.0, d = 0.5, x₀ = [3.0, 0.0], kw...)

function _bc_domains(tank)
  model = G.build_model(tank)
  model, G.get_integration_domains(G.build_triangulations(tank, model))
end
_meas(dom, key) = sum(∫(1.0)dom[key])

@testset "BodyCutout descriptor" begin
  b = _box2()
  @test b.domain_symbol == :Γ_b && b.W === nothing
  @test G._cutout_wetted_area(b) ≈ 3.0
  @test G._cutout_wetted_area(G.BodyCutout(L = 2.0, W = 1.0, d = 0.5, x₀ = [0.0, 0.0, 0.0])) ≈ 2.0 + 3.0
  @test_throws AssertionError G.BodyCutout(L = -1.0, d = 0.5, x₀ = [3.0, 0.0])
  @test_throws AssertionError G.BodyCutout(L = 1.0, d = 0.0, x₀ = [3.0, 0.0])
  @test_throws AssertionError G.BodyCutout(L = 1.0, d = 0.5, x₀ = [3.0, 0.0, 0.0])            # 3D needs W
  @test_throws AssertionError G.BodyCutout(L = 1.0, W = 1.0, d = 0.5, x₀ = [3.0, 0.0])        # W only in 3D
  # validated against the tank
  @test_throws ErrorException _bc_tank3(body_cutouts = [_box2()])
  @test_throws ErrorException _bc_tank2(body_cutouts = [_box2(), _box2(x₀ = [6.0, 0.0])])   # same symbol
  @test_throws ErrorException _bc_tank2(body_cutouts = [_box2()],
    structure_domains = [G.StructureDomain(L = 1.0, x₀ = [6.0, 0.0], domain_symbol = :Γ_b)])
end

@testset "BodyCutout 2D: cut mesh, tags and measures" begin
  @test G.build_model(_bc_tank2()) isa CartesianDiscreteModel   # no cut-outs: unchanged
  tank = _bc_tank2(body_cutouts = [_box2()])
  @test :body_cutouts in propertynames(tank)
  model, dom = _bc_domains(tank)
  @test model isa Gridap.Geometry.UnstructuredDiscreteModel
  @test num_cells(model) == 32 * 8 - 8 * 2
  @test _meas(dom, :dΩ) ≈ 16.0 - 1.0
  @test _meas(dom, :dΓ_b) ≈ 3.0
  @test _meas(dom, :dΓfs) ≈ 6.0                 # free surface over the body removed
  @test _meas(dom, :dΓin) ≈ 2.0 && _meas(dom, :dΓout) ≈ 2.0 && _meas(dom, :dΓbot) ≈ 8.0
  # fluid normal points into the body: ∫n = (0, B) (walls cancel), and
  # ∫ n_z over the bottom only, ∫ x n_x = −2·d·… from the two walls
  n = dom[:nΓ_b]
  @test sum(∫(n ⋅ VectorValue(1.0, 0.0))dom[:dΓ_b]) ≈ 0.0 atol = 1e-12
  @test sum(∫(n ⋅ VectorValue(0.0, 1.0))dom[:dΓ_b]) ≈ 2.0
  xΓ = get_physical_coordinate(get_triangulation(dom[:dΓ_b].quad))
  @test sum(∫((x -> x[1] - 4.0) ∘ xΓ * (n ⋅ VectorValue(1.0, 0.0)))dom[:dΓ_b]) ≈ -2 * 1.0 * 0.5
  # FE spaces build on the cut model
  V = TestFESpace(model, ReferenceFE(lagrangian, Float64, 2))
  @test num_free_dofs(V) > 0
  # the cut model keeps its tags when labelled again
  G.build_triangulations(tank, model)
  @test "Γ_b" in get_face_labeling(model).tag_to_name
end

@testset "BodyCutout 2D: several bodies and surface zones" begin
  tank = _bc_tank2(body_cutouts = [_box2(), _box2(L = 1.0, d = 1.0, x₀ = [6.0, 0.0], domain_symbol = :Γ_c)],
                   structure_domains = [G.StructureDomain(L = 1.0, x₀ = [1.0, 0.0], domain_symbol = :Γ_s)])
  model = G.build_model(tank)
  tr = G.build_triangulations(tank, model)
  dom = G.get_integration_domains(tr)
  @test _meas(dom, :dΓ_b) ≈ 3.0 && _meas(dom, :dΓ_c) ≈ 3.0
  @test _meas(dom, :dΓfs) ≈ 8.0 - 2.0 - 1.0 - 1.0     # bodies and structure removed
  @test sum(∫(1.0)Measure(tr[:Γ_s], 2)) ≈ 1.0
  @test tr[:body_symbols] == [:Γ_b, :Γ_c]
end

@testset "BodyCutout 2D: invalid boxes" begin
  # box face off the mesh lines
  @test_throws ErrorException G.build_model(_bc_tank2(body_cutouts = [_box2(L = 2.1)]))
  @test_throws ErrorException G.build_model(_bc_tank2(body_cutouts = [_box2(d = 0.6)]))
  # submerged box (top below the free surface)
  @test_throws ErrorException G.build_model(_bc_tank2(body_cutouts = [_box2(x₀ = [3.0, -0.5])]))
  # box thinner than a cell removes nothing
  @test_throws ErrorException G.build_model(_bc_tank2(body_cutouts = [_box2(L = 0.1)]))
  # triangulations need the cut model
  tank = _bc_tank2(body_cutouts = [_box2()])
  plain = G.build_model(_bc_tank2())
  @test_throws ErrorException G.build_triangulations(tank, plain)
end

@testset "BodyCutout 3D" begin
  b = G.BodyCutout(L = 2.0, W = 2.0, d = 0.5, x₀ = [3.0, 1.0, 0.0])
  model, dom = _bc_domains(_bc_tank3(body_cutouts = [b]))
  @test _meas(dom, :dΩ) ≈ 64.0 - 2.0
  @test _meas(dom, :dΓ_b) ≈ 8.0
  @test _meas(dom, :dΓfs) ≈ 32.0 - 4.0
  @test _meas(dom, :dΓlateral) ≈ 2 * 16.0
  @test sum(∫(dom[:nΓ_b] ⋅ VectorValue(0.0, 0.0, 1.0))dom[:dΓ_b]) ≈ 4.0
  @test_throws ErrorException G.build_model(_bc_tank3(body_cutouts = [G.BodyCutout(L = 2.0, W = 1.7, d = 0.5, x₀ = [3.0, 1.0, 0.0])]))
end
