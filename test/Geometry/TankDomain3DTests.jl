using Test
using Gridap
import HydroElasticFEM.Geometry as G

# =========================================================================
# TankDomain{3}: floating structures, damping zones, structure boundaries
#
# Reference tank: [0,8] × [0,4] × [0,1], 8 × 8 × 2 cells  (Δx = 1, Δy = 0.5)
# Plate:          x ∈ [2,4], y ∈ [1,3] at z = 1  → 2 × 4 = 8 surface cells
# Inlet damping:  x ∈ [0,1], full width          → 1 × 8 = 8 surface cells
# =========================================================================

const _NX, _NY, _NZ = 8, 8, 2
_plate(; kw...) = G.StructureDomain(; L = 2.0, W = 2.0, x₀ = [2.0, 1.0, 1.0], kw...)
_inlet() = G.DampingZone(L = 1.0, x₀ = [0.0, 0.0, 1.0], domain_symbol = :Γ_d_in)
_tank3(; kw...) = G.TankDomain(; L = 8.0, W = 4.0, H = 1.0, nx = _NX, ny = _NY, nz = _NZ, kw...)

@testset "Surface-zone descriptors infer dimensions from x₀" begin
  s2 = G.StructureDomain(L = 1.0, x₀ = [1.0, 1.0])
  @test (s2.ambient_dim, s2.manifold_dim) == (2, 1)
  @test s2.W === nothing
  @test s2.boundary_tag == "Γ_s_boundary"

  s3 = _plate(domain_symbol = :Γ_p)
  @test (s3.ambient_dim, s3.manifold_dim) == (3, 2)
  @test s3.boundary_tag == "Γ_p_boundary"
  @test s3 isa G.AbstractSurfaceZone
  @test _inlet() isa G.AbstractSurfaceZone
end

@testset "surface_mask — 3D plate and full-width damping zone" begin
  pm = G.surface_mask(_plate())
  @test pm([VectorValue(2.5, 1.5, 1.0)])
  @test !pm([VectorValue(2.5, 3.5, 1.0)])   # outside in y
  @test !pm([VectorValue(4.5, 1.5, 1.0)])   # outside in x
  @test !pm([VectorValue(2.5, 1.5, 0.5)])   # not on the surface

  dm = G.surface_mask(_inlet())
  @test dm([VectorValue(0.5, 3.9, 1.0)])
  @test dm([VectorValue(0.5, -100.0, 1.0)]) # W = nothing: no y constraint
  @test !dm([VectorValue(1.5, 0.5, 1.0)])
end

@testset "TankDomain{3} input validation" begin
  # 3D structure without W
  @test_throws ErrorException _tank3(structure_domains = [G.StructureDomain(L = 2.0, x₀ = [2.0, 1.0, 1.0])])
  # dimension mismatch
  @test_throws ErrorException _tank3(structure_domains = [G.StructureDomain(L = 2.0, x₀ = [2.0, 1.0])])
  @test_throws ErrorException G.TankDomain(L = 4.0, H = 1.0, nx = 4, ny = 2,
    damping_zones = [_inlet()])
  # non-positive extent
  @test_throws ErrorException _tank3(structure_domains = [_plate(W = 0.0)])
  # joints are 2D only
  @test_throws ErrorException _tank3(structure_domains = [_plate()],
    joint_domains = [G.JointDomain(location = [3.0, 2.0, 1.0], domain_symbol = :dΛj, normal_symbol = :nΛj)])
  # structure that hits no cell is reported at build time
  tank = _tank3(structure_domains = [_plate(x₀ = [2.0, 1.0, 0.5])])
  @test_throws ErrorException G.build_triangulations(tank, G.build_model(tank))
end

@testset "build_triangulations — 3D plate + damping partition" begin
  tank = _tank3(structure_domains = [_plate()], damping_zones = [_inlet()])
  model = G.build_model(tank)
  tr = G.build_triangulations(tank, model)

  n_surface = _NX * _NY
  @test num_cells(tr[:Γ]) == n_surface
  @test num_cells(tr[:Γη]) == 8
  @test num_cells(tr[:Γ_s]) == 8
  @test num_cells(tr[:Γ_d_in]) == 8
  @test num_cells(tr[:Γfs]) == n_surface - 16
  @test num_cells(tr[:Γκ]) == n_surface - 8
  @test num_cells(tr[:Ω]) == _NX * _NY * _NZ
  @test num_cells(tr[:Γbot]) == n_surface
  @test num_cells(tr[:Γin]) == _NY * _NZ
  @test num_cells(tr[:Γout]) == _NY * _NZ
  @test num_cells(tr[:Γlateral]) == 2 * _NX * _NZ

  # C/DG skeleton of a 2×4 cell patch: 1·4 + 2·3 interior edges
  @test num_cells(tr[:Λη]) == 10
  @test num_cells(tr[:Λ_structures][1]) == 10
  # structure boundary: 2·(2 + 4) edges
  @test num_cells(tr[:∂Γ_structures][1]) == 12
  @test num_cells(tr[:∂Γη]) == 12

  @test num_cells(G.get_boundary(tank, "structure")) == 8
  @test num_cells(G.get_boundary(tank, "surface")) == n_surface

  # Rebuilding on the same model is a no-op for labels and edge tags
  tr2 = G.build_triangulations(tank, model)
  @test num_cells(tr2[:Γη]) == 8
end

@testset "build_triangulations — free surface at z = 0, centred in y" begin
  # Same layout as the Yago / tmp_3D scripts: y ∈ [-W/2, W/2], z ∈ [-H, 0]
  cart = G.CartesianDomain(mins = (0.0, -2.0, -1.0), maxs = (8.0, 2.0, 0.0), parts = (_NX, _NY, _NZ))
  tank = G.TankDomain(cart;
    structure_domains = [G.StructureDomain(L = 2.0, W = 2.0, x₀ = [2.0, -1.0, 0.0])],
    damping_zones = [G.DampingZone(L = 1.0, x₀ = [0.0, -2.0, 0.0])],
  )
  tr = G.build_triangulations(tank, G.build_model(tank))
  @test num_cells(tr[:Γη]) == 8
  @test num_cells(tr[:Γ_d]) == 8
  @test num_cells(tr[:Γfs]) == _NX * _NY - 16
  zc = [sum(n[3] for n in c) / length(c) for c in get_cell_coordinates(tr[:Γη])]
  @test all(abs.(zc) .< 1e-12)
end

@testset "Structure boundary tags — Dirichlet on plate edges" begin
  tank = _tank3(structure_domains = [_plate()])
  model = G.build_model(tank)
  tr = G.build_triangulations(tank, model)
  labels = get_face_labeling(model)
  @test "Γ_s_boundary" in labels.tag_to_name
  @test "structure_boundary" in labels.tag_to_name

  # Relabelling must not change existing tags
  @test num_cells(Boundary(model, tags = "surface")) == _NX * _NY
  @test num_cells(Boundary(model, tags = "boundary")) == num_cells(Boundary(model))

  # Q2 on a 2×4 plate: 5·9 = 45 nodes, 2·(4 + 8) = 24 on the edge
  reffe = ReferenceFE(lagrangian, Float64, 2)
  V0 = TestFESpace(tr[:Γη], reffe; conformity = :H1)
  V = TestFESpace(tr[:Γη], reffe; conformity = :H1, dirichlet_tags = "Γ_s_boundary")
  @test num_free_dofs(V0) == 45
  @test num_dirichlet_dofs(V) == 24
  @test num_free_dofs(V) == 21
end

@testset "Structure boundary tags — adjacent plates share an edge" begin
  s1 = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [2.0, 1.0, 1.0], domain_symbol = :Γ_s1)
  s2 = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [4.0, 1.0, 1.0], domain_symbol = :Γ_s2)
  tank = _tank3(structure_domains = [s1, s2])
  model = G.build_model(tank)
  tr = G.build_triangulations(tank, model)

  reffe = ReferenceFE(lagrangian, Float64, 1)
  # Q1 on each 2×4 patch: 3·5 = 15 nodes, 12 on its own boundary
  V1 = TestFESpace(tr[:Γ_s1], reffe; conformity = :H1, dirichlet_tags = "Γ_s1_boundary")
  V2 = TestFESpace(tr[:Γ_s2], reffe; conformity = :H1, dirichlet_tags = "Γ_s2_boundary")
  @test num_dirichlet_dofs(V1) == 12
  @test num_dirichlet_dofs(V2) == 12
  # Union boundary tag on Γη (4×4 patch, 5·5 nodes): both outlines, shared edge x=4 counted once
  Vη = TestFESpace(tr[:Γη], reffe; conformity = :H1, dirichlet_tags = "structure_boundary")
  @test num_dirichlet_dofs(Vη) == 12 + 12 - 5
  # Geometric boundary of the union excludes the shared edge
  @test num_cells(tr[:∂Γη]) == 2 * (4 + 4)
  # Union skeleton includes the shared edge, per-structure skeletons do not
  @test num_cells(tr[:Λη]) == num_cells(tr[:Λ_structures][1]) + num_cells(tr[:Λ_structures][2]) + 4
end

@testset "Structure boundary — 2D beam end points" begin
  tank = G.TankDomain(L = 4.0, H = 1.0, nx = 8, ny = 2,
    structure_domains = [G.StructureDomain(L = 1.0, x₀ = [1.0, 1.0])])
  model = G.build_model(tank)
  tr = G.build_triangulations(tank, model)
  @test num_cells(tr[:∂Γ_structures][1]) == 2
  @test num_cells(Boundary(model, tags = "surface")) == 8
  V = TestFESpace(tr[:Γη], ReferenceFE(lagrangian, Float64, 2); conformity = :H1,
    dirichlet_tags = "Γ_s_boundary")
  @test num_dirichlet_dofs(V) == 2
end

@testset "get_integration_domains — 3D structure keys" begin
  r1 = G.ResonatorDomain(location = [3.0, 2.0, 1.0])
  tank = _tank3(structure_domains = [_plate()], damping_zones = [_inlet()], resonator_domains = [r1])
  dom = G.get_integration_domains(G.build_triangulations(tank, G.build_model(tank)); degree = 4)

  for key in (:dΓη, :dΓη_1, :nΓη_1, :dΛη, :n_Λ_η, :h_η, :dΛη_1, :n_Λ_η_1, :h_η_1,
              :dΛ∂η, :n_Λ∂η, :dΛ∂η_1, :n_Λ∂η_1, :dΓd_1, :nΓd_1, :dΓκ, :dΓfs, :δ_p)
    @test haskey(dom, key)
  end

  # h is a length, not an area: sqrt(Δx·Δy)
  @test dom[:h_η] ≈ sqrt(1.0 * 0.5)
  @test dom[:h_η_1] ≈ dom[:h_η]

  @test sum(∫(1.0)dom[:dΓη]) ≈ 4.0
  @test sum(∫(1.0)dom[:dΛ∂η]) ≈ 8.0
  @test sum(∫(1.0)dom[:dΓd_1]) ≈ 4.0
  # Outward conormal: ∮ x·n_x = area (divergence theorem on the plate)
  n = dom[:n_Λ∂η_1]
  xcoord(x) = x[1]
  @test sum(∫(xcoord * (n ⋅ VectorValue(1.0, 0.0, 0.0)))dom[:dΛ∂η_1]) ≈ 4.0 atol = 1e-10
end

@testset "get_integration_domains — 2D h_η unchanged" begin
  tank = G.TankDomain(L = 4.0, H = 1.0, nx = 8, ny = 2,
    structure_domains = [G.StructureDomain(L = 1.0, x₀ = [1.0, 1.0])])
  dom = G.get_integration_domains(G.build_triangulations(tank, G.build_model(tank)))
  @test dom[:h_η] ≈ 0.5
  @test haskey(dom, :dΛ∂η_1)
  @test sum(∫(1.0)dom[:dΛ∂η_1]) ≈ 2.0   # two end points, unit weight each
end

# =========================================================================
# StructureConnection: interface between two plates with separate fields
# =========================================================================

_pa(; x = 2.0, y = 1.0) = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [x, y, 1.0], domain_symbol = :Γ_a)
_pb(; x = 4.0, y = 1.0) = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [x, y, 1.0], domain_symbol = :Γ_b)
_conn(; a = :Γ_a, b = :Γ_b) = G.StructureConnection(a = a, b = b, domain_symbol = :dΛ_ab, normal_symbol = :n_Λ_ab)

@testset "StructureConnection — interface measure, orientation and traces" begin
  ex, ey = VectorValue(1.0, 0.0, 0.0), VectorValue(0.0, 1.0, 0.0)
  for (pa, pb, conn, len, nref) in (
      (_pa(), _pb(), _conn(), 2.0, ex),                         # a left of b, edge x = 4
      (_pa(), _pb(), _conn(a = :Γ_b, b = :Γ_a), 2.0, -ex),      # reversed roles
      (_pa(), _pb(x = 2.0, y = 3.0), _conn(), 2.0, ey))         # b above a in y, edge y = 3
    tank = _tank3(structure_domains = [pa, pb], structure_connections = [conn])
    @test tank.structure_connections == [conn]
    tr = G.build_triangulations(tank, G.build_model(tank))
    @test haskey(tr, :dΛ_ab)
    dom = G.get_integration_domains(tr; degree = 4)
    dΛ, n = dom[:dΛ_ab], dom[:n_Λ_ab]
    @test sum(∫(1.0)dΛ) ≈ len
    @test sum(∫(n.⁺ ⋅ nref)dΛ) ≈ len          # n⁺ points out of structure a
    @test sum(∫(n.⁻ ⋅ nref)dΛ) ≈ -len

    # separate FE fields on the two plates: plus trace from a, minus from b
    lag = ReferenceFE(lagrangian, Float64, 2)
    fa(x) = 1.0 + x[1] * x[2]
    fb(x) = 5.0 - x[1]
    Ta, Tb = tr[conn.a], tr[conn.b]
    ηa = interpolate(fa, TestFESpace(Ta, lag))
    ηb = interpolate(fb, TestFESpace(Tb, lag))
    Λ = tr[:dΛ_ab]
    @test sum(∫(ηa.⁺)dΛ) ≈ sum(∫(CellField(fa, Λ))dΛ)
    @test sum(∫(ηb.⁻)dΛ) ≈ sum(∫(CellField(fb, Λ))dΛ)
  end
end

@testset "StructureConnection — validation" begin
  @test_throws ErrorException _tank3(structure_domains = [_pa(), _pb()],
                                     structure_connections = [_conn(b = :Γ_zz)])
  @test_throws ErrorException _tank3(structure_domains = [_pa(), _pb()],
                                     structure_connections = [_conn(b = :Γ_a)])
  # plates that do not touch
  t = _tank3(structure_domains = [_pa(), _pb(x = 5.0)], structure_connections = [_conn()])
  @test_throws ErrorException G.build_triangulations(t, G.build_model(t))
end
    
# JointLineDomain: line joints on a 3D plate
#
# Plate x ∈ [2,6], y ∈ [1,3] → 4 × 4 surface cells (Δx = 1, Δy = 0.5),
# 24 interior skeleton facets.  A 2 × 2 floater grid of 2 m × 1 m floaters
# has the connection lines x = 4 (length 2, 4 facets) and y = 2 (length 4,
# 4 facets).
# =========================================================================

_big_plate() = G.StructureDomain(L = 4.0, W = 2.0, x₀ = [2.0, 1.0, 1.0])
_hinges(; kw...) = G.hinge_grid(; x₀ = [2.0, 1.0, 1.0], a = 2.0, b = 1.0, nfx = 2, nfy = 2,
                                domain_symbol = :dΛh, normal_symbol = :n_Λh, kw...)

@testset "hinge_grid — interior connection lines" begin
  j = _hinges()
  @test j isa G.JointLineDomain
  @test length(j.segments) == 2
  @test j.segments[1] == ([4.0, 1.0, 1.0], [4.0, 3.0, 1.0])
  @test j.segments[2] == ([2.0, 2.0, 1.0], [6.0, 2.0, 1.0])
  @test isempty(G.hinge_grid(x₀ = [0.0, 0.0, 0.0], a = 1.0, nfx = 1, nfy = 1,
                             domain_symbol = :d, normal_symbol = :n).segments)
end

@testset "JointLineDomain — skeleton partition and measure" begin
  tank = _tank3(structure_domains = [_big_plate()], joint_domains = [_hinges()])
  trians = G.build_triangulations(tank, G.build_model(tank))
  @test num_cells(trians[:Λ_joints][1]) == 8
  @test num_cells(trians[:Λη]) == 24 - 8
  @test num_cells(trians[:Λ_structures][1]) == 24 - 8
  @test haskey(trians, :dΛh)

  dom = G.get_integration_domains(trians; degree = 4)
  @test haskey(dom, :dΛh) && haskey(dom, :n_Λh)
  @test sum(∫(1.0)dom[:dΛh]) ≈ 6.0
  # Joint facets carry the in-plane normal to the line (no vertical part)
  ez = VectorValue(0.0, 0.0, 1.0)
  @test sum(∫(abs(dom[:n_Λh].⁺ ⋅ ez))dom[:dΛh]) ≈ 0.0 atol = 1e-12
  # Total skeleton length is unchanged: joint + remaining C/DG facets
  @test sum(∫(1.0)dom[:dΛh]) + sum(∫(1.0)dom[:dΛη]) ≈ sum(∫(1.0)Measure(Skeleton(trians[:Γη]), 2))
end

@testset "JointLineDomain — validation" begin
  mk(segs) = G.JointLineDomain(segments = segs, domain_symbol = :dΛh, normal_symbol = :n_Λh)
  build(j) = (t = _tank3(structure_domains = [_big_plate()], joint_domains = [j]);
              G.build_triangulations(t, G.build_model(t)))

  # line between mesh faces: no facet matches
  @test_throws ErrorException build(mk([([4.3, 1.0, 1.0], [4.3, 3.0, 1.0])]))
  # line running past the plate: selected length ≠ segment length
  @test_throws ErrorException build(mk([([4.0, 0.0, 1.0], [4.0, 4.0, 1.0])]))
  # overlapping segments: covered length counted once
  @test_throws ErrorException build(mk([([4.0, 1.0, 1.0], [4.0, 3.0, 1.0]),
                                        ([4.0, 1.0, 1.0], [4.0, 2.0, 1.0])]))
  # JointLineDomain is 3D only, JointDomain is 2D only
  @test_throws ErrorException G.TankDomain(L = 4.0, H = 1.0, nx = 8, ny = 2,
    structure_domains = [G.StructureDomain(L = 1.0, x₀ = [1.0, 1.0])],
    joint_domains = [mk([([1.5, 1.0], [1.5, 1.0])])])
end
