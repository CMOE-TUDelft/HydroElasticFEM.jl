using Test
using Gridap
using LinearAlgebra
using Gridap.FESpaces: ConstantFESpace
using Gridap.ODEs: TransientTrialFESpace
using Gridap.Geometry: DiscreteModelPortion, get_face_labeling, add_tag_from_tags!
import HydroElasticFEM.Geometry as G

# ==========================================================================
# Spike (rigid bodies, PR 0): global rigid-body DOFs coupled to a fluid
#
# Questions:
#  S1  Can a vector-valued ConstantFESpace (one VectorValue{n} unknown per
#      body) defined on Ω be integrated against the generalized normal
#      N(x) = Operation(f)(x, n) on a boundary triangulation, in a complex
#      multi-field system with the potential ϕ?
#  S2  Does a constant matrix term ∫ ζ⋅(A⋅ξ)/|Γ| dΓ assemble to exactly A?
#  S3  Does the same space work as a TransientTrialFESpace?
#  S4  Can a finite-draft body be represented without Gmsh, as a box of
#      cells removed from a Cartesian tank, with a correct wetted surface
#      and fluid normal?
#
# Outcome: GO on all four.  Notes for the implementation:
#  - The body DOF space can live on Ω and be integrated on any boundary
#    triangulation; N = Operation(f)(get_physical_coordinate(Γ), n_Γ) works.
#  - Constant matrices assemble exactly via ∫ ζ⋅(A⋅ξ)/|Γ_b| dΓ_b.
#  - Gridap's boundary normal points out of the fluid (into the body).
#  - Finite-draft boxes work without Gmsh through DiscreteModelPortion;
#    build the spaces from the portion *model* (Triangulation(portion)
#    fails Gridap's active-model check in TestFESpace).
# ==========================================================================

# 2D generalized normal for modes (surge, heave, pitch) about x_r, with the
# rotation about y: u = (ξ₁ + α r_z, ξ₃ − α r_x), r = x − x_r.
_spike_gen_normal(x_r) = (x, n) -> VectorValue(n[1], n[2], (x[2] - x_r[2]) * n[1] - (x[1] - x_r[1]) * n[2])

@testset "S1–S3: vector ConstantFESpace on a zero-draft body" begin
  tank = G.TankDomain(L = 8.0, H = 1.0, nx = 16, ny = 2,
    structure_domains = [G.StructureDomain(L = 2.0, x₀ = [3.0, 1.0], domain_symbol = :Γ_b)])
  tr = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(tr; degree = Dict(:dΩ => 4, :Γ_b => 4))
  Ω, Γb, dΓb = tr[:Ω], tr[:Γ_b], dom[:dΓ_b]
  x_r = VectorValue(4.0, 1.0)

  Vϕ = TestFESpace(Ω, ReferenceFE(lagrangian, Float64, 1); vector_type = Vector{ComplexF64})
  Vξ = ConstantFESpace(Ω; field_type = VectorValue{3,Float64}, vector_type = Vector{ComplexF64})
  @test num_free_dofs(Vξ) == 3
  Y = MultiFieldFESpace([Vϕ, Vξ])
  X = MultiFieldFESpace([TrialFESpace(Vϕ), TrialFESpace(Vξ)])

  nΓ = get_normal_vector(Γb)
  xΓ = get_physical_coordinate(Γb)
  N = Operation(_spike_gen_normal(x_r))(xΓ, nΓ)

  # S1: fluid row of the kinematic coupling  −∫ w (ξ·N)
  a1((ϕ, ξ), (w, ζ)) = ∫(-w * (ξ ⋅ N))dΓb
  A1 = assemble_matrix(a1, X, Y)
  nϕ = num_free_dofs(Vϕ)
  # reference: −∫ w N_j for each mode j, assembled on the ϕ space alone
  for j in 1:3
    ej = VectorValue(ntuple(k -> float(k == j), 3))
    bj = assemble_vector(w -> ∫(-w * (N ⋅ ej))dΓb, Vϕ)
    @test A1[1:nϕ, nϕ + j] ≈ bj
  end
  # top boundary: n = +e_z (out of the fluid), so surge column is zero,
  # heave column integrates to −L, pitch column to +∫(x − x_r) = 0
  @test norm(A1[1:nϕ, nϕ + 1]) < 1e-12
  @test sum(A1[1:nϕ, nϕ + 2]) ≈ -2.0
  @test abs(sum(A1[1:nϕ, nϕ + 3])) < 1e-12

  # S1: body row (pressure force)  ∫ (ζ·N) ϕ  is the transpose structure
  a2((ϕ, ξ), (w, ζ)) = ∫((ζ ⋅ N) * ϕ)dΓb
  A2 = assemble_matrix(a2, X, Y)
  @test A2[nϕ+1:end, 1:nϕ] ≈ -transpose(A1[1:nϕ, nϕ+1:end])

  # S1: zero-draft hydrostatics ∫ g (ζ·N)(ξ·N) gives C33 = g L, C35 = 0,
  # C55 = g L³/12 (per unit width, normalized by ρ)
  g = 9.81
  a3((ϕ, ξ), (w, ζ)) = ∫(g * (ζ ⋅ N) * (ξ ⋅ N))dΓb
  C = Matrix(assemble_matrix(a3, X, Y)[nϕ+1:end, nϕ+1:end])
  @test C[2, 2] ≈ g * 2.0
  @test abs(C[2, 3]) < 1e-10
  @test C[3, 3] ≈ g * 2.0^3 / 12
  @test norm(C[1, :]) < 1e-12

  # S2: constant matrix term assembles to exactly A
  Am = [2.0 0.1 0.3; 0.1 5.0 -0.2; 0.3 -0.2 7.0]
  At = TensorValue{3,3}(Tuple(Am))   # column-major, same as Am
  area = sum(∫(1.0)dΓb)
  a4((ϕ, ξ), (w, ζ)) = ∫(ζ ⋅ (At ⋅ ξ) / area)dΓb
  @test Matrix(assemble_matrix(a4, X, Y)[nϕ+1:end, nϕ+1:end]) ≈ Am

  # Two bodies in one system
  Vξ2 = ConstantFESpace(Ω; field_type = VectorValue{3,Float64}, vector_type = Vector{ComplexF64})
  Y2 = MultiFieldFESpace([Vϕ, Vξ, Vξ2])
  X2 = MultiFieldFESpace([TrialFESpace(Vϕ), TrialFESpace(Vξ), TrialFESpace(Vξ2)])
  a5((ϕ, ξ, χ), (w, ζ, ψ)) = ∫(ζ ⋅ (At ⋅ ξ) / area + ψ ⋅ χ / area - ζ ⋅ χ / area)dΓb
  A5 = Matrix(assemble_matrix(a5, X2, Y2))
  @test A5[nϕ+1:nϕ+3, nϕ+1:nϕ+3] ≈ Am
  @test A5[nϕ+1:nϕ+3, nϕ+4:nϕ+6] ≈ -I(3)

  # S3: transient trial space (real, time domain)
  Vϕr = TestFESpace(Ω, ReferenceFE(lagrangian, Float64, 1))
  Vξr = ConstantFESpace(Ω; field_type = VectorValue{3,Float64})
  Yt = MultiFieldFESpace([Vϕr, Vξr])
  Xt = MultiFieldFESpace([TransientTrialFESpace(Vϕr), TransientTrialFESpace(Vξr)])
  At_ = assemble_matrix(a1, Xt(0.0), Yt)
  @test size(At_) == (num_free_dofs(Vϕr) + 3, num_free_dofs(Vϕr) + 3)
end

@testset "S4: finite-draft box cut out of a Cartesian tank" begin
  # Tank (0,8)×(0,2), cells 0.25 × 0.25; body box x ∈ (3,5), z ∈ (1.5,2)
  # (beam B = 2, draft d = 0.5, free surface at z = 2).
  model = CartesianDiscreteModel((0.0, 8.0, 0.0, 2.0), (32, 8))
  Ωall = Triangulation(model)
  cc = vec([sum(c) / length(c) for c in get_cell_coordinates(Ωall)])
  in_body(x) = 3.0 < x[1] < 5.0 && x[2] > 1.5
  fluid_cells = findall(!in_body, cc)
  fmodel = DiscreteModelPortion(model, fluid_cells)
  Ω = Triangulation(fmodel)
  Γall = Boundary(fmodel)
  fc = [sum(c) / length(c) for c in get_cell_coordinates(Γall)]
  on_body(x) = 2.999 < x[1] < 5.001 && 1.499 < x[2] < 2.001 &&
               !(isapprox(x[2], 2.0) )
  body_faces = findall(on_body, fc)
  Γb = view(Γall, body_faces)
  dΓb = Measure(Γb, 2)
  nΓ = get_normal_vector(Γb)

  # wetted length: two sides (0.5 each) + bottom (2) = 3
  @test sum(∫(1.0)dΓb) ≈ 3.0
  # fluid outward normal points into the body: on the bottom n = +e_z
  @test sum(∫(nΓ ⋅ VectorValue(0.0, 1.0))dΓb) ≈ 2.0
  @test abs(sum(∫(nΓ ⋅ VectorValue(1.0, 0.0))dΓb)) < 1e-12
  # waterplane and displaced volume from the wetted surface (body normal = −nΓ):
  #   A_wp = −∫ n_b,z = ∫ n_z,  V = −∫ (z − z_fs) n_b,z = ∫ (z − z_fs) n_z
  z_fs = 2.0
  xΓ = get_physical_coordinate(Γb)
  @test sum(∫(nΓ ⋅ VectorValue(0.0, 1.0))dΓb) ≈ 2.0                     # A_wp = B
  V = -sum(∫((x -> x[2] - z_fs) ∘ xΓ * (nΓ ⋅ VectorValue(0.0, 1.0)))dΓb)
  @test V ≈ 2.0 * 0.5                                                      # B d

  # Spaces on the cut fluid model (built from the model: Triangulation(fmodel)
  # trips an active-model check in Gridap), coupled on the wetted surface
  Vϕ = TestFESpace(fmodel, ReferenceFE(lagrangian, Float64, 1); vector_type = Vector{ComplexF64})
  Vξ = ConstantFESpace(fmodel; field_type = VectorValue{3,Float64}, vector_type = Vector{ComplexF64})
  Y = MultiFieldFESpace([Vϕ, Vξ]); X = MultiFieldFESpace([TrialFESpace(Vϕ), TrialFESpace(Vξ)])
  N = Operation(_spike_gen_normal(VectorValue(4.0, 2.0)))(xΓ, nΓ)
  A = assemble_matrix(((ϕ, ξ), (w, ζ)) -> ∫(-w * (ξ ⋅ N))dΓb, X, Y)
  nϕ = num_free_dofs(Vϕ)
  # Σ_i of the column = −∫ N_j: surge 0 (sides cancel), heave −B,
  # pitch −∫(r_z n_x − r_x n_z) = −(0 + 0) by symmetry
  @test abs(sum(A[1:nϕ, nϕ + 1])) < 1e-12
  @test sum(A[1:nϕ, nϕ + 2]) ≈ -2.0
  @test abs(sum(A[1:nϕ, nϕ + 3])) < 1e-12
end
