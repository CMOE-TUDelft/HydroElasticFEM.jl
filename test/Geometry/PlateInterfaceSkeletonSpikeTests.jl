using Test
using Gridap
using LinearAlgebra
using Gridap.Geometry: CompositeTriangulation, get_active_model, get_grid_topology
using Gridap.ReferenceFEs: get_faces
import HydroElasticFEM.Geometry as G

# ==========================================================================
# Spike (Phase C1): interface skeleton between two separate plates
#
# Question: can two plates with their own FE spaces (η_a on Γ_a, η_b on Γ_b)
# be coupled across their common edge with a Gridap SkeletonTriangulation
# whose plus side always lies in Γ_a and minus side in Γ_b?
#
# Construction: on the active model M of Γη (the union of the plates), take
# the facets with one cell in each plate and build
#     SkeletonTriangulation(BoundaryTriangulation(M, faces, lcell_plus),
#                           BoundaryTriangulation(M, faces, lcell_minus))
# with per-facet local cell indices, wrapped as CompositeTriangulation(Γη, ·).
#
# Outcome: GO. Traces, gradients, normals and multi-field assembly are
# correct.  One caveat: on this hand-built skeleton, combining a minus-side
# trace with a constant (e.g. `∇(η_b).⁻ ⋅ e₁`) evaluates to zero, while
# `(∇(η_b) ⋅ e₁).⁻` and `∇(η_b).⁻ ⋅ n.⁻` are correct (Gridap's native
# Skeleton does not show this).  Interface forms must therefore combine
# fields with constants before taking the ⁺/⁻ trace, or use the normals.
# ==========================================================================

# Interface skeleton between the Γη cells for which `in_a` is true / false.
function _spike_interface_skeleton(Γη, in_a)
  M = get_active_model(Γη)
  Dc = num_cell_dims(M)
  f2c = get_faces(get_grid_topology(M), Dc - 1, Dc)
  cc = [sum(c) / length(c) for c in get_cell_coordinates(Γη)]
  nf = num_faces(M, Dc - 1)
  faces = Int[]
  lcp, lcm = ones(Int, nf), ones(Int, nf)
  for f in 1:nf
    cells = f2c[f]
    length(cells) == 2 || continue
    a1, a2 = in_a(cc[cells[1]]), in_a(cc[cells[2]])
    a1 == a2 && continue
    push!(faces, f)
    lcp[f], lcm[f] = a1 ? (1, 2) : (2, 1)
  end
  plus = BoundaryTriangulation(M, faces, lcp)
  minus = BoundaryTriangulation(M, faces, lcm)
  CompositeTriangulation(Γη, SkeletonTriangulation(plus, minus)), length(faces)
end

@testset "Spike: plate–plate interface skeleton" begin
  s1 = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [2.0, 1.0, 1.0], domain_symbol = :Γ_a)
  s2 = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [4.0, 1.0, 1.0], domain_symbol = :Γ_b)
  tank = G.TankDomain(L = 8.0, W = 4.0, H = 1.0, nx = 8, ny = 8, nz = 2, structure_domains = [s1, s2])
  tr = G.build_triangulations(tank, G.build_model(tank))
  Γa, Γb = tr[:Γ_a], tr[:Γ_b]
  Λ, nfaces = _spike_interface_skeleton(tr[:Γη], c -> c[1] < 4.0)
  dΛ = Measure(Λ, 4)
  n = get_normal_vector(Λ)
  e1, e2 = VectorValue(1.0, 0.0, 0.0), VectorValue(0.0, 1.0, 0.0)

  # geometry: the common edge x = 4, y ∈ [1, 3] (4 facets of Δy = 0.5)
  @test nfaces == 4
  @test sum(∫(1.0)dΛ) ≈ 2.0
  @test sum(∫(n.⁺ ⋅ e1)dΛ) ≈ 2.0            # plus normal points out of plate a
  @test sum(∫(n.⁻ ⋅ e1)dΛ) ≈ -2.0

  Va = TestFESpace(Γa, ReferenceFE(lagrangian, Float64, 2))
  Vb = TestFESpace(Γb, ReferenceFE(lagrangian, Float64, 2))
  f(x) = x[1]^2 + 3x[2] * x[1]              # ∇f = (2x + 3y, 3x)
  ηa, ηb = interpolate(f, Va), interpolate(f, Vb)

  # traces: each side sees its own field, and both agree on the edge
  If = sum(∫(CellField(f, Λ))dΛ)
  @test sum(∫(ηa.⁺)dΛ) ≈ If
  @test sum(∫(ηb.⁻)dΛ) ≈ If
  @test sum(∫((ηa.⁺ - ηb.⁻) * (ηa.⁺ - ηb.⁻))dΛ) < 1e-20

  # gradients (contract before taking the trace, see caveat)
  Igx = sum(∫(CellField(x -> 2x[1] + 3x[2], Λ))dΛ)
  Igy = sum(∫(CellField(x -> 3x[1], Λ))dΛ)
  @test sum(∫((∇(ηa) ⋅ e1).⁺)dΛ) ≈ Igx
  @test sum(∫((∇(ηb) ⋅ e1).⁻)dΛ) ≈ Igx
  @test sum(∫((∇(ηa) ⋅ e2).⁺)dΛ) ≈ Igy
  @test sum(∫((∇(ηb) ⋅ e2).⁻)dΛ) ≈ Igy
  @test sum(∫(∇(ηa).⁺ ⋅ n.⁺)dΛ) ≈ Igx
  @test sum(∫(∇(ηb).⁻ ⋅ n.⁻)dΛ) ≈ -Igx
  @test_broken sum(∫(∇(ηb).⁻ ⋅ e1)dΛ) ≈ Igx   # caveat: minus trace ⋅ constant

  # multi-field assembly across the interface
  Y = MultiFieldFESpace([Va, Vb])
  X = MultiFieldFESpace([TrialFESpace(Va), TrialFESpace(Vb)])
  A = assemble_matrix(((ua, ub), (va, vb)) -> ∫((va.⁺ - vb.⁻) * (ua.⁺ - ub.⁻))dΛ, X, Y)
  @test size(A) == (num_free_dofs(Va) + num_free_dofs(Vb), num_free_dofs(Va) + num_free_dofs(Vb))
  @test norm(A) > 0
  @test norm(A * ones(size(A, 2))) < 1e-12          # continuous (constant) field: no jump
  # slope jump ∇u_a⁺⋅n⁺ + ∇u_b⁻⋅n⁻ vanishes for a globally linear field
  S = assemble_matrix(((ua, ub), (va, vb)) ->
        ∫((∇(va).⁺ ⋅ n.⁺ + ∇(vb).⁻ ⋅ n.⁻) * (∇(ua).⁺ ⋅ n.⁺ + ∇(ub).⁻ ⋅ n.⁻))dΛ, X, Y)
  lin(x) = 2x[1] - x[2]
  xl = [get_free_dof_values(interpolate(lin, Va)); get_free_dof_values(interpolate(lin, Vb))]
  @test norm(S * xl) < 1e-10 * norm(S)
  xk = [get_free_dof_values(interpolate(lin, Va)); zeros(num_free_dofs(Vb))]
  @test norm(S * xk) > 1e-3 * norm(S)               # a kink is penalised
end
