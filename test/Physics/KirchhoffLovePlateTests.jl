using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces

import HydroElasticFEM.Physics as P
import HydroElasticFEM.Geometry as D
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.ParameterHandler as FES

# =========================================================================
# KirchhoffLovePlate weak-form integration tests
#
# All problems are structural-only (no fluid).  The plate is normalised
# with ρ_fluid = 1.0 so that plate.C[1,1,1,1] == D (bending stiffness).
#
# Reference: Timoshenko & Woinowsky-Krieger, "Theory of Plates and Shells",
#   2nd ed., McGraw-Hill.  §29: w_max = 0.00416 q a^4 / D  (ν = 0.3).
# =========================================================================

# -----------------------------------------------------------------------
# Helper: build 2D plate mesh, FE spaces, and IntegrationDomains
# -----------------------------------------------------------------------

function _build_plate_problem(; L, n, order)
  model = CartesianDiscreteModel((0.0, L, 0.0, L), (n, n))
  Ω  = Triangulation(model)
  Λ  = Skeleton(Ω)
  hm = L / n

  dom = D.IntegrationDomains(
    dΓη    = Measure(Ω, 2 * order + 2),
    dΛη    = Measure(Λ, 2 * order + 2),
    n_Λ_η  = get_normal_vector(Λ),
    h_η    = hm,
  )

  return (; model, dom)
end

# -----------------------------------------------------------------------
# Helper: build FE spaces with Dirichlet η = 0 on all four edges
# (simply-supported: zero deflection, zero moment is natural)
# -----------------------------------------------------------------------

function _build_plate_fe_spaces(model, order)
  reffe = ReferenceFE(lagrangian, Float64, order)
  V = TestFESpace(model, reffe;
    conformity   = :H1,
    dirichlet_tags = "boundary",
    vector_type  = Vector{Float64})
  U = TrialFESpace(V, 0.0)
  Y = MultiFieldFESpace([V])
  X = MultiFieldFESpace([U])
  return X, Y
end

# -----------------------------------------------------------------------
# Helper: solve static plate problem and return centre deflection
# -----------------------------------------------------------------------

function _run_kl_plate_simply_supported(; E, nu, h_plate, L, q, n, order)
  prob  = _build_plate_problem(L=L, n=n, order=order)
  X, Y  = _build_plate_fe_spaces(prob.model, order)

  # ρ_fluid = 1.0: no FSI normalisation; C[1,1,1,1] = D (bending stiffness)
  plate = P.KirchhoffLovePlate(
    E            = E,
    ν            = nu,
    hb           = h_plate,
    ρ            = 1.0,
    ρb           = 1.0,
    g            = 0.0,
    ambient_dim  = 2,
    manifold_dim = 2,
    fe           = FES.FESpaceConfig(order=order, vector_type=Vector{Float64}),
  )

  sym  = P.variable_symbol(plate)
  fmap = Dict(sym => 1)

  a((u,), (v,)) = P.stiffness(plate, prob.dom,
                               FO.FieldMap((u,), fmap),
                               FO.FieldMap((v,), fmap))

  src(x)  = q
  l((v,)) = P.rhs(plate, prob.dom,
                  FO.FieldMap((src,), fmap),
                  FO.FieldMap((v,), fmap))

  op  = AffineFEOperator(a, l, X, Y)
  uh  = solve(LUSolver(), op)
  uh1 = uh[1]   # scalar deflection field

  return uh1(Point(L / 2, L / 2))
end

# =========================================================================
# Tests
# =========================================================================

@testset "KirchhoffLovePlate — struct and traits" begin
  ρw   = 1025.0
  plate = P.KirchhoffLovePlate(
    E = 11.9e9, ν = 0.13, hb = 2.0, ρ = ρw,
  )

  @test plate isa P.Structure
  @test P.variable_symbol(plate) == :η
  @test plate.space_domain_symbol == :Γη
  @test P.ambient_dimension(plate) == 3
  @test P.manifold_dimension(plate) == 2

  # Stiffness-proportional damping is not implemented for the plate
  @test P.has_damping_form(plate) == false
  # Mass and stiffness forms are present
  @test P.has_mass_form(plate) == true
  @test P.has_stiffness_form(plate) == true
  @test P.has_rhs_form(plate) == true
end

@testset "KirchhoffLovePlate — constitutive tensor symmetry" begin
  # 2D ambient tensor: build_kl_tensor(ambient=2, manifold=2)
  C2 = P.build_kl_tensor(2, 2, 10.92e6, 0.3, 0.01, 1.0)

  for i in 1:2, j in 1:2, k in 1:2, l in 1:2
    @test C2[i, j, k, l] ≈ C2[k, l, i, j]  atol=1e-10
    @test C2[i, j, k, l] ≈ C2[j, i, k, l]  atol=1e-10
    @test C2[i, j, k, l] ≈ C2[i, j, l, k]  atol=1e-10
  end

  # D/ρ at (1,1,1,1) — bending stiffness / ρ
  E = 10.92e6;  nu = 0.3;  h = 0.01;  rho = 1.0
  D_rho = E * h^3 / (12 * (1 - nu^2) * rho)
  @test C2[1, 1, 1, 1] ≈ D_rho  atol=1e-10 * D_rho

  # 3D ambient tensor via build_KL_tensor
  # Checks major symmetry C[i,j,k,l] = C[k,l,i,j] for all i,j,k,l ∈ 1:3
  C3 = P.build_KL_tensor(10.92e6, 0.3, 0.01, 1025.0)
  for i in 1:3, j in 1:3, k in 1:3, l in 1:3
    @test C3[i, j, k, l] ≈ C3[k, l, i, j]  atol=1e-10
  end
end

@testset "KirchhoffLovePlate — simply-supported square, uniform load" begin
  # Timoshenko §29: w_max = 0.00416 · q · L^4 / D  (for ν = 0.3)
  L       = 1.0
  q       = 1.0
  E       = 10.92e6
  nu      = 0.3
  h_plate = 0.01
  D       = E * h_plate^3 / (12 * (1 - nu^2))
  w_ref   = 0.00416 * q * L^4 / D

  w_h = _run_kl_plate_simply_supported(
    E=E, nu=nu, h_plate=h_plate, L=L, q=q, n=10, order=2,
  )

  # 5 % relative tolerance: coarse 10×10 mesh with quadratic C/DG
  # converges to the analytical value but is under-resolved at this resolution.
  @test abs(w_h - w_ref) / w_ref < 0.05
end
