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

@testset "KirchhoffLovePlate — 1D manifold (beam reduction)" begin
  # manifold_dim=1: only C[1,1,1,1] = E*I/ρ is active; all other entries 0.
  # Reproduces the Euler-Bernoulli bending rigidity per unit fluid density.
  E = 10.92e6;  h = 0.01;  ρ = 1.0
  C1 = P.build_kl_tensor(2, 1, E, 0.3, h, ρ)
  @test C1[1, 1, 1, 1] ≈ E * h^3 / (12 * ρ)  rtol=1e-12
  @test C1[1, 2, 1, 2] == 0.0
  @test C1[2, 2, 2, 2] == 0.0
end

@testset "KirchhoffLovePlate — build_kl_tensor argument validation" begin
  @test_throws ErrorException P.build_kl_tensor(4, 2, 1e9, 0.3, 0.01, 1.0)  # bad ambient_dim
  @test_throws ErrorException P.build_kl_tensor(2, 0, 1e9, 0.3, 0.01, 1.0)  # bad manifold_dim
  @test_throws ErrorException P.build_kl_tensor(2, 3, 1e9, 0.3, 0.01, 1.0)  # manifold > ambient
end

@testset "KirchhoffLovePlate — constitutive tensor ρ scaling" begin
  # Doubling ρ halves every tensor component.
  E = 10.92e6;  ν = 0.3;  h = 0.01
  C1 = P.build_kl_tensor(2, 2, E, ν, h, 1.0)
  C2 = P.build_kl_tensor(2, 2, E, ν, h, 2.0)
  for i in 1:2, j in 1:2, k in 1:2, l in 1:2
    @test C1[i, j, k, l] ≈ 2 * C2[i, j, k, l]  atol=1e-12
  end
end

@testset "KirchhoffLovePlate — check_major_symmetry" begin
  C = P.build_KL_tensor(10.92e6, 0.3, 0.01, 1025.0)
  @test P.check_major_symmetry(C; dim=3)
  @test P.check_major_symmetry(C; dim=2)
end

@testset "KirchhoffLovePlate — equivalent_beam_rigidity" begin
  E = 10.92e6;  ν = 0.3;  h = 0.01;  ρ = 1025.0
  plate   = P.KirchhoffLovePlate(E=E, ν=ν, hb=h, ρ=ρ)
  D_rho   = E * h^3 / (12 * (1 - ν^2) * ρ)
  @test P.equivalent_beam_rigidity(plate)            ≈ D_rho      rtol=1e-12
  @test P.equivalent_beam_rigidity(plate; width=3.0) ≈ 3 * D_rho  rtol=1e-12
end

@testset "KirchhoffLovePlate — mass form integral" begin
  # Assemble l(v) = mass(plate, dom, ηₜₜ=1, y) on a free (no Dirichlet) space.
  # By Lagrange partition of unity: Σᵢ bᵢ = ∫ m_ρ · 1 dΓ = m_ρ · L² (exact).
  L = 1.0;  n = 4;  order = 2
  E = 10.92e6;  ν = 0.3;  h = 0.01;  ρ = 1.0;  ρb = 500.0
  m_ρ = ρb * h / ρ

  prob  = _build_plate_problem(L=L, n=n, order=order)
  reffe = ReferenceFE(lagrangian, Float64, order)
  V     = TestFESpace(prob.model, reffe;
                      conformity=:H1, vector_type=Vector{Float64})
  Y     = MultiFieldFESpace([V])

  plate = P.KirchhoffLovePlate(
    E=E, ν=ν, hb=h, ρ=ρ, ρb=ρb, g=0.0,
    ambient_dim=2, manifold_dim=2,
    fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}),
  )

  sym  = P.variable_symbol(plate)
  fmap = Dict(sym => 1)
  one_fn(x) = 1.0

  l((v,)) = P.mass(plate, prob.dom,
                   FO.FieldMap((one_fn,), fmap),
                   FO.FieldMap((v,), fmap))

  b = assemble_vector(l, Y)
  @test sum(b) ≈ m_ρ * L^2  rtol=1e-10
end

@testset "KirchhoffLovePlate — print_parameters" begin
  plate = P.KirchhoffLovePlate(E=11.9e9, ν=0.13, hb=2.0, ρ=1025.0)
  @test_nowarn P.print_parameters(plate)
end

# =========================================================================
# Hessian SIPG form — targeted tests
# =========================================================================

@testset "KirchhoffLovePlate — Hessian SIPG: stiffness matrix symmetry" begin
  # The Hessian SIPG form a(η,v) is symmetric by construction.
  # Verify ||K - Kᵀ||∞ / ||K||∞ < 1e-10 on a small unconstrained mesh.
  # A failure here means the skeleton terms mean(C⊙∇∇η)⋅n are mis-assembled.
  L = 1.0;  n = 3;  order = 2
  E = 10.92e6;  ν = 0.3;  h_plate = 0.01

  prob  = _build_plate_problem(L=L, n=n, order=order)
  reffe = ReferenceFE(lagrangian, Float64, order)
  V  = TestFESpace(prob.model, reffe; conformity=:H1,
                   vector_type=Vector{Float64})
  U  = TrialFESpace(V)
  Y  = MultiFieldFESpace([V])
  X  = MultiFieldFESpace([U])

  plate = P.KirchhoffLovePlate(
    E=E, ν=ν, hb=h_plate, ρ=1.0, ρb=1.0, g=0.0,
    ambient_dim=2, manifold_dim=2,
    fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}),
  )
  sym  = P.variable_symbol(plate)
  fmap = Dict(sym => 1)

  a((u,), (v,)) = P.stiffness(plate, prob.dom,
                               FO.FieldMap((u,), fmap),
                               FO.FieldMap((v,), fmap))

  K  = assemble_matrix(a, X, Y)
  Kd = Matrix(K)
  sym_err = maximum(abs.(Kd .- Kd')) / maximum(abs.(Kd))
  @test sym_err < 1e-10
end

@testset "KirchhoffLovePlate — Hessian SIPG: ν=0 simply-supported" begin
  # For ν=0 the Poisson coupling coefficient λ=0; only the shear term μ is
  # active in C.  The SS plate solution is governed by D=Eh³/12 (ν=0) and
  # the Timoshenko reference still applies with the same α≈0.00416.
  L       = 1.0
  q       = 1.0
  E       = 10.92e6
  h_plate = 0.01
  D       = E * h_plate^3 / 12    # ν=0 → D = Eh³/12
  w_ref   = 0.00416 * q * L^4 / D

  w_h = _run_kl_plate_simply_supported(
    E=E, nu=0.0, h_plate=h_plate, L=L, q=q, n=10, order=2,
  )
  @test abs(w_h - w_ref) / w_ref < 0.05
end

@testset "KirchhoffLovePlate — Hessian SIPG: p-refinement" begin
  # On the same coarse 5×5 mesh, order=3 should give smaller error than order=2.
  L = 1.0;  q = 1.0
  E = 10.92e6;  ν = 0.3;  h_plate = 0.01
  D   = E * h_plate^3 / (12 * (1 - ν^2))
  w_ref = 0.00416 * q * L^4 / D

  w2 = _run_kl_plate_simply_supported(E=E, nu=ν, h_plate=h_plate,
                                       L=L, q=q, n=5, order=2)
  w3 = _run_kl_plate_simply_supported(E=E, nu=ν, h_plate=h_plate,
                                       L=L, q=q, n=5, order=3)

  err2 = abs(w2 - w_ref) / w_ref
  err3 = abs(w3 - w_ref) / w_ref

  @test err3 < err2           # order=3 is more accurate than order=2
  @test err3 < 0.05           # order=3 on n=5 mesh is within 5%
end
