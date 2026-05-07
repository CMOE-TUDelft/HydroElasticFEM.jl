using Test
using Gridap
using LinearAlgebra

import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.ParameterHandler as PH

# =========================================================================
# 3D sloshing tests
#
# Validates TankDomain{3} geometry, 3D PotentialFlow + FreeSurface assembly,
# and the fundamental sloshing eigenfrequency against the analytical
# dispersion relation:
#
#   ω² = g · k · tanh(k · H),   k = π / L   (mode m=1, n=0)
#
# All units SI: meters, kg, seconds, radians, Pascals.
# =========================================================================

@testset "TankDomain{3} geometry" begin

  tank = G.TankDomain(L = 10.0, W = 5.0, H = 1.0, nx = 4, ny = 2, nz = 2)

  @test G.ambient_dimension(tank) == 3
  @test G.manifold_dimension(tank) == 3

  btags = G.boundary_tags(tank)
  for tag in G.STANDARD_TAGS
    @test haskey(btags, tag)
  end
  @test haskey(btags, "lateral_walls")

  model  = G.build_model(tank)
  trians = G.build_triangulations(tank, model)

  # nx*ny*nz hexahedral cells in the interior
  @test num_cells(trians[:Ω]) == 4 * 2 * 2

  # top and bottom faces: nx*ny quads each
  @test num_cells(trians[:Γfs])  == 4 * 2
  @test num_cells(trians[:Γbot]) == 4 * 2

  # inlet / outlet faces: ny*nz quads each
  @test num_cells(trians[:Γin])  == 2 * 2
  @test num_cells(trians[:Γout]) == 2 * 2

  # lateral walls: 2 faces × nx*nz quads
  @test num_cells(trians[:Γlateral]) == 2 * 4 * 2

  # structure surface is empty (no structures in TankDomain{3})
  @test num_cells(trians[:Γη]) == 0
end

@testset "PotentialFlow + FreeSurface 3D assembly" begin

  tank = G.TankDomain(L = 10.0, W = 5.0, H = 1.0, nx = 4, ny = 2, nz = 2)

  fe3 = PH.FESpaceConfig(order = 1, vector_type = Vector{ComplexF64})
  pf  = P.PotentialFlow(dim = 3, fe = fe3)
  fs  = P.FreeSurface(dim = 3, fe = fe3)

  @test P.ambient_dimension(pf) == 3
  @test P.ambient_dimension(fs) == 3

  # Dimension mismatch must fail early in build_problem.
  @test_throws ErrorException begin
    pf2d = P.PotentialFlow(fe = fe3)  # default dim=2
    fs3d = P.FreeSurface(dim = 3, fe = fe3)
    SM.build_problem(tank, P.PhysicsParameters[pf2d, fs3d], SM.FreqDomainConfig(ω = 1.5))
  end

  # Solve at ω = 1.5 rad/s (non-resonant for this coarse mesh)
  config  = SM.FreqDomainConfig(ω = 1.5)
  problem = SM.build_problem(tank, P.PhysicsParameters[pf, fs], config)

  result = @test_nowarn SM.simulate(problem)
  result = SM.simulate(problem)

  # Two-field solution (ϕ, κ)
  @test length(result.solution) == 2

  # Solution is finite
  dom = SM.get_integration_domains(problem)
  dΩ  = dom[:dΩ]
  ϕₕ  = result.solution[1]
  l2  = real(sum(∫(ϕₕ * conj(ϕₕ))dΩ))
  @test isfinite(l2)
end

@testset "3D sloshing eigenfrequency within 0.5%" begin

  # Physical constants (SI)
  g = 9.81
  L = 10.0; W = 5.0; H = 1.0

  # Analytical fundamental sloshing frequency — mode (1,0) in x
  # ω² = g · k · tanh(k · H),  k = π / L
  k10    = π / L
  ω_anal = sqrt(g * k10 * tanh(k10 * H))

  # Fine enough mesh for ≤ 0.5 % FEM eigenvalue error with P1 elements.
  # Error estimate (1D analysis): Δω/ω ≈ (k·h_x)² / 12 ≈ 0.2 % with h_x = 0.5 m.
  tank   = G.TankDomain(L = L, W = W, H = H, nx = 20, ny = 2, nz = 4)
  model  = G.build_model(tank)
  trians = G.build_triangulations(tank, model)

  dΩ     = Measure(trians[:Ω],   2)
  dΓtop  = Measure(trians[:Γfs], 2)

  # P1 Lagrange test/trial spaces (pure Neumann problem; no Dirichlet BCs)
  reffe = ReferenceFE(lagrangian, Float64, 1)
  V = TestFESpace(model, reffe; conformity = :H1)
  U = TrialFESpace(V)

  # Stiffness: K_ij = ∫ ∇ψᵢ · ∇ψⱼ dΩ
  a_k(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ

  # Free-surface mass: M_ij = ∫ ψᵢ · ψⱼ dΓtop
  a_m(u, v) = ∫(v * u)dΓtop

  K = assemble_matrix(a_k, U, V)
  M = assemble_matrix(a_m, U, V)

  # K is PSD (null space = constants); M is PSD (zero rows for non-top DOFs).
  # Regularise M → M_reg = M + ε·I to obtain a positive-definite denominator.
  # Sloshing eigenvalues:  λ ≈ k·tanh(kH)  (small, O(0.1))
  # Regularisation artefacts:  λ_art ≈ K_diag / ε  (very large, O(1/ε))
  # Choose ε small enough that artefacts are above 1e6 and sloshing modes
  # are unaffected to machine precision.
  n   = size(K, 1)
  ε   = 1.0e-10
  K_d = Matrix(K)
  M_r = Matrix(M) + ε * Matrix(I, n, n)

  vals = eigvals(Symmetric(K_d), Symmetric(M_r))

  # Retain only genuine sloshing eigenvalues:
  # - strip null-space zeros  (λ ≤ 1e-5)
  # - strip regularisation artefacts  (λ ≥ 1e6)
  λ_sloshing = sort([v for v in real.(vals) if 1.0e-5 < v < 1.0e6])

  @test !isempty(λ_sloshing)

  ω_fem = sqrt(g * λ_sloshing[1])

  # Quantitative tolerance: FEM eigenfrequency within 0.5 % of analytical value
  @test abs(ω_fem - ω_anal) / ω_anal < 0.005
end
