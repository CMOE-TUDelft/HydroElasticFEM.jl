using Test

using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces

import HydroElasticFEM.PhysicsCore.Entities as E
import HydroElasticFEM.Geometry as D
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.Simulation.FESpaceAssembly as FEA
import HydroElasticFEM.PhysicsCore.FESpaces as FES

# =========================================================================
# EulerBernoulliBeam weak form integration tests
#
# Solve static beam problems (stiffness = rhs) and verify against
# analytical Euler-Bernoulli solutions.
#
# PDE:  EIᵨ·u'''' = q
# =========================================================================

@testset "EulerBernoulliBeam weak forms" begin

  # -----------------------------------------------------------------------
  # Helper: build 1D beam mesh, FE spaces, and IntegrationDomains
  # -----------------------------------------------------------------------

  function build_beam_problem(; L, nel, order)
    model = CartesianDiscreteModel((0, L), (nel,))
    Ω  = Triangulation(model)
    Λ  = Skeleton(Ω)
    Λb = Boundary(model, tags="boundary")

    h = L / nel

    dom = D.IntegrationDomains(
      dΓ_s   = Measure(Ω, 2 * order + 2),
      dΛ_s   = Measure(Λ, 2 * order + 2),
      n_Λ_s  = get_normal_vector(Λ),
      h_s    = h,
      dΛ_sb  = Measure(Λb, 2 * order + 2),
      n_Λ_sb = get_normal_vector(Λb),
    )

    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(model, reffe, dirichlet_tags="boundary",
                    vector_type=Vector{Float64})
    U = TrialFESpace(V, 0.0)

    # Wrap in MultiFieldFESpace so Gridap passes 1-tuples to closures,
    # allowing FieldMap indexing to work.
    Y = MultiFieldFESpace([V])
    X = MultiFieldFESpace([U])

    return (; model, Ω, dom, X, Y, h)
  end

  # -----------------------------------------------------------------------
  # Helper: solve and extract max deflection
  # -----------------------------------------------------------------------

  function solve_beam(beam, nel, order, q)
    prob = build_beam_problem(L=beam.L, nel=nel, order=order)
    sym  = E.variable_symbol(beam)
    fmap = Dict(sym => 1)

    a((u,), (v,)) = E.stiffness(beam, prob.dom,
                               FO.FieldMap((u,), fmap), FO.FieldMap((v,), fmap))

    src(x) = q
    l((v,)) = E.rhs(beam, prob.dom,
                   FO.FieldMap((src,), fmap), FO.FieldMap((v,), fmap))

    op = AffineFEOperator(a, l, prob.X, prob.Y)
    uh = solve(LUSolver(), op)

    # Extract single-field solution from multi-field wrapper
    uh1 = uh[1]

    return (; uh=uh1, Ω=prob.Ω, h=prob.h)
  end

  # -----------------------------------------------------------------------
  # Test 1: Simply-supported beam (FreeBoundary), uniform load
  #
  #   Analytical: w(x) = q/(24·EIᵨ) · x·(L³ - 2L·x² + x³)
  #               w_max = 5·q·L⁴ / (384·EIᵨ)   at x = L/2
  # -----------------------------------------------------------------------

  @testset "Simply-supported beam — uniform load" begin
    L   = 1.0
    EIᵨ = 100.0
    q   = 1.0

    beam = E.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
                              fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))

    w_exact_max = 5 * q * L^4 / (384 * EIᵨ)

    # Refine mesh to get close to analytical solution
    sol = solve_beam(beam, 40, 2, q)

    # Evaluate at midpoint
    x_mid = Point(L / 2)
    w_mid = sol.uh(x_mid)

    @test w_mid ≈ w_exact_max  rtol=1e-2

    # Boundary conditions: u(0) = u(L) = 0
    @test abs(sol.uh(Point(0.0))) < 1e-12
    @test abs(sol.uh(Point(L)))   < 1e-12

    # Symmetry: w(L/4) ≈ w(3L/4)
    @test sol.uh(Point(L / 4)) ≈ sol.uh(Point(3L / 4))  rtol=1e-6
  end

  # -----------------------------------------------------------------------
  # Test 2: Clamped-clamped beam (FixedBoundary), uniform load
  #
  #   Analytical: w_max = q·L⁴ / (384·EIᵨ)   at x = L/2
  # -----------------------------------------------------------------------

  @testset "Clamped-clamped beam — uniform load" begin
    L   = 1.0
    EIᵨ = 100.0
    q   = 1.0

    beam = E.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
                              fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))

    w_exact_max = q * L^4 / (384 * EIᵨ)

    sol = solve_beam(beam, 40, 2, q)

    x_mid = Point(L / 2)
    w_mid = sol.uh(x_mid)

    @test w_mid ≈ w_exact_max  rtol=1e-2

    # Boundary conditions: u(0) = u(L) = 0
    @test abs(sol.uh(Point(0.0))) < 1e-12
    @test abs(sol.uh(Point(L)))   < 1e-12

    # Symmetry
    @test sol.uh(Point(L / 4)) ≈ sol.uh(Point(3L / 4))  rtol=1e-6
  end

  # -----------------------------------------------------------------------
  # Test 3: Scaling — deflection scales as L⁴/EIᵨ
  #
  #   Doubling L should multiply w_max by 16.
  #   Doubling EIᵨ should halve w_max.
  # -----------------------------------------------------------------------

  @testset "Deflection scaling with L and EIᵨ" begin
    q = 1.0

    beam1 = E.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=100.0,
                               fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))
    beam2 = E.EulerBernoulliBeam(L=2.0, mᵨ=1.0, EIᵨ=100.0,
                               fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))
    beam3 = E.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=200.0,
                               fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))

    sol1 = solve_beam(beam1, 40, 2, q)
    sol2 = solve_beam(beam2, 80, 2, q)
    sol3 = solve_beam(beam3, 40, 2, q)

    w1 = sol1.uh(Point(0.5))
    w2 = sol2.uh(Point(1.0))
    w3 = sol3.uh(Point(0.5))

    # w ∝ L⁴ ⟹ w2/w1 ≈ 16
    @test w2 / w1 ≈ 16.0  rtol=1e-2

    # w ∝ 1/EIᵨ ⟹ w3/w1 ≈ 0.5
    @test w3 / w1 ≈ 0.5  rtol=1e-2
  end

  # -----------------------------------------------------------------------
  # Test 4: Mesh convergence — error decreases with refinement
  # -----------------------------------------------------------------------

  @testset "Mesh convergence" begin
    L   = 1.0
    EIᵨ = 100.0
    q   = 1.0

    beam = E.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
                              fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))
    w_exact = 5 * q * L^4 / (384 * EIᵨ)

    errors = Float64[]
    nels   = [10, 20, 40]

    for nel in nels
      sol = solve_beam(beam, nel, 2, q)
      w_mid = sol.uh(Point(L / 2))
      push!(errors, abs(w_mid - w_exact))
    end

    # Each doubling of elements should reduce error
    @test errors[2] < errors[1]
    @test errors[3] < errors[2]

    # Finest mesh should be within 1% of analytical
    @test errors[end] / abs(w_exact) < 0.01
  end

  # -----------------------------------------------------------------------
  # Test 5: Zero load produces zero deflection
  # -----------------------------------------------------------------------

  @testset "Zero load — zero deflection" begin
    beam = E.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=100.0, g=0.0,
                              fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))

    sol = solve_beam(beam, 20, 2, 0.0)
    w_mid = sol.uh(Point(0.5))

    @test abs(w_mid) < 1e-14
  end

end
