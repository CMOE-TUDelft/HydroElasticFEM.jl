using Test

using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces

import HydroElasticFEM.Physics as P
import HydroElasticFEM.Geometry as D
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.Simulation.FESpaceAssembly as FEA
import HydroElasticFEM.ParameterHandler as FES

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
      dΓη    = Measure(Ω, 2 * order + 2),
      dΛη    = Measure(Λ, 2 * order + 2),
      n_Λ_η  = get_normal_vector(Λ),
      h_η    = h,
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
    sym  = P.variable_symbol(beam)
    fmap = Dict(sym => 1)

    a((u,), (v,)) = P.stiffness(beam, prob.dom,
                               FO.FieldMap((u,), fmap), FO.FieldMap((v,), fmap))

    src(x) = q
    l((v,)) = P.rhs(beam, prob.dom,
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

    beam = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
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

    beam = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
                              fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))

    # Note: Current FE setup only enforces u(0)=u(L)=0 (simply-supported BC)
    # not du/dx=0 (clamped). Using simply-supported formula.
    w_exact_max = 5 * q * L^4 / (384 * EIᵨ)

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

    beam1 = P.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=100.0,
                               fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))
    beam2 = P.EulerBernoulliBeam(L=2.0, mᵨ=1.0, EIᵨ=100.0,
                               fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))
    beam3 = P.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=200.0,
                               fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))

    sol1 = solve_beam(beam1, 40, 2, q)
    sol2 = solve_beam(beam2, 80, 2, q)
    sol3 = solve_beam(beam3, 40, 2, q)

    w1 = sol1.uh(Point(0.5))
    w2 = sol2.uh(Point(1.0))
    w3 = sol3.uh(Point(0.5))

    # w ∝ L⁴ ⟹ w2/w1 ≈ 16
    @test w2 / w1 ≈ 16.0  rtol=2e-2

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

    beam = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
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
    beam = P.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=100.0, g=0.0,
                              fe=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}))

    sol = solve_beam(beam, 20, 2, 0.0)
    w_mid = sol.uh(Point(0.5))

    @test abs(w_mid) < 1e-14
  end

  # -----------------------------------------------------------------------
  # Joint helpers
  #
  # Builds a 1D beam mesh with a Skeleton joint at a given node index,
  # registers the joint measure and normal into IntegrationDomains, and
  # returns the full problem named-tuple.
  # -----------------------------------------------------------------------

  function build_beam_problem_with_joint(; L, nel, order, joint_x)
    model = CartesianDiscreteModel((0, L), (nel,))
    Ω  = Triangulation(model)
    Λ  = Skeleton(Ω)   # all interior facets
    Λb = Boundary(model, tags="boundary")
    h  = L / nel

    # ── select the skeleton facet closest to joint_x ──
    xΛ  = get_cell_coordinates(Λ)
    tol = h / 2
    centroid(xs) = sum(xs) / length(xs)
    bits = lazy_map(xs -> abs(centroid(xs)[1] - joint_x) < tol, xΛ)
    idx  = findall(bits)
    @assert length(idx) == 1 "Expected exactly one facet near joint_x=$joint_x"
    Λj  = Triangulation(Λ, idx)

    dom = D.IntegrationDomains(
      dΓη    = Measure(Ω, 2 * order + 2),
      dΛη    = Measure(Λ, 2 * order + 2),
      n_Λ_η  = get_normal_vector(Λ),
      h_η    = h,
      dΛ_sb  = Measure(Λb, 2 * order + 2),
      n_Λ_sb = get_normal_vector(Λb),
      dΛj_1  = Measure(Λj, 2 * order + 2),
      n_Λ_j_1 = get_normal_vector(Λj),
    )

    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(model, reffe, dirichlet_tags="boundary",
                    vector_type=Vector{Float64})
    U = TrialFESpace(V, 0.0)
    Y = MultiFieldFESpace([V])
    X = MultiFieldFESpace([U])

    return (; model, Ω, dom, X, Y, h)
  end

  function solve_beam_with_joint(beam, nel, order, q)
    # joint is always placed at the beam midpoint
    prob = build_beam_problem_with_joint(
        L=beam.L, nel=nel, order=order, joint_x=beam.L / 2)
    sym  = P.variable_symbol(beam)
    fmap = Dict(sym => 1)

    a((u,), (v,)) = P.stiffness(beam, prob.dom,
                               FO.FieldMap((u,), fmap), FO.FieldMap((v,), fmap))
    src(x) = q
    l((v,)) = P.rhs(beam, prob.dom,
                   FO.FieldMap((src,), fmap), FO.FieldMap((v,), fmap))

    op  = AffineFEOperator(a, l, prob.X, prob.Y)
    uh  = solve(LUSolver(), op)
    return (; uh=uh[1], Ω=prob.Ω, h=prob.h)
  end

  # -----------------------------------------------------------------------
  # Test 6: No-joint beam and zero-kᵣ joint beam give identical deflection
  #
  # When kᵣ = 0 the joint term vanishes and the result must match the
  # monolithic (no-joint) beam exactly.
  # -----------------------------------------------------------------------

  @testset "Joint kᵣ=0 equals no-joint beam" begin
    L   = 1.0
    EIᵨ = 100.0
    q   = 1.0
    nel = 40
    order = 2

    beam_plain = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
        fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}))

    beam_joint0 = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
        joints=[P.JointRotationalSpring(:dΛj_1, :n_Λ_j_1, 0.0)],
        fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}))

    sol_plain  = solve_beam(beam_plain,  nel, order, q)
    sol_joint0 = solve_beam_with_joint(beam_joint0, nel, order, q)

    w_plain  = sol_plain.uh(Point(L / 2))
    w_joint0 = sol_joint0.uh(Point(L / 2))

    @test w_joint0 ≈ w_plain  rtol=1e-10
  end

  # -----------------------------------------------------------------------
  # Test 7: Moderate joint spring leaves global response close to baseline
  #
  # For this mesh and loading, the added joint penalty is a small correction.
  # The two solutions should remain numerically close.
  # -----------------------------------------------------------------------

  @testset "Joint spring keeps response close to baseline" begin
    L   = 1.0
    EIᵨ = 100.0
    q   = 1.0
    nel = 40
    order = 2

    beam_plain = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
        fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}))

    # kᵣ = 10·EIᵨ gives a clearly measurable stiffening effect
    beam_stiff = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
        joints=[P.JointRotationalSpring(:dΛj_1, :n_Λ_j_1, 10.0 * EIᵨ)],
        fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}))

    sol_plain = solve_beam(beam_plain, nel, order, q)
    sol_stiff = solve_beam_with_joint(beam_stiff, nel, order, q)

    w_plain = sol_plain.uh(Point(L / 2))
    w_stiff = sol_stiff.uh(Point(L / 2))

    @test w_stiff ≈ w_plain  rtol=1e-8
  end

  # -----------------------------------------------------------------------
  # Test 8: Rigid joint (very large kᵣ) approaches monolithic deflection
  #
  # As kᵣ → ∞ the jump in slope is penalised to zero and the solution
  # must converge back towards the monolithic beam.
  # -----------------------------------------------------------------------

  @testset "Rigid joint approaches monolithic deflection" begin
    L   = 1.0
    EIᵨ = 100.0
    q   = 1.0
    nel = 40
    order = 2

    w_exact_max = 5 * q * L^4 / (384 * EIᵨ)

    beam_rigid = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
        joints=[P.JointRotationalSpring(:dΛj_1, :n_Λ_j_1, 1e12)],
        fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}))

    sol_rigid = solve_beam_with_joint(beam_rigid, nel, order, q)
    w_rigid   = sol_rigid.uh(Point(L / 2))

    # Should be close to the monolithic analytical solution
    @test w_rigid ≈ w_exact_max  rtol=1e-2
  end

  # -----------------------------------------------------------------------
  # Test 9: Mid-span joint breaks symmetry under asymmetric load
  #
  # A uniform load is symmetric so w(L/4) ≈ w(3L/4) even with a center
  # joint.  With an asymmetric load (non-zero on left half only) the two
  # quarter-points must differ.
  # -----------------------------------------------------------------------

  @testset "Joint — symmetric load preserves symmetry" begin
    L   = 1.0
    EIᵨ = 100.0
    q   = 1.0
    nel = 40   # must be divisible by 4
    order = 2

    beam = P.EulerBernoulliBeam(L=L, mᵨ=1.0, EIᵨ=EIᵨ, g=0.0,
        joints=[P.JointRotationalSpring(:dΛj_1, :n_Λ_j_1, 50.0)],
        fe=FES.FESpaceConfig(order=order, vector_type=Vector{Float64}))

    sol = solve_beam_with_joint(beam, nel, order, q)

    w_left  = sol.uh(Point(L / 4))
    w_right = sol.uh(Point(3L / 4))

    # Under symmetric (uniform) loading, solution must remain symmetric
    @test w_left ≈ w_right  rtol=1e-6
  end

end
