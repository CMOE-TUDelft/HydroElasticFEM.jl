using Test
using Gridap
using Gridap.ODEs
using LinearAlgebra
import HydroElasticFEM
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.AssemblyContexts as AC

# Residual of the wrapped operator vs Gridap's re-assembled residual at state us
function _residual_pair(op_alg, t, us)
  odeop = Gridap.FESpaces.get_algebraic_operator(op_alg)
  cache = allocate_odeopcache(odeop, t, us)
  update_odeopcache!(cache, odeop, t)
  r_alg = Gridap.Algebra.allocate_residual(odeop, t, us, cache)
  r_std = similar(r_alg)
  Gridap.Algebra.residual!(r_alg, odeop, t, us, cache)
  Gridap.Algebra.residual!(r_std, odeop.inner, t, us, cache)
  r_alg, r_std, cache
end

@testset "Algebraic linear residual" begin
  # Membrane between two damping zones, inlet forcing: all three forms and a
  # time-dependent RHS are active.
  H0, Lm, Ld = 10.0, 10.0, 10.0
  LΩ = 2 * Ld + 3 * Lm
  βₕ, Δt = 0.5, 0.1
  αₕ = 0.5 / (0.25 * Δt) / 9.81 * (1.0 - βₕ) / βₕ
  η_in(x, t) = 0.05 * cos(0.2 * x[1] - t)
  vz_in(x, t) = 0.05 * sin(0.2 * x[1] - t)
  tank = G.TankDomain(L = LΩ, H = H0, nx = 20, ny = 4,
    structure_domains = [G.StructureDomain(L = Lm, x₀ = [Ld + Lm / 2, H0], domain_symbol = :Γm)],
    damping_zones = [G.DampingZone(L = Ld, x₀ = [0.0, H0], domain_symbol = :Γd_in),
                     G.DampingZone(L = Ld, x₀ = [LΩ - Ld, H0], domain_symbol = :Γd_out)])
  fe = PH.FESpaceConfig(vector_type = Vector{Float64})
  fluid = P.PotentialFlow(fe = fe, boundary_conditions = [
    P.PrescribedInletPotentialBC(domain = :dΓin, quantity = :traction,
                                 forcing = (t -> (x -> 0.1 * cos(0.2 * x[1] - t)))),
    P.DampingZoneBC(domain = :dΓd_1, μ₁ = (x -> 2.0), μ₂ = (x -> 1.0),
                    η_in = (t -> (x -> η_in(x, t))), vz_in = (t -> (x -> vz_in(x, t)))),
    P.DampingZoneBC(domain = :dΓd_2, μ₁ = (x -> 2.0), μ₂ = (x -> 1.0),
                    η_in = (x -> 0.0), vz_in = (x -> 0.0))])
  physics = P.PhysicsParameters[fluid, P.FreeSurface(βₕ = βₕ, fe = fe),
                                P.Membrane(L = Lm, mᵨ = 0.9, Tᵨ = 98.1, fe = fe)]
  cfg = PH.TimeDomainConfig(t₀ = 0.0, tf = 1.0)
  tc = PH.TimeConfig(Δt = Δt, tf = 1.0, αₕ = αₕ, u0 = zeros(3), u0t = zeros(3), u0tt = zeros(3))
  prob = SM.build_problem(tank, physics, cfg; tconfig = tc)
  op_alg = SM.get_fe_operator(prob)
  @test op_alg isa FO.AlgebraicLinearTFEOperator
  @test op_alg isa Gridap.ODEs.TransientFEOperator

  X, Y, fmap = SM.get_trial_fe_space(prob), SM.get_test_fe_space(prob), SM.get_field_map(prob)
  ctx = SM.get_assembly_context(prob)
  op_std = FO.build_time_fe_operator(physics, ctx, fmap, X, Y; algebraic_residual = false)
  @test !(op_std isa FO.AlgebraicLinearTFEOperator)

  @testset "residual equals Gridap's re-assembled residual" begin
    n = num_free_dofs(X(0.0))
    for t in (0.0, 0.37, 1.3)
      us = (randn(n), randn(n), randn(n))
      r_alg, r_std, cache = _residual_pair(op_alg, t, us)
      @test FO._algebraic_residual_applies(cache)
      @test norm(r_alg - r_std) <= 1e-12 * norm(r_std)
    end
  end

  @testset "time integration matches the stock operator" begin
    X0 = X(0.0)
    u0 = (zero(X0), zero(X0), zero(X0))
    solver = GeneralizedAlpha2(LUSolver(), Δt, 0.8)
    s_alg = [copy(get_free_dof_values(u)) for (_, u) in solve(solver, op_alg, 0.0, 1.0, u0)]
    s_std = [copy(get_free_dof_values(u)) for (_, u) in solve(solver, op_std, 0.0, 1.0, u0)]
    @test length(s_alg) == length(s_std) == 10
    @test maximum(abs, s_std[end]) > 0
    for (a, b) in zip(s_alg, s_std)
      @test norm(a - b) <= 1e-10 * norm(b)
    end
  end

  @testset "falls back when Dirichlet values are non-zero" begin
    model = CartesianDiscreteModel((0, 1, 0, 1), (4, 4))
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 2)
    V = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1); dirichlet_tags = "boundary")
    U = TransientTrialFESpace(V, t -> (x -> 1.0 + t))
    m(t, u, v) = ∫(u * v)dΩ
    c(t, u, v) = ∫(0.1 * u * v)dΩ
    a(t, u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
    l(t, v) = ∫(sin(t) * v)dΩ
    op = FO.AlgebraicLinearTFEOperator(
      TransientLinearFEOperator((a, c, m), l, U, V; constant_forms = (true, true, true)))
    n = num_free_dofs(V)
    r_alg, r_std, cache = _residual_pair(op, 0.3, (randn(n), randn(n), randn(n)))
    @test !FO._algebraic_residual_applies(cache)
    @test r_alg ≈ r_std
  end
end
