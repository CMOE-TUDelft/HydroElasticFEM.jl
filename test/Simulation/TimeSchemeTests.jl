using Test
using Gridap
using Gridap.ODEs
import HydroElasticFEM
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.AssemblyContexts as AC

@testset "TimeConfig — scheme and αₕ" begin

  @testset "time_integration_parameters" begin
    ga = PH.time_integration_parameters(PH.TimeConfig(Δt = 0.1, tf = 1.0))
    @test (ga.γ, ga.β, ga.αf, ga.αm) == (0.5, 0.25, 0.5, 0.5)
    ga8 = PH.time_integration_parameters(PH.TimeConfig(Δt = 0.1, tf = 1.0, ρ∞ = 0.8))
    @test ga8.γ ≈ 0.5 - ga8.αm + ga8.αf
    nm = PH.time_integration_parameters(PH.TimeConfig(Δt = 0.1, tf = 1.0, scheme = :newmark,
                                                      γ = 0.6, β = 0.3025))
    @test (nm.γ, nm.β, nm.αf, nm.αm) == (0.6, 0.3025, 0.0, 0.0)
    @test_throws AssertionError PH.TimeConfig(Δt = 0.1, tf = 1.0, scheme = :euler)
    @test_throws AssertionError PH.TimeConfig(Δt = 0.1, tf = 1.0, αₕ = :magic)
  end

  @testset "stabilization_αₕ matches the hand-written formula" begin
    Δt, g, βₕ = 0.05, 9.81, 0.5
    tc = PH.TimeConfig(Δt = Δt, tf = 1.0, scheme = :newmark)
    @test PH.stabilization_αₕ(tc, g, βₕ) ≈ 0.5 / (0.25 * Δt) / g * (1 - βₕ) / βₕ
  end

  # 2D tank with one damping zone: αₕ is required there
  tank = G.TankDomain(L = 20.0, H = 5.0, nx = 10, ny = 2,
    damping_zones = [G.DampingZone(L = 5.0, x₀ = [15.0, 5.0], domain_symbol = :Γd)])
  trians = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(trians; degree = 2)
  fe = PH.FESpaceConfig(vector_type = Vector{Float64})
  dz = P.DampingZoneBC(domain = :dΓd_1, μ₁ = (x -> 1.0), μ₂ = (x -> 1.0),
                       η_in = (x -> 0.0), vz_in = (x -> 0.0))
  fs = P.FreeSurface(g = 9.81, βₕ = 0.4, fe = fe)
  pf_dz = P.PotentialFlow(g = 9.81, fe = fe, boundary_conditions = [dz])
  pf = P.PotentialFlow(g = 9.81, fe = fe)
  config = PH.TimeDomainConfig(t₀ = 0.0, tf = 1.0)

  @testset "automatic αₕ in build_time_context" begin
    tc = PH.TimeConfig(Δt = 0.1, tf = 1.0, scheme = :newmark)
    expected = PH.stabilization_αₕ(tc, 9.81, 0.4)

    ctx = SM.build_time_context(dom, P.PhysicsParameters[pf_dz, fs], config, tc)
    @test AC.stabilization_parameter(ctx) ≈ expected

    # no damping zone and αₕ = nothing: free surface stays unstabilised
    ctx = SM.build_time_context(dom, P.PhysicsParameters[pf, fs], config, tc)
    @test !AC.has_stabilization(ctx)

    # :auto opts in without damping zones
    tc_auto = PH.TimeConfig(Δt = 0.1, tf = 1.0, scheme = :newmark, αₕ = :auto)
    ctx = SM.build_time_context(dom, P.PhysicsParameters[pf, fs], config, tc_auto)
    @test AC.stabilization_parameter(ctx) ≈ expected

    # an explicit value wins
    tc_val = PH.TimeConfig(Δt = 0.1, tf = 1.0, αₕ = 1.25)
    ctx = SM.build_time_context(dom, P.PhysicsParameters[pf_dz, fs], config, tc_val)
    @test AC.stabilization_parameter(ctx) == 1.25

    # automatic αₕ needs g and βₕ from a FreeSurface
    @test_throws ErrorException SM.build_time_context(dom, P.PhysicsParameters[pf_dz], config, tc)
  end

  @testset "simulate with Newmark from rest matches a hand-built Gridap solve" begin
    Δt, tf = 0.1, 0.5
    inlet = P.PrescribedInletPotentialBC(domain = :dΓin, quantity = :traction,
                                         forcing = (t -> (x -> 0.1 * sin(2t))))
    pf_in = P.PotentialFlow(g = 9.81, fe = fe, boundary_conditions = [inlet, dz])
    physics = P.PhysicsParameters[pf_in, fs]
    tc = PH.TimeConfig(Δt = Δt, tf = tf, scheme = :newmark, γ = 0.5, β = 0.25)
    cfg = PH.TimeDomainConfig(t₀ = 0.0, tf = tf)
    prob = SM.build_problem(tank, physics, cfg; tconfig = tc)
    res = SM.simulate(prob, tc)

    op = SM.get_fe_operator(prob)
    X0 = SM.get_test_fe_space(prob)(0.0)
    ref = solve(Newmark(LUSolver(), Δt, 0.5, 0.25), op, 0.0, tf,
                (zero(X0), zero(X0), zero(X0)))

    nsteps = 0
    for ((t1, u1), (t2, u2)) in zip(res.solution, ref)
      nsteps += 1
      @test t1 ≈ t2
      @test get_free_dof_values(u1) ≈ get_free_dof_values(u2)
    end
    @test nsteps == round(Int, tf / Δt)
    @test maximum(abs, get_free_dof_values(last(collect(res.solution))[2])) > 0
  end
end
