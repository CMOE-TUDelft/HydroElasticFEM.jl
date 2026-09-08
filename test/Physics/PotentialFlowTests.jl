using Test
using Gridap
using WaveSpec
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.ParameterHandler as PH

function _single_frequency_state(; H=0.2, T=5.0, h=10.0, θ=0.0)
  spec = WaveSpec.ContinuousSpectrums.RegularWave(H, T)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(θ)
  ω = [2π / T]
  k = [WaveSpec.AiryWaves.solve_wavenumber(ω[1], h)]
  θ_vec = [θ]
  WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, ω, k, θ_vec, h, 1)
end

function _multi_frequency_state()
  single = _single_frequency_state()
  ω = [single.ω[1], 1.25 * single.ω[1]]
  k = [single.k[1], WaveSpec.AiryWaves.solve_wavenumber(ω[2], single.h)]
  WaveSpec.AiryWaves.AiryState(single.spectrum, single.spread, length(ω), single.nθ, ω, k, single.θ, single.h, single.seed)
end

function _potential_flow_problem(pf; ω=2π / 5.0)
  tank = G.TankDomain(L=20.0, H=10.0, nx=12, ny=4)
  config = SM.FreqDomainConfig(ω=ω)
  SM.build_problem(tank, P.PhysicsParameters[pf], config)
end

function _potential_flow_time_problem(pf)
  tank = G.TankDomain(L=20.0, H=10.0, nx=12, ny=4)
  config = SM.TimeDomainConfig(t₀=0.0, tf=0.1)
  SM.build_problem(tank, P.PhysicsParameters[pf], config)
end

@testset "PotentialFlow struct" begin
  pf = P.PotentialFlow()
  @test pf isa P.PhysicsParameters
  @test pf.ρw == 1025.0
  @test pf.g  == 9.81
  @test isempty(pf.boundary_conditions)

  # Custom values
  pf2 = P.PotentialFlow(ρw=1000.0, g=9.80)
  @test pf2.ρw == 1000.0
  @test pf2.g  == 9.80

  # variable_symbol
  @test P.variable_symbol(pf) == :ϕ

  # trait queries
  @test P.has_mass_form(pf)    == false
  @test P.has_damping_form(pf) == false
  @test P.has_stiffness_form(pf) == true
  @test P.has_rhs_form(pf)      == true

  # print_parameters should not throw
  @test_nowarn P.print_parameters(pf)
end

@testset "PotentialFlow boundary conditions" begin
  state = _single_frequency_state()
  ω = state.ω[1]

  pf_rad = P.PotentialFlow(
    sea_state=state,
    boundary_conditions=[P.RadiationBC()],
    fe=PH.FESpaceConfig(order=1, vector_type=Vector{ComplexF64}),
  )
  problem_rad = _potential_flow_problem(pf_rad; ω=ω)
  @test_nowarn SM.simulate(problem_rad)

  pf_multi = P.PotentialFlow(
    sea_state=_multi_frequency_state(),
    boundary_conditions=[P.RadiationBC()],
    fe=PH.FESpaceConfig(order=1, vector_type=Vector{ComplexF64}),
  )
  @test_throws ErrorException _potential_flow_problem(pf_multi; ω=ω)

  @test_nowarn _potential_flow_time_problem(pf_rad)
  @test_throws ErrorException _potential_flow_time_problem(pf_multi)

  pf_inlet = P.PotentialFlow(
    sea_state=state,
    boundary_conditions=[
      P.RadiationBC(domain=:dΓin),
      P.RadiationBC(domain=:dΓout),
      P.PrescribedInletPotentialBC(forcing=(x -> 1.0 + 0.0im), quantity=:potential),
    ],
    fe=PH.FESpaceConfig(order=1, vector_type=Vector{ComplexF64}),
  )
  problem_inlet = _potential_flow_problem(pf_inlet; ω=ω)
  result_inlet = SM.simulate(problem_inlet)
  ϕₕ = result_inlet.solution[1]
  dΩ = SM.get_integration_domains(problem_inlet)[:dΩ]
  ϕ_l2 = sqrt(abs(sum(∫(ϕₕ * conj(ϕₕ))dΩ)))
  @test ϕ_l2 > 0.0
  @test isfinite(ϕ_l2)

  pf_no_inlet_time = P.PotentialFlow(
    sea_state=state,
    boundary_conditions=[],
    fe=PH.FESpaceConfig(order=1, vector_type=Vector{ComplexF64}),
  )

  pf_inlet_time = P.PotentialFlow(
    sea_state=state,
    boundary_conditions=[
      P.PrescribedInletPotentialBC(forcing=(x -> 1.0 + 0.0im), quantity=:normal_gradient),
    ],
    fe=PH.FESpaceConfig(order=1, vector_type=Vector{ComplexF64}),
  )

  problem_no_inlet_time = _potential_flow_time_problem(pf_no_inlet_time)
  problem_inlet_time = _potential_flow_time_problem(pf_inlet_time)

  tconfig = PH.TimeConfig(Δt=0.05, tf=0.1, u0=[0.0])
  result_no_inlet_time = SM.simulate(problem_no_inlet_time, tconfig)
  result_inlet_time = SM.simulate(problem_inlet_time, tconfig)

  sol_no = collect(result_no_inlet_time.solution)
  sol_in = collect(result_inlet_time.solution)

  _, uh_no = sol_no[end]
  _, uh_in = sol_in[end]

  ϕ_no = uh_no[1]
  ϕ_in = uh_in[1]

  dΩ = SM.get_integration_domains(problem_inlet_time)[:dΩ]
  V = SM.get_trial_fe_space(problem_inlet_time)
  v_no = get_free_dof_values(ϕ_no)
  v_in = get_free_dof_values(ϕ_in)
  ϕ_diff_mf = FEFunction(V, v_in - v_no)
  ϕ_diff = ϕ_diff_mf[1]
  diff_l2 = sqrt(abs(sum(∫(ϕ_diff * conj(ϕ_diff))dΩ)))

  @test isfinite(diff_l2)
  @test diff_l2 > 0.0
end
