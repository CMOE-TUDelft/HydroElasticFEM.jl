using Gridap
using Gridap.ODEs
using Gridap.CellData
using WaveSpec

import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as FES
import HydroElasticFEM.Simulation.FESpaceAssembly as FA
import HydroElasticFEM.Geometry as G

include("FEOperatorsTests.jl")
include("FESpaceAssemblyTests.jl")

function _single_frequency_state(; H=0.2, T=5.0, h=10.0, θ=0.0)
  spec = WaveSpec.ContinuousSpectrums.RegularWave(H, T)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(θ)
  ω = [2π / T]
  k = [WaveSpec.AiryWaves.solve_wavenumber(ω[1], h)]
  θ_vec = [θ]
  WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, ω, k, θ_vec, h, 1)
end

# =========================================================================
# Shared mesh setup: 50m × 10m tank, membrane from x=15 to x=35
# Uses CartesianGeometry functions from the Geometry module.
# =========================================================================

@testset "Simulation" begin

  order = 1

  mem_domain = G.StructureDomain1D(L=20.0, x₀=[15.0, 10.0])
  tank = G.TankDomain2D(L=50.0, H=10.0, nx=20, ny=4,
      structure_domains=[mem_domain])

  model = G.build_model(tank)
  trian = G.build_triangulations(tank, model)
  dom   = G.get_integration_domains(trian; degree=2*order)

  # Extra measure on the inlet for the RHS forcing
  dΓin = Measure(trian[:Γin], 2*order)
  dΩ   = Measure(trian[:Ω], 2*order)

  fluid = P.PotentialFlow(ρw=1025.0, g=9.81,
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}), space_domain_symbol=:Ω)
  fsurf = P.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5,
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}), space_domain_symbol=:Γκ)
  mem   = P.Membrane2D(L=20.0, mᵨ=0.9, Tᵨ=98.1,
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}), space_domain_symbol=:Γη)

  # =========================================================================
  # detect_couplings
  # =========================================================================

  @testset "detect_couplings" begin
    entities = [fluid, fsurf, mem]
    pairs = SM.detect_couplings(entities)

    # PotentialFlow ↔ FreeSurface: has mass and damping coupling
    @test (fluid, fsurf) in pairs

    # PotentialFlow ↔ Membrane2D (PhysicsParameters): has damping coupling
    @test (fluid, mem) in pairs

    # No self-coupling
    @test !any(p -> p[1] === p[2], pairs)

    # FreeSurface ↔ Membrane2D: no coupling defined
    @test !((fsurf, mem) in pairs)
    @test !((mem, fsurf) in pairs)
  end

  # =========================================================================
  # build_fe_operator — frequency domain
  # =========================================================================

  @testset "build_fe_operator — frequency domain" begin
    entities = [fluid, fsurf, mem]
    ω = 2.0

    X, Y, fmap = FA.build_fe_spaces(entities, trian, PH.FreqDomainConfig(ω=ω))

    # With explicit rhs_fn
    rhs_fn(y) = [1.0, 0.0, 0.0]  # corresponds to ϕ, κ, η order in fmap
    op = SM.build_fe_operator(entities, dom, ω, fmap, X, Y;
                              rhs_fn=rhs_fn)
    @test op isa Gridap.FESpaces.AffineFEOperator

    uh = solve(LUSolver(), op)
    ϕₕ, κₕ, ηₕ = uh
    ϕ_l2 = sqrt(abs(sum(∫(ϕₕ * conj(ϕₕ))dΩ)))
    @test ϕ_l2 > 0.0
    @test isfinite(ϕ_l2)

    # With rhs_fn=nothing (zero RHS)
    op0 = SM.build_fe_operator(entities, dom, ω, fmap, X, Y)
    @test op0 isa Gridap.FESpaces.AffineFEOperator

    uh0 = solve(LUSolver(), op0)
    ϕ0, _, _ = uh0
    ϕ0_l2 = sqrt(abs(sum(∫(ϕ0 * conj(ϕ0))dΩ)))
    @test ϕ0_l2 ≈ 0.0 atol=1e-12
  end

  # =========================================================================
  # build_fe_operator — time domain
  # =========================================================================

  @testset "build_fe_operator — time domain" begin

    fluid_real = P.PotentialFlow(ρw=1025.0, g=9.81,
      fe=FES.FESpaceConfig(order=order), space_domain_symbol=:Ω)
    fsurf_real = P.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5,
      fe=FES.FESpaceConfig(order=order), space_domain_symbol=:Γκ)
    mem_real   = P.Membrane2D(L=20.0, mᵨ=0.9, Tᵨ=98.1,
      fe=FES.FESpaceConfig(order=order), space_domain_symbol=:Γη)

    entities = [fluid_real, fsurf_real, mem_real]

    X, Y, fmap = FA.build_fe_spaces(entities, trian, PH.TimeDomainConfig())

    # With explicit rhs_fn
    ω_f = 2.0
    rhs_fn(t, y) = [cos(ω_f * t), 0.0, 0.0]  # corresponds to ϕ, κ, η order in fmap
    op = SM.build_fe_operator(entities, dom, fmap, X, Y;
                              rhs_fn=rhs_fn)
    @test op isa Gridap.ODEs.TransientFEOperator

    # With rhs_fn=nothing (zero RHS)
    op0 = SM.build_fe_operator(entities, dom, fmap, X, Y)
    @test op0 isa Gridap.ODEs.TransientFEOperator
  end

  # =========================================================================
  # Frequency-domain simulate
  # =========================================================================

  @testset "simulate — frequency domain" begin
    ω = 2.0
    config = SM.FreqDomainConfig(ω=ω)

    # Simple RHS: unit forcing on inlet
    rhs_fn(y) = [1.0, 0.0, 0.0]  # corresponds to ϕ, κ, η order in fmap

    entities = [fluid, fsurf, mem]
    problem = SM.build_problem(tank, entities, config; rhs_fn=rhs_fn)
    result = SM.simulate(problem)

    @test result isa SM.SimResult
    @test result.fmap[:ϕ] == 1
    @test result.fmap[:κ] == 2
    @test result.fmap[:η_m] == 3

    # Operator was created
    @test result.op isa Gridap.FESpaces.AffineFEOperator

    # Solution is a valid FE function
    uh = result.solution
    ϕₕ, κₕ, ηₕ = uh

    measures = SM.get_integration_domains(problem)
    dΩ_ = measures[:dΩ]

    # Non-trivial solution (forcing is non-zero)
    ϕ_l2 = sqrt(abs(sum(∫(ϕₕ * conj(ϕₕ))dΩ_)))
    @test ϕ_l2 > 0.0
    @test isfinite(ϕ_l2)
  end

  @testset "simulate — potential flow boundary conditions" begin
    state = _single_frequency_state(h=10.0)
    ω = state.ω[1]
    config = SM.FreqDomainConfig(ω=ω)

    fluid_bc_only = P.PotentialFlow(
      ρw=1025.0,
      g=9.81,
      sea_state=state,
      boundary_conditions=[
        P.PrescribedInletPotentialBC(domain=:dΓin, forcing=(x -> 1.0 + 0.0im), quantity=:traction),
      ],
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}),
      space_domain_symbol=:Ω,
    )

    problem_bc_only = SM.build_problem(tank, P.PhysicsParameters[fluid_bc_only], config)
    result_bc_only = SM.simulate(problem_bc_only)
    ϕ_bc_only = result_bc_only.solution[1]
    dΩ_bc_only = SM.get_integration_domains(problem_bc_only)[:dΩ]
    ϕ_l2_bc_only = sqrt(abs(sum(∫(ϕ_bc_only * conj(ϕ_bc_only))dΩ_bc_only)))
    @test ϕ_l2_bc_only > 0.0
    @test isfinite(ϕ_l2_bc_only)

    fluid_mixed = P.PotentialFlow(
      ρw=1025.0,
      g=9.81,
      sea_state=state,
      boundary_conditions=[
        P.RadiationBC(domain=:dΓin),
        P.RadiationBC(domain=:dΓout),
        P.PrescribedInletPotentialBC(domain=:dΓin, forcing=(x -> 1.0 + 0.0im), quantity=:potential),
      ],
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}),
      space_domain_symbol=:Ω,
    )

    problem_mixed = SM.build_problem(tank, P.PhysicsParameters[fluid_mixed], config)
    result_mixed = SM.simulate(problem_mixed)
    ϕ_mixed = result_mixed.solution[1]
    dΩ_mixed = SM.get_integration_domains(problem_mixed)[:dΩ]
    ϕ_l2_mixed = sqrt(abs(sum(∫(ϕ_mixed * conj(ϕ_mixed))dΩ_mixed)))
    @test ϕ_l2_mixed > 0.0
    @test isfinite(ϕ_l2_mixed)
  end

  @testset "simulate — damping-zone frequency domain" begin
    state = _single_frequency_state(h=10.0)
    ω = state.ω[1]
    config = SM.FreqDomainConfig(ω=ω)

    tank_damp = G.TankDomain2D(
      L=50.0,
      H=10.0,
      nx=20,
      ny=4,
      damping_zones=[
        G.DampingZone1D(L=5.0, x₀=[0.0, 10.0], domain_symbol=:Γ_d_in),
        G.DampingZone1D(L=5.0, x₀=[45.0, 10.0], domain_symbol=:Γ_d_out),
      ],
    )

    fluid_active = P.PotentialFlow(
      ρw=1025.0,
      g=9.81,
      sea_state=state,
      boundary_conditions=[
        P.RadiationBC(domain=:dΓin),
        P.RadiationBC(domain=:dΓout),
        P.DampingZoneBC(
          domain=:dΓd_1,
          μ₁=(x -> 1.5 + 0.0im),
          μ₂=(x -> 0.8 + 0.0im),
          η_in=(x -> 0.1 + 0.0im),
          vz_in=(x -> 0.05 + 0.0im),
        ),
      ],
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}),
      space_domain_symbol=:Ω,
    )

    fsurf_active = P.FreeSurface(
      ρw=1025.0,
      g=9.81,
      βₕ=0.5,
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}),
      space_domain_symbol=:Γκ,
    )

    problem_active = SM.build_problem(tank_damp, P.PhysicsParameters[fluid_active, fsurf_active], config)
    dom_active = SM.get_integration_domains(problem_active)
    @test haskey(dom_active, :dΓfs)
    @test haskey(dom_active, :nΓd_1)
    @test haskey(dom_active, :nΓd_2)

    result_active = SM.simulate(problem_active)
    ϕ_active, κ_active = result_active.solution
    dΩ_active = dom_active[:dΩ]
    ϕ_l2_active = sqrt(abs(sum(∫(ϕ_active * conj(ϕ_active))dΩ_active)))
    κ_l2_active = sqrt(abs(sum(∫(κ_active * conj(κ_active))dom_active[:dΓκ])))
    @test ϕ_l2_active > 0.0
    @test κ_l2_active > 0.0
    @test isfinite(ϕ_l2_active)
    @test isfinite(κ_l2_active)

    fluid_disabled = P.PotentialFlow(
      ρw=1025.0,
      g=9.81,
      sea_state=state,
      boundary_conditions=[
        P.RadiationBC(domain=:dΓin),
        P.RadiationBC(domain=:dΓout),
        P.DampingZoneBC(
          domain=:dΓd_1,
          μ₁=(x -> 1.5 + 0.0im),
          μ₂=(x -> 0.8 + 0.0im),
          η_in=(x -> 0.1 + 0.0im),
          vz_in=(x -> 0.05 + 0.0im),
          enabled=false,
        ),
      ],
      fe=FES.FESpaceConfig(order=order, vector_type=Vector{ComplexF64}),
      space_domain_symbol=:Ω,
    )

    problem_disabled = SM.build_problem(tank_damp, P.PhysicsParameters[fluid_disabled, fsurf_active], config)
    result_disabled = SM.simulate(problem_disabled)
    ϕ_disabled = result_disabled.solution[1]
    dΩ_disabled = SM.get_integration_domains(problem_disabled)[:dΩ]
    ϕ_l2_disabled = sqrt(abs(sum(∫(ϕ_disabled * conj(ϕ_disabled))dΩ_disabled)))
    @test abs(ϕ_l2_active - ϕ_l2_disabled) > 1e-8
  end

  # =========================================================================
  # SimConfig validation
  # =========================================================================

  @testset "SimConfig validation" begin
    # Frequency-domain requires ω
    config_ω = SM.FreqDomainConfig(ω=2.0)
    @test config_ω.ω == 2.0
    @test config_ω.solver === nothing
    @test isa(config_ω, SM.FreqDomainConfig)

    # Time-domain without ω is fine
    config_t = SM.TimeDomainConfig()
    @test config_t.t₀ == 0.0
  end

  # =========================================================================
  # TimeConfig validation
  # =========================================================================

  @testset "TimeConfig validation" begin
    @test_throws AssertionError SM.TimeConfig(Δt=-0.1, tf=1.0)
    @test_throws AssertionError SM.TimeConfig(Δt=0.1, t₀=2.0, tf=1.0)

    tc = SM.TimeConfig(Δt=0.01, tf=0.1, u0=[0.0, 0.0, 0.0])
    @test tc.Δt == 0.01
    @test tc.t₀ == 0.0
    @test tc.ρ∞ == 1.0
  end

  # =========================================================================
  # Time-domain simulate
  # =========================================================================

  @testset "simulate — time domain" begin
    config = SM.TimeDomainConfig(t₀=0.0, tf=0.3)
    tconfig = SM.TimeConfig(
        Δt  = 0.1,
        t₀  = 0.0,
        tf  = 0.3,
        ρ∞  = 1.0,
        u0  = [0.0, 0.0, 0.0],
        u0t = [0.0, 0.0, 0.0],
        u0tt = [0.0, 0.0, 0.0],
    )

    # Time-dependent rhs: simple harmonic forcing on inlet
    ω_f = 2.0
    rhs_fn(t, y) = [1.0 * t, 0.0, 0.0]  # corresponds to ϕ, κ, η order in fmap

    fluid_real = P.PotentialFlow(ρw=1025.0, g=9.81,
      fe=FES.FESpaceConfig(order=order), space_domain_symbol=:Ω)
    fsurf_real = P.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5,
      fe=FES.FESpaceConfig(order=order), space_domain_symbol=:Γκ)
    mem_real   = P.Membrane2D(L=20.0, mᵨ=0.9, Tᵨ=98.1,
      fe=FES.FESpaceConfig(order=order), space_domain_symbol=:Γη)

    entities = [fluid_real, fsurf_real, mem_real]
    problem = SM.build_problem(tank, entities, config; rhs_fn=rhs_fn)

    result = SM.simulate(problem, tconfig)

    @test result isa SM.SimResult
    @test result.fmap[:ϕ] == 1
    @test result.fmap[:κ] == 2
    @test result.fmap[:η_m] == 3

    # Iterate over ODE solution
    step_count = 0
    for (t, uh) in result.solution
      step_count += 1
      ϕₕ, κₕ, ηₕ = uh
      # Check solution is finite
      @test all(isfinite, get_free_dof_values(ϕₕ))
    end

    expected_steps = round(Int, (tconfig.tf - tconfig.t₀) / tconfig.Δt)
    @test step_count == expected_steps
  end

  # =========================================================================
  # Custom degree dictionary for integration domains
  # =========================================================================
  @testset "integration domains with custom degree dict" begin
    degree_dict = Dict(:dΩ => 3, :dΓκ => 2, :dΓη => 5, :dΓin => 1, :dΓout => 1, :dΓbot => 1)
    dom_custom = G.get_integration_domains(trian; degree=degree_dict)

    # Run a simple frequency-domain simulation with the custom dom
    ω = 2.0
    config = SM.FreqDomainConfig(ω=ω)
    rhs_fn(y) = [1.0, 0.0, 0.0]  # corresponds to ϕ, κ, η order in fmap
    entities = [fluid, fsurf, mem]
    problem = SM.build_problem(tank, entities, config; rhs_fn=rhs_fn)
    result = SM.simulate(problem)
    measures = SM.get_integration_domains(problem)
    dΩ_ = measures[:dΩ]
    @test result isa SM.SimResult
    @test result.op isa Gridap.FESpaces.AffineFEOperator
    uh = result.solution
    ϕₕ, κₕ, ηₕ = uh
    ϕ_l2 = sqrt(abs(sum(∫(ϕₕ * conj(ϕₕ))dΩ_)))
    @test ϕ_l2 > 0.0
    @test isfinite(ϕ_l2)
  end

end
