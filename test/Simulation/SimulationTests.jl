using Gridap
using Gridap.ODEs
using Gridap.CellData

import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as FES
import HydroElasticFEM.Simulation.FESpaceAssembly as FA
import HydroElasticFEM.Geometry as G

include("FEOperatorsTests.jl")
include("FESpaceAssemblyTests.jl")

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
