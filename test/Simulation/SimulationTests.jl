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
  tri   = G.build_triangulations(tank, model)
  dom   = G.get_integration_domains(tri; degree=2*order)

  # Extra measure on the inlet for the RHS forcing
  dΓin = Measure(tri.Γin, 2*order)
  dΩ   = Measure(tri.Ω, 2*order)

  fluid = P.PotentialFlow(ρw=1025.0, g=9.81,
      fe=FES.FESpaceConfig(order=order))
  fsurf = P.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5,
      fe=FES.FESpaceConfig(order=order))
  mem   = P.Membrane2D(L=20.0, mᵨ=0.9, Tᵨ=98.1,
      fe=FES.FESpaceConfig(order=order))

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
    coupling_pairs = SM.detect_couplings(entities)
    ω = 2.0

    X, Y, fmap = FA.build_fe_spaces(
        fluid => tri.Ω, fsurf => tri.Γκ, mem => tri.Γη)

    # With explicit rhs_fn
    rhs_fn(y) = ∫(y[:ϕ] * 1.0)dΓin
    op = SM.build_fe_operator(entities, coupling_pairs, dom, ω, fmap, X, Y;
                              rhs_fn=rhs_fn)
    @test op isa Gridap.FESpaces.AffineFEOperator

    uh = solve(LUSolver(), op)
    ϕₕ, κₕ, ηₕ = uh
    ϕ_l2 = sqrt(abs(sum(∫(ϕₕ * conj(ϕₕ))dΩ)))
    @test ϕ_l2 > 0.0
    @test isfinite(ϕ_l2)

    # With rhs_fn=nothing (zero RHS)
    op0 = SM.build_fe_operator(entities, coupling_pairs, dom, ω, fmap, X, Y)
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
    entities = [fluid, fsurf, mem]
    coupling_pairs = SM.detect_couplings(entities)

    X, Y, fmap = FA.build_fe_spaces(
        fluid => tri.Ω, fsurf => tri.Γκ, mem => tri.Γη; transient=true)

    # With explicit rhs_fn
    ω_f = 2.0
    rhs_fn(t, y) = ∫(y[:ϕ] * cos(ω_f * t))dΓin
    op = SM.build_fe_operator(entities, coupling_pairs, dom, fmap, X, Y;
                              rhs_fn=rhs_fn)
    @test op isa Gridap.ODEs.TransientFEOperator

    # With rhs_fn=nothing (zero RHS)
    op0 = SM.build_fe_operator(entities, coupling_pairs, dom, fmap, X, Y)
    @test op0 isa Gridap.ODEs.TransientFEOperator
  end

  # =========================================================================
  # Frequency-domain simulate
  # =========================================================================

  @testset "simulate — frequency domain" begin
    ω = 2.0
    config = SM.SimConfig(domain=:frequency, ω=ω)

    # Simple RHS: unit forcing on inlet
    rhs_fn(y) = ∫(y[:ϕ] * 1.0)dΓin

    result = SM.simulate(config,
        fluid => tri.Ω,
        fsurf => tri.Γκ,
        mem   => tri.Γη;
        dom=dom,
        rhs_fn=rhs_fn)

    @test result isa SM.SimResult
    @test result.fmap[:ϕ] == 1
    @test result.fmap[:κ] == 2
    @test result.fmap[:η_m] == 3

    # Operator was created
    @test result.op isa Gridap.FESpaces.AffineFEOperator

    # Solution is a valid FE function
    uh = result.solution
    ϕₕ, κₕ, ηₕ = uh

    # Non-trivial solution (forcing is non-zero)
    ϕ_l2 = sqrt(abs(sum(∫(ϕₕ * conj(ϕₕ))dΩ)))
    @test ϕ_l2 > 0.0
    @test isfinite(ϕ_l2)
  end

  # =========================================================================
  # SimConfig validation
  # =========================================================================

  @testset "SimConfig validation" begin
    # Frequency-domain requires ω
    @test_throws AssertionError SM.SimConfig(domain=:frequency)

    # Invalid domain
    @test_throws AssertionError SM.SimConfig(domain=:invalid)

    # Time-domain without ω is fine
    config_t = SM.SimConfig(domain=:time)
    @test config_t.domain == :time
    @test config_t.ω === nothing
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
    config = SM.SimConfig(domain=:time)
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
    rhs_fn(t, y) = ∫(y[:ϕ] * cos(ω_f * t))dΓin

    result = SM.simulate(config, tconfig,
        fluid => tri.Ω,
        fsurf => tri.Γκ,
        mem   => tri.Γη;
        dom=dom,
        rhs_fn=rhs_fn)

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

end
