using Test
using Gridap

import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Geometry as G

@testset "Time-domain membrane+damping zone package API" begin
  order = 1
  H0 = 10.0
  Lm = 10.0
  Ld = 10.0
  LΩ = 2 * Ld + 3 * Lm

  βₕ = 0.5
  γₜ = 0.5
  βₜ = 0.25
  Δt = 0.1
  αₕ = γₜ / (βₜ * Δt) / 9.81 * (1.0 - βₕ) / βₕ

  ω = 1.0
  η₀ = 0.05
  ph0 = π / 2

  η_in(x, t) = η₀ * cos(0.2 * x[1] - ω * t + ph0)
  vz_in(x, t) = η₀ * ω * sin(0.2 * x[1] - ω * t + ph0)
  inlet_v(x, t) = 0.1 * cos(0.2 * x[1] - ω * t + ph0)

  tank = G.TankDomain2D(
    L=LΩ,
    H=H0,
    nx=20,
    ny=4,
    structure_domains=[
      G.StructureDomain1D(L=Lm, x₀=[Ld + Lm / 2, H0], domain_symbol=:Γm),
    ],
    damping_zones=[
      G.DampingZone1D(L=Ld, x₀=[0.0, H0], domain_symbol=:Γd_in),
      G.DampingZone1D(L=Ld, x₀=[LΩ - Ld, H0], domain_symbol=:Γd_out),
    ],
  )

  fluid = P.PotentialFlow(
    ρw=1025.0,
    g=9.81,
    boundary_conditions=[
      P.PrescribedInletPotentialBC(domain=:dΓin, forcing=(t -> (x -> inlet_v(x, t))), quantity=:traction),
      P.DampingZoneBC(
        domain=:dΓd_1,
        μ₁=(x -> 2.0),
        μ₂=(x -> 1.0),
        η_in=(t -> (x -> η_in(x, t))),
        vz_in=(t -> (x -> vz_in(x, t))),
      ),
      P.DampingZoneBC(
        domain=:dΓd_2,
        μ₁=(x -> 2.0),
        μ₂=(x -> 1.0),
        η_in=(x -> 0.0),
        vz_in=(x -> 0.0),
      ),
    ],
    fe=PH.FESpaceConfig(order=order),
    space_domain_symbol=:Ω,
  )

  fsurf = P.FreeSurface(
    ρw=1025.0,
    g=9.81,
    βₕ=βₕ,
    fe=PH.FESpaceConfig(order=order),
    space_domain_symbol=:Γκ,
  )

  membrane = P.Membrane2D(
    L=Lm,
    mᵨ=0.9,
    Tᵨ=98.1,
    τ=0.0,
    g=9.81,
    fe=PH.FESpaceConfig(order=order),
    space_domain_symbol=:Γη,
  )

  entities = P.PhysicsParameters[fluid, fsurf, membrane]
  config = SM.TimeDomainConfig(t₀=0.0, tf=0.3)
  tconfig = SM.TimeConfig(
    Δt=Δt,
    t₀=0.0,
    tf=0.3,
    ρ∞=1.0,
    αₕ=αₕ,
    u0=[0.0, 0.0, 0.0],
    u0t=[0.0, 0.0, 0.0],
    u0tt=[0.0, 0.0, 0.0],
  )

  @test_throws ErrorException SM.build_problem(tank, entities, config)

  problem = SM.build_problem(tank, entities, config; tconfig=tconfig)
  dom = SM.get_integration_domains(problem)
  @test haskey(dom, :dΓfs)
  @test haskey(dom, :nΓd_1)
  @test haskey(dom, :nΓd_2)
  @test dom[:αₕ] == αₕ

  result = SM.simulate(problem, tconfig)

  step_count = 0
  ϕ_norms = Float64[]
  η_norms = Float64[]

  for (_, uh) in result.solution
    step_count += 1
    ϕₕ, κₕ, ηₕ = uh
    push!(ϕ_norms, maximum(abs.(get_free_dof_values(ϕₕ))))
    push!(η_norms, maximum(abs.(get_free_dof_values(ηₕ))))
    @test all(isfinite, get_free_dof_values(ϕₕ))
    @test all(isfinite, get_free_dof_values(κₕ))
    @test all(isfinite, get_free_dof_values(ηₕ))
  end

  @test step_count == round(Int, (tconfig.tf - tconfig.t₀) / tconfig.Δt)
  @test maximum(ϕ_norms) > 0.0
  @test maximum(η_norms) > 0.0
end
