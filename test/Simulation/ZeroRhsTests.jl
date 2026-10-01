using Test
using Gridap
import HydroElasticFEM
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.AssemblyContexts as AC

@testset "Zero forcing assembles no integrals" begin
  tank = G.TankDomain(L = 10.0, H = 2.0, nx = 10, ny = 2,
    structure_domains = [G.StructureDomain(L = 4.0, x₀ = [5.0, 2.0])])
  trians = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(trians; degree = 2)
  fe = PH.FESpaceConfig(vector_type = Vector{Float64})
  pf = P.PotentialFlow(fe = fe)
  mem = P.Membrane(L = 4.0, mᵨ = 0.9, Tᵨ = 98.1, fe = fe)
  sη = P.variable_symbol(mem)
  V = TestFESpace(trians[:Ω], ReferenceFE(lagrangian, Float64, 1); conformity = :H1)
  Vη = TestFESpace(trians[:Γη], ReferenceFE(lagrangian, Float64, 1); conformity = :H1)

  # Exact numeric zero → no contribution; any other forcing is integrated
  @test P._forcing_contribution(0.0, nothing, nothing) === nothing
  @test P._forcing_contribution(0.0 + 0.0im, nothing, nothing) === nothing
  vη = get_fe_basis(Vη)
  @test P.rhs(mem, dom, Dict(sη => 0.0), Dict(sη => vη)) === nothing
  @test P.rhs(pf, dom, Dict(:ϕ => 0.0), Dict(:ϕ => get_fe_basis(V))) === nothing
  f = assemble_vector(v -> P.rhs(mem, dom, Dict(sη => 2.0), Dict(sη => v)), Vη)
  @test sum(f) ≈ 2.0 * 4.0
  fx = assemble_vector(v -> P.rhs(mem, dom, Dict(sη => (x -> x[1])), Dict(sη => v)), Vη)
  @test sum(fx) ≈ (9.0^2 - 5.0^2) / 2   # membrane spans x ∈ [5, 9]

  # Whole-problem RHS without any forcing: empty form, zero load vector
  fs = P.FreeSurface(fe = fe)
  physics = P.PhysicsParameters[pf, fs, mem]
  X, Y, fmap = SM.build_fe_spaces(physics, trians, PH.TimeDomainConfig())
  ctx = AC.TimeAssemblyContext(dom, 0.0, nothing)
  pairs = SM.detect_couplings(physics, ctx)
  l(y) = FO._assemble_rhs_total(physics, pairs, ctx, fmap, zeros(length(fmap)), y)
  @test l(get_fe_basis(Y)) isa Gridap.CellData.DomainContribution
  @test iszero(assemble_vector(l, Y))

  # Time-domain simulate still runs from a non-trivial initial state
  cfg = PH.TimeDomainConfig(t₀ = 0.0, tf = 0.3)
  tc = PH.TimeConfig(Δt = 0.1, tf = 0.3, u0 = [0.0, x -> 0.01 * cos(π * x[1] / 10), 0.0],
                     u0t = [0.0, 0.0, 0.0], u0tt = [0.0, 0.0, 0.0])
  prob = SM.build_problem(tank, physics, cfg; tconfig = tc)
  n = 0
  for (_, uh) in SM.simulate(prob, tc).solution
    n += 1
    @test all(isfinite, get_free_dof_values(uh))
  end
  @test n == 3
end
