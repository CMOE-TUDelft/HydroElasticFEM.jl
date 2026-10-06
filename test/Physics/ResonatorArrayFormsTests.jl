using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO

# ResonatorArray self and resonator↔structure forms must equal the sum of the
# point (Dirac) terms only: no host-surface integral is involved.
@testset "ResonatorArray forms are pure point terms" begin
  tank = G.TankDomain(L = 10.0, H = 2.0, nx = 20, ny = 2,
    structure_domains = [G.StructureDomain(L = 4.0, x₀ = [3.0, 2.0])],
    resonator_domains = [G.ResonatorDomain(location = [4.1, 2.0], delta_symbol = :δ1),
                         G.ResonatorDomain(location = [5.7, 2.0], delta_symbol = :δ2)])
  trians = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(trians; degree = 2)
  fe = PH.FESpaceConfig(vector_type = Vector{Float64})
  mem = P.Membrane(L = 4.0, mᵨ = 0.9, Tᵨ = 98.1, fe = fe)
  M, K, C, ρw = [100.0, 150.0], [500.0, 700.0], [5.0, 8.0], 1025.0
  ra = P.resonator_array(2, M, K, C; delta_domain_symbols = [:δ1, :δ2], ρw = ρw, fe = fe)
  physics = P.PhysicsParameters[mem, ra]
  X, Y, fmap = SM.build_fe_spaces(physics, trians, PH.TimeDomainConfig())
  U = X(0.0)
  sη = P.variable_symbol(mem)
  qs = [P.variable_symbol(r) for r in ra.resonators]
  δ = [dom[:δ1], dom[:δ2]]
  e1 = VectorValue(1.0)
  wrap(f) = (u, v) -> f(FO.FieldMap(u, fmap), FO.FieldMap(v, fmap))
  A(f) = assemble_matrix(wrap(f), U, Y)

  # hand-written point terms
  m_ref(x, y) = sum(M[i] * δ[i](x[qs[i]] ⋅ y[qs[i]]) for i in 1:2)
  c_ref(x, y) = sum(C[i] * δ[i](x[qs[i]] ⋅ y[qs[i]]) for i in 1:2)
  k_ref(x, y) = sum(K[i] * δ[i](x[qs[i]] ⋅ y[qs[i]]) for i in 1:2)
  kc_ref(x, y) = sum((-K[i] / ρw) * δ[i](y[sη] * ((x[qs[i]] ⋅ e1) - x[sη])) -
                     K[i] * δ[i]((y[qs[i]] ⋅ e1) * x[sη]) for i in 1:2)
  cc_ref(x, y) = sum((-C[i] / ρw) * δ[i](y[sη] * ((x[qs[i]] ⋅ e1) - x[sη])) -
                     C[i] * δ[i]((y[qs[i]] ⋅ e1) * x[sη]) for i in 1:2)

  for (name, pkg, ref) in (
      ("mass",               (x, y) -> P.mass(ra, dom, x, y),           m_ref),
      ("damping",            (x, y) -> P.damping(ra, dom, x, y),        c_ref),
      ("stiffness",          (x, y) -> P.stiffness(ra, dom, x, y),      k_ref),
      ("coupling stiffness", (x, y) -> P.stiffness(ra, mem, dom, x, y), kc_ref),
      ("coupling damping",   (x, y) -> P.damping(ra, mem, dom, x, y),   cc_ref))
    @testset "$name" begin
      Ap, Ar = A(pkg), A(ref)
      @test norm(Ar) > 0
      @test norm(Ap - Ar) <= 1e-13 * norm(Ar)
    end
  end
end

# End to end: package resonators on a floating membrane (host :Γη) assemble
# and respond to an inlet forcing through build_problem + simulate.
@testset "ResonatorArray on a membrane — time-domain simulate" begin
  tank = G.TankDomain(L = 10.0, H = 2.0, nx = 20, ny = 2,
    structure_domains = [G.StructureDomain(L = 4.0, x₀ = [3.0, 2.0])],
    resonator_domains = [G.ResonatorDomain(location = [4.1, 2.0], delta_symbol = :δ1),
                         G.ResonatorDomain(location = [5.7, 2.0], delta_symbol = :δ2)])
  fe = PH.FESpaceConfig(vector_type = Vector{Float64})
  inlet = P.PrescribedInletPotentialBC(domain = :dΓin, quantity = :traction,
                                       forcing = (t -> (x -> 0.1 * sin(2t))))
  physics = P.PhysicsParameters[P.PotentialFlow(fe = fe, boundary_conditions = [inlet]),
                                P.FreeSurface(fe = fe),
                                P.Membrane(L = 4.0, mᵨ = 0.9, Tᵨ = 98.1, fe = fe),
                                P.resonator_array(2, 100.0, 500.0, 5.0;
                                                  delta_domain_symbols = [:δ1, :δ2], fe = fe)]
  cfg = PH.TimeDomainConfig(t₀ = 0.0, tf = 0.5)
  tc = PH.TimeConfig(Δt = 0.1, tf = 0.5, u0 = zeros(5), u0t = zeros(5), u0tt = zeros(5))
  prob = SM.build_problem(tank, physics, cfg; tconfig = tc)
  fmap = SM.get_field_map(prob)
  qmax = 0.0
  for (_, uh) in SM.simulate(prob, tc).solution
    @test all(isfinite, get_free_dof_values(uh))
    qmax = max(qmax, maximum(abs, get_free_dof_values(uh[fmap[:q1]])))
  end
  @test qmax > 0
end

# The default zero forcing must stay type-consistent with the complex resonator
# matrix terms in the frequency domain.
@testset "ResonatorArray on a membrane — frequency-domain default forcing" begin
  tank = G.TankDomain(L = 10.0, H = 2.0, nx = 20, ny = 2,
    structure_domains = [G.StructureDomain(L = 4.0, x₀ = [3.0, 2.0])],
    resonator_domains = [G.ResonatorDomain(location = [4.1, 2.0], delta_symbol = :δ1),
                         G.ResonatorDomain(location = [5.7, 2.0], delta_symbol = :δ2)])
  fe = PH.FESpaceConfig(vector_type = Vector{ComplexF64})
  physics = P.PhysicsParameters[P.PotentialFlow(fe = fe),
                                P.FreeSurface(fe = fe),
                                P.Membrane(L = 4.0, mᵨ = 0.9, Tᵨ = 98.1, fe = fe),
                                P.resonator_array(2, 100.0, 500.0, 5.0;
                                                  delta_domain_symbols = [:δ1, :δ2], fe = fe)]
  prob = SM.build_problem(tank, physics, SM.FreqDomainConfig(ω = 2.0))
  @test prob !== nothing
end
