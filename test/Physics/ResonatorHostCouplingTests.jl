using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO

# Two floating membranes with separate fields; one resonator on each.  Every
# resonator must couple only to the structure on its host domain.
@testset "Resonator ↔ structure coupling follows the host domain" begin
  tank = G.TankDomain(L = 20.0, H = 2.0, nx = 40, ny = 2,
    structure_domains = [G.StructureDomain(L = 4.0, x₀ = [4.0, 2.0], domain_symbol = :Γ_a),
                         G.StructureDomain(L = 4.0, x₀ = [12.0, 2.0], domain_symbol = :Γ_b)],
    resonator_domains = [G.ResonatorDomain(location = [5.1, 2.0], trian_symbol = :Γ_a, delta_symbol = :δ1),
                         G.ResonatorDomain(location = [13.3, 2.0], trian_symbol = :Γ_b, delta_symbol = :δ2)])
  fe = PH.FESpaceConfig(vector_type = Vector{Float64})
  mem_a = P.Membrane(L = 4.0, mᵨ = 0.9, Tᵨ = 98.1, fe = fe, symbol = :η_a, space_domain_symbol = :Γ_a)
  mem_b = P.Membrane(L = 4.0, mᵨ = 0.9, Tᵨ = 98.1, fe = fe, symbol = :η_b, space_domain_symbol = :Γ_b)
  ra = P.ResonatorArray([
    P.ResonatorSingle(M = 100.0, K = 500.0, C = 5.0, symbol = :q1, host_domain_symbol = :Γ_a,
                      delta_domain_symbol = :δ1, fe = fe),
    P.ResonatorSingle(M = 150.0, K = 700.0, C = 8.0, symbol = :q2, host_domain_symbol = :Γ_b,
                      delta_domain_symbol = :δ2, fe = fe)])

  @test [r.symbol for r in P.attached_resonators(ra, mem_a)] == [:q1]
  @test [r.symbol for r in P.attached_resonators(ra, mem_b)] == [:q2]
  @test P.has_stiffness_form(ra, mem_a) && P.has_damping_form(ra, mem_b)
  mem_η = P.Membrane(L = 4.0, mᵨ = 0.9, Tᵨ = 98.1, fe = fe)          # on :Γη
  @test !P.has_stiffness_form(ra, mem_η)

  inlet = P.PrescribedInletPotentialBC(domain = :dΓin, quantity = :traction,
                                       forcing = (t -> (x -> 0.1 * sin(2t))))
  physics = P.PhysicsParameters[P.PotentialFlow(fe = fe, boundary_conditions = [inlet]),
                                P.FreeSurface(fe = fe), mem_a, mem_b, ra]
  tr = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(tr; degree = 2)
  X, Y, fmap = SM.build_fe_spaces(physics, tr, PH.TimeDomainConfig())
  U = X(0.0)
  Kc(s) = assemble_matrix((x, y) -> P.stiffness(ra, s, dom, FO.FieldMap(x, fmap), FO.FieldMap(y, fmap)), U, Y)

  # Block structure: (ra, mem_a) touches only q1 and η_a
  offs = cumsum([0; [num_free_dofs(V) for V in Y.spaces]])
  blk(A, i, j) = A[offs[i]+1:offs[i+1], offs[j]+1:offs[j+1]]
  Ka = Kc(mem_a)
  for (i, j) in ((fmap[:q2], fmap[:q2]), (fmap[:q2], fmap[:η_a]), (fmap[:η_b], fmap[:q1]), (fmap[:η_b], fmap[:η_b]))
    @test iszero(blk(Ka, i, j))
  end
  @test !iszero(blk(Ka, fmap[:q1], fmap[:η_a]))
  @test !iszero(blk(Ka, fmap[:η_a], fmap[:q1]))
  Kb = Kc(mem_b)
  @test iszero(blk(Kb, fmap[:q1], fmap[:η_a]))
  @test !iszero(blk(Kb, fmap[:q2], fmap[:η_b]))

  # End to end: both resonators respond to the incoming forcing
  cfg = PH.TimeDomainConfig(t₀ = 0.0, tf = 0.5)
  tc = PH.TimeConfig(Δt = 0.1, tf = 0.5, u0 = zeros(6), u0t = zeros(6), u0tt = zeros(6))
  prob = SM.build_problem(tank, physics, cfg; tconfig = tc)
  pairs = SM.detect_couplings(physics)
  @test count(p -> p[1] === ra, pairs) == 2
  q1max = q2max = 0.0
  for (_, uh) in SM.simulate(prob, tc).solution
    q1max = max(q1max, maximum(abs, get_free_dof_values(uh[fmap[:q1]])))
    q2max = max(q2max, maximum(abs, get_free_dof_values(uh[fmap[:q2]])))
  end
  @test q1max > 0 && q2max > 0
end
