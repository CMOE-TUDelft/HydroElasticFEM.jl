using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO

# Two plates a (x ∈ [2,4]) and b (x ∈ [4,6]), y ∈ [1,3], connected along x = 4,
# compared with a single continuous plate x ∈ [2,6] on the same mesh.
@testset "PlateConnection" begin
  order = 2
  fe = PH.FESpaceConfig(order = order, vector_type = Vector{Float64}, γ = 6.0)
  E, ν, hb = 1.0e8, 0.3, 0.2
  plate(sym, dsym) = P.KirchhoffLovePlate(E = E, ν = ν, hb = hb, fe = fe, symbol = sym, space_domain_symbol = dsym)
  pa, pb = plate(:η_a, :Γ_a), plate(:η_b, :Γ_b)
  sa = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [2.0, 1.0, 1.0], domain_symbol = :Γ_a)
  sb = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [4.0, 1.0, 1.0], domain_symbol = :Γ_b)
  conn = G.StructureConnection(a = :Γ_a, b = :Γ_b, domain_symbol = :dΛ_ab, normal_symbol = :n_Λ_ab)
  mk_tank(structs; kw...) = G.TankDomain(L = 8.0, W = 4.0, H = 1.0, nx = 16, ny = 16, nz = 1,
                                         structure_domains = structs; kw...)
  tank = mk_tank([sa, sb]; structure_connections = [conn])
  tr = G.build_triangulations(tank, G.build_model(tank))
  deg = Dict(:Γ_a => 4, :Γ_b => 4, :default => 4)
  dom = G.get_integration_domains(tr; degree = deg)
  lag = ReferenceFE(lagrangian, Float64, order)
  Va, Vb = TestFESpace(tr[:Γ_a], lag), TestFESpace(tr[:Γ_b], lag)
  Y = MultiFieldFESpace([Va, Vb]); X = MultiFieldFESpace([TrialFESpace(Va), TrialFESpace(Vb)])
  fmap = Dict(:η_a => 1, :η_b => 2)
  q(x) = exp(-4 * ((x[1] - 3.0)^2 + (x[2] - 2.0)^2))        # load on plate a

  function solve_pair(c)
    ents = c === nothing ? (pa, pb) : (pa, pb, c)
    a(u, v) = (xu = FO.FieldMap(u, fmap); yv = FO.FieldMap(v, fmap);
               reduce(P._add_contribution, [P.stiffness(e, dom, xu, yv) for e in ents]))
    l((va, vb)) = ∫(q * va)dom[:dΓ_a]                        # plate b is unloaded
    solve(AffineFEOperator(a, l, X, Y))
  end

  # Continuous reference plate on the same mesh
  sc = G.StructureDomain(L = 4.0, W = 2.0, x₀ = [2.0, 1.0, 1.0], domain_symbol = :Γ_c)
  tank_c = mk_tank([sc]); tr_c = G.build_triangulations(tank_c, G.build_model(tank_c))
  dom_c = G.get_integration_domains(tr_c; degree = Dict(:Γ_c => 4, :default => 4))
  pc = plate(:η, :Γ_c)
  Vc = TestFESpace(tr_c[:Γ_c], lag)
  qa(x) = x[1] < 4.0 ? q(x) : 0.0
  ηc = solve(AffineFEOperator((u, v) -> P.stiffness(pc, dom_c, Dict(:η => u), Dict(:η => v)),
                              v -> ∫(qa * v)dom_c[:dΓ_c], TrialFESpace(Vc), Vc))
  pts_a = [Point(x, y, 1.0) for x in 2.2:0.3:3.9, y in 1.2:0.3:2.9]
  pts_b = [Point(x, y, 1.0) for x in 4.1:0.3:5.9, y in 1.2:0.3:2.9]
  ref = [ηc.(vec(pts_a)); ηc.(vec(pts_b))]
  sample(sol) = [sol[1].(vec(pts_a)); sol[2].(vec(pts_b))]
  rel(sol) = norm(sample(sol) - ref) / norm(ref)

  conn_(; kw...) = P.PlateConnection(; plate_a = pa, plate_b = pb, interface = :dΛ_ab, normal = :n_Λ_ab, kw...)

  @testset "rigid shear + rigid rotation ≈ continuous plate" begin
    e_rigid = rel(solve_pair(conn_()))
    e_hinge = rel(solve_pair(conn_(rotation = :free)))
    @test e_rigid < 1e-3
    @test e_hinge > 1e-3                 # the hinge changes the deflection
    @test e_rigid < 1e-2 * e_hinge
  end

  @testset "rotational spring interpolates between hinge and rigid" begin
    D_ρ = pa.C[1, 1, 1, 1]
    es = [rel(solve_pair(conn_(rotation = k))) for k in (0.0, 1e-1 * D_ρ, 1e1 * D_ρ, 1e4 * D_ρ)]
    @test issorted(es; rev = true)
    @test es[1] ≈ rel(solve_pair(conn_(rotation = :free))) rtol = 1e-10
    @test es[end] < 1e-2
  end

  @testset "no shear, free rotation: decoupled plates" begin
    sol = solve_pair(conn_(shear = 0.0, rotation = :free))
    @test norm(get_free_dof_values(sol[2])) < 1e-12 * norm(get_free_dof_values(sol[1]))
    sol0 = solve_pair(nothing)
    @test get_free_dof_values(sol[1]) ≈ get_free_dof_values(sol0[1])
  end

  @testset "field-less entity in build_problem" begin
    fe1 = PH.FESpaceConfig(vector_type = Vector{Float64})
    c = conn_(rotation = 1.0e3)
    physics = P.PhysicsParameters[P.PotentialFlow(dim = 3, fe = fe1), P.FreeSurface(dim = 3, fe = fe1), pa, pb, c]
    X2, Y2, fmap2 = SM.build_fe_spaces(physics, tr, PH.TimeDomainConfig())
    @test length(Y2.spaces) == 4
    @test sort(collect(keys(fmap2))) == sort([:ϕ, :κ, :η_a, :η_b])
    cfg = PH.TimeDomainConfig(tf = 0.2)
    tc = PH.TimeConfig(Δt = 0.1, tf = 0.2, u0 = [0.0, x -> 0.01 * cos(π * x[1] / 8), 0.0, 0.0],
                       u0t = zeros(4), u0tt = zeros(4))
    prob = SM.build_problem(mk_tank([sa, sb]; structure_connections = [conn]), physics, cfg; tconfig = tc)
    n = 0
    for (_, uh) in SM.simulate(prob, tc).solution
      n += 1
      @test all(isfinite, get_free_dof_values(uh))
    end
    @test n == 2
  end
end
