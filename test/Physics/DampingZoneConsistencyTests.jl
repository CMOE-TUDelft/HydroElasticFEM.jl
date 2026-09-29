using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.AssemblyContexts as AC
import HydroElasticFEM.Simulation.FEOperators as FO

# A DampingZoneBC relaxes (∂ₙϕ, κ) towards (vz_in, η_in).  When the discrete
# fields already equal the target, the zone terms of the stiffness forms must
# cancel the zone terms of the right-hand side exactly (in both the u- and the
# αₕ·w-weighted parts).
@testset "DampingZoneBC target consistency" begin
  H, L, Ld = 5.0, 20.0, 10.0
  tank = G.TankDomain(L = L, H = H, nx = 20, ny = 4,
    damping_zones = [G.DampingZone(L = Ld, x₀ = [0.0, H], domain_symbol = :Γd)])
  trians = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(trians; degree = 2)

  βₕ, g, αₕ = 0.5, 9.81, 0.7
  c, η₀ = 0.3, 0.05
  bc = P.DampingZoneBC(domain = :dΓd_1, μ₁ = (x -> 2.0), μ₂ = (x -> 1.5),
                       η_in = (x -> η₀), vz_in = (x -> c), βₕ = βₕ)
  pf = P.PotentialFlow(g = g, boundary_conditions = [bc])
  fs = P.FreeSurface(g = g, βₕ = βₕ)

  reffe = ReferenceFE(lagrangian, Float64, 1)
  Vϕ = TestFESpace(trians[:Ω], reffe; conformity = :H1)
  Vκ = TestFESpace(trians[:Γκ], reffe; conformity = :H1)
  Y = MultiFieldFESpace([Vϕ, Vκ])
  X = MultiFieldFESpace([TrialFESpace(Vϕ), TrialFESpace(Vκ)])
  fmap = Dict(:ϕ => 1, :κ => 2)

  # ∂ₙϕ = ∂ϕ/∂y = c on the top surface, κ = η₀: exactly the relaxation target
  xh = interpolate_everywhere([x -> c * x[2], x -> η₀], X)
  x = FO.FieldMap(xh, fmap)
  ctx = AC.TimeAssemblyContext(dom, 0.0, αₕ)

  function zone_residual(yb)
    y = FO.FieldMap(yb, fmap)
    w, κ = y[:ϕ], x[:κ]
    # PF/FS stiffness minus its free-surface stabilisation part = zone terms
    stab = ∫(βₕ * g * αₕ * w * κ)dom[:dΓκ]
    lhs = P._stiffness_bc_contribution(pf, bc, ctx, x[:ϕ], w) +
          P.stiffness(pf, fs, ctx, x, y) - stab
    rhs = P._rhs_bc_contribution(pf, bc, ctx, w) + P.rhs(pf, fs, ctx, nothing, y)
    lhs - rhs
  end
  rhs_only(yb) = (y = FO.FieldMap(yb, fmap);
    P._rhs_bc_contribution(pf, bc, ctx, y[:ϕ]) + P.rhs(pf, fs, ctx, nothing, y))

  r = assemble_vector(zone_residual, Y)
  f = assemble_vector(rhs_only, Y)
  @test norm(f) > 0
  @test norm(r) <= 1e-12 * norm(f)
end
