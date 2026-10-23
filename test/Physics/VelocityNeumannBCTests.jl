using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.AssemblyContexts as AC

# PrescribedInletPotentialBC(quantity = :velocity) must reproduce the
# hand-written Neumann flux ∫ w (v_in ⋅ nΓ) dΓ, with an optional mask.
@testset "PrescribedInletPotentialBC — velocity forcing (3D)" begin
  L, W, H = 8.0, 4.0, 1.0
  tank = G.TankDomain(L = L, W = W, H = H, nx = 8, ny = 4, nz = 2)
  trians = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(trians; degree = 2)
  V = TestFESpace(trians[:Ω], ReferenceFE(lagrangian, Float64, 1); conformity = :H1)

  a, b, c = 0.7, -0.4, 0.2
  vin(x, t) = VectorValue(a * (1 + t), b + x[1], c * x[3])
  χ(x) = x[1] < 4.0 ? 1.0 : 0.0
  t = 0.3
  ctx = AC.TimeAssemblyContext(dom, t, 0.0)

  nin = get_normal_vector(trians[:Γin])
  nlat = get_normal_vector(trians[:Γlateral])
  dΓlat = Measure(trians[:Γlateral], 2)

  @testset "inlet (normal from the wall triangulation)" begin
    bc = P.PrescribedInletPotentialBC(domain = :dΓin, forcing = (t -> (x -> vin(x, t))),
                                      quantity = :velocity)
    pf = P.PotentialFlow(dim = 3, boundary_conditions = [bc])
    f_pkg = assemble_vector(w -> P._rhs_bc_contribution(pf, bc, ctx, w), V)
    f_ref = assemble_vector(w -> ∫(w * ((x -> vin(x, t)) ⋅ nin))dom[:dΓin], V)
    @test norm(f_pkg - f_ref) <= 1e-12 * norm(f_ref)
    # Inlet at x = 0 has outward normal -e_x: ∫ v⋅n = -a(1+t) W H
    @test sum(f_pkg) ≈ -a * (1 + t) * W * H
  end

  @testset "lateral walls with mask" begin
    # Standard 3D measure :dΓlateral with its stored normal :nΓlateral
    bc = P.PrescribedInletPotentialBC(domain = :dΓlateral, forcing = (t -> (x -> vin(x, t))),
                                      quantity = :velocity, mask = χ)
    pf = P.PotentialFlow(dim = 3, boundary_conditions = [bc])
    f_pkg = assemble_vector(w -> P._rhs_bc_contribution(pf, bc, ctx, w), V)
    f_ref = assemble_vector(w -> ∫(w * χ * ((x -> vin(x, t)) ⋅ nlat))dΓlat, V)
    @test norm(f_ref) > 0
    @test norm(f_pkg - f_ref) <= 1e-12 * norm(f_ref)
  end

  @testset "explicit normal key takes precedence" begin
    dom_n = G.IntegrationDomains(merge(dom.data,
      Dict{Symbol,Any}(:dΓlat => dΓlat, :nΓlat => VectorValue(0.0, 1.0, 0.0))))
    ctx_n = AC.TimeAssemblyContext(dom_n, t, 0.0)
    bc = P.PrescribedInletPotentialBC(domain = :dΓlat, forcing = (x -> VectorValue(0.0, 2.0, 0.0)),
                                      quantity = :velocity)
    pf = P.PotentialFlow(dim = 3, boundary_conditions = [bc])
    f_pkg = assemble_vector(w -> P._rhs_bc_contribution(pf, bc, ctx_n, w), V)
    @test sum(f_pkg) ≈ 2.0 * 2 * L * H   # both walls, same (fixed) normal
  end

  @testset "mask on a scalar Neumann forcing" begin
    bc = P.PrescribedInletPotentialBC(domain = :dΓbot, forcing = (x -> 1.5),
                                      quantity = :normal_gradient, mask = χ)
    pf = P.PotentialFlow(dim = 3, boundary_conditions = [bc])
    f_pkg = assemble_vector(w -> P._rhs_bc_contribution(pf, bc, ctx, w), V)
    @test sum(f_pkg) ≈ 1.5 * 4.0 * W
  end
end
