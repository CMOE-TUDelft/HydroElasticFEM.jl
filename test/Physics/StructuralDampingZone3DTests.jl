using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM

# =========================================================================
# StructuralDampingZone on a KirchhoffLovePlate in 3D.
#
# A 3D channel one (square) cell wide, driven by a width-uniform inlet flux,
# covered by a plate with ν = 0 and a ramped zone (dashpot + Kelvin-Voigt).
# With ν = 0 the free side edges carry no moment for a width-uniform field,
# so the 3D solution is uniform across the width and, on the tensor-product
# mesh, equal to the 2D solution with an EulerBernoulliBeam of rigidity
# EI = D/ρ, the same mass and the same zone (square cells give the same C/DG
# element size h). This checks the plate's dashpot and its weighted
# Hessian C/DG operator (bulk, skeleton, penalty) to round-off.
# (With ν ≠ 0 the free edges bend anticlastically and the fields differ.)
# =========================================================================

@testset "3D plate zone == 2D beam zone (width-uniform, ν = 0)" begin
  g, H, m, D, ω = 9.81, 5.0, 0.01, 49.0, 2.0
  λ = 2π / 0.3849
  L, nx, nz = 3λ, 24, 4
  dx = L / nx
  hb = 0.01
  zone = P.StructuralDampingZone(
    profile=P.ramp_profile(x_start=L - λ, x_end=L), cᵨ=3.0, τ=0.1)
  fe(; kw...) = PH.FESpaceConfig(; order=2, vector_type=Vector{ComplexF64},
                                 kw...)
  inlet = P.PrescribedInletPotentialBC(domain=:dΓin, quantity=:normal_gradient,
                                       forcing=x -> 1.0 + 0im)
  xs = collect(range(0.3λ, L, length=15))   # includes the zone

  tank2 = G.TankDomain(L=L, H=H, nx=nx, ny=nz,
    structure_domains=[G.StructureDomain(L=L, x₀=[0.0, H])])
  beam = P.EulerBernoulliBeam(L=L, mᵨ=m, EIᵨ=D, g=g, fe=fe(γ=6.0),
                              damping_zones=[zone])
  r2 = SM.simulate(SM.build_problem(tank2,
    P.PhysicsParameters[P.PotentialFlow(g=g, fe=fe(),
                                        boundary_conditions=[inlet]), beam],
    SM.FreqDomainConfig(ω=ω)))
  η2 = r2.solution[r2.fmap[:η_b]]
  v2 = ComplexF64[η2(Point(x, H)) for x in xs]

  tank3 = G.TankDomain(L=L, W=dx, H=H, nx=nx, ny=1, nz=nz,
    structure_domains=[G.StructureDomain(L=L, W=dx, x₀=[0.0, 0.0, H])])
  plate = P.KirchhoffLovePlate(E=12 * D / hb^3, ν=0.0, hb=hb, ρ=1.0,
                               ρb=m / hb, g=g, fe=fe(γ=6.0),
                               damping_zones=[zone])
  @test plate.C[1, 1, 1, 1] ≈ D
  @test P.mass_density(plate) ≈ m
  @test P.has_damping_form(plate)
  r3 = SM.simulate(SM.build_problem(tank3,
    P.PhysicsParameters[P.PotentialFlow(g=g, dim=3, fe=fe(),
                                        boundary_conditions=[inlet]), plate],
    SM.FreqDomainConfig(ω=ω)))
  η3 = r3.solution[r3.fmap[:η]]
  v3 = ComplexF64[η3(Point(x, dx / 2, H)) for x in xs]
  v3_edge = ComplexF64[η3(Point(x, 0.1 * dx, H)) for x in xs]

  @test norm(v3_edge - v3) / norm(v3) < 1e-10   # uniform across the width
  @test norm(v3 - v2) / norm(v2) < 1e-10
end
