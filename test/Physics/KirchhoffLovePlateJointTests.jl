using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH

# KirchhoffLovePlate with line joints (JointLineDomain + JointRotationalSpring)
# against the hand-written C/DG plate + hinge form of the 3D FPV script.
#
# Plate x ∈ [2,6], y ∈ [1,3] on an 8 × 8 surface mesh (Δx = 1, Δy = 0.5),
# split by the lines x = 4 and y = 2 into a 2 × 2 floater grid.
@testset "KirchhoffLovePlate line joints" begin
  order = 2
  hinges = G.hinge_grid(x₀ = [2.0, 1.0, 1.0], a = 2.0, b = 1.0, nfx = 2, nfy = 2,
                        domain_symbol = :dΛh, normal_symbol = :n_Λh)
  tank = G.TankDomain(L = 8.0, W = 4.0, H = 1.0, nx = 8, ny = 8, nz = 2,
    structure_domains = [G.StructureDomain(L = 4.0, W = 2.0, x₀ = [2.0, 1.0, 1.0])],
    joint_domains = [hinges])
  trians = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(trians; degree = 2 * order)
  Γη = trians[:Γη]
  V = TestFESpace(Γη, ReferenceFE(lagrangian, Float64, order); conformity = :H1)
  U = TrialFESpace(V)

  E, ν, hb, ρ, g, γ = 1.0e9, 0.3, 0.2, 1025.0, 9.81, 6.0
  kᵣ = 3.0e4
  fe = PH.FESpaceConfig(order = order, vector_type = Vector{Float64}, γ = γ)
  plate(joints) = P.KirchhoffLovePlate(E = E, ν = ν, hb = hb, ρ = ρ, g = g, fe = fe,
                                       joints = joints)
  K(p) = assemble_matrix((η, v) -> P.stiffness(p, dom, Dict(:η => η), Dict(:η => v)), U, V)

  # Hand-written reference: C/DG terms only on facets inside floaters (:dΛη),
  # plus the hinge spring kᵣ/ρ [[∇v]]⋅[[∇η]] on the connection lines.
  p0 = plate(P.JointRotationalSpring[])
  C, D_ρ, h, n_Λ = p0.C, p0.C[1, 1, 1, 1], dom[:h_η], dom[:n_Λ_η]
  a_plate(η, v) = ∫(∇∇(v) ⊙ (C ⊙ ∇∇(η)) + g * v * η)dom[:dΓη] +
    ∫(-jump(∇(v)) ⊙ (mean(C ⊙ ∇∇(η)) ⋅ n_Λ.⁺) - (mean(C ⊙ ∇∇(v)) ⋅ n_Λ.⁺) ⊙ jump(∇(η)) +
      D_ρ * γ / h * jump(∇(v)) ⊙ jump(∇(η)))dom[:dΛη]
  a_hinge(η, v) = ∫(kᵣ / ρ * jump(∇(v)) ⋅ jump(∇(η)))dom[:dΛh]

  K_free_ref = assemble_matrix(a_plate, U, V)
  K_spring_ref = assemble_matrix((η, v) -> a_plate(η, v) + a_hinge(η, v), U, V)

  @testset "free hinge: joint facets excluded from C/DG" begin
    K_free = K(p0)
    @test norm(K_free - K_free_ref) <= 1e-12 * norm(K_free_ref)
    # a zero-stiffness spring is the same free hinge
    K_zero = K(plate([P.JointRotationalSpring(:dΛh, :n_Λh, 0.0)]))
    @test norm(K_zero - K_free_ref) <= 1e-12 * norm(K_free_ref)
  end

  @testset "rotational spring matches the script hinge form" begin
    K_spring = K(plate([P.JointRotationalSpring(:dΛh, :n_Λh, kᵣ / ρ)]))
    @test norm(K_spring - K_free_ref) > 1e-6 * norm(K_free_ref)
    @test norm(K_spring - K_spring_ref) <= 1e-10 * norm(K_spring_ref)
  end

  @testset "rigid limit recovers the continuous plate" begin
    # Continuous plate: same mesh without joints, C/DG on every facet.
    tank_c = G.TankDomain(L = 8.0, W = 4.0, H = 1.0, nx = 8, ny = 8, nz = 2,
      structure_domains = [G.StructureDomain(L = 4.0, W = 2.0, x₀ = [2.0, 1.0, 1.0])])
    trians_c = G.build_triangulations(tank_c, G.build_model(tank_c))
    dom_c = G.get_integration_domains(trians_c; degree = 2 * order)
    Vc = TestFESpace(trians_c[:Γη], ReferenceFE(lagrangian, Float64, order); conformity = :H1)
    Uc = TrialFESpace(Vc)
    # Load concentrated on the floater x ∈ [2,4], y ∈ [1,2] (a linear load
    # would give a curvature-free deflection that never activates the joints)
    q(x) = exp(-4 * ((x[1] - 3.0)^2 + (x[2] - 1.5)^2))
    η_static(p, d, U_, V_) = begin
      op = AffineFEOperator((η, v) -> P.stiffness(p, d, Dict(:η => η), Dict(:η => v)),
                            v -> ∫(q * v)d[:dΓη], U_, V_)
      solve(op)
    end
    # Identical meshes and plate cells → identical DOF numbering
    uc = get_free_dof_values(η_static(p0, dom_c, Uc, Vc))
    e(kr) = begin
      uj = get_free_dof_values(η_static(plate([P.JointRotationalSpring(:dΛh, :n_Λh, kr)]), dom, U, V))
      length(uj) == length(uc) || error("DOF layouts differ")
      norm(uj - uc) / norm(uc)
    end
    e_free, e_stiff = e(0.0), e(1.0e6 * D_ρ)
    @test e_free > 1e-4                 # the hinge makes a visible difference
    @test e_stiff < 1e-3 * e_free       # a very stiff spring closes the hinge
  end
end
