using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH

@testset "KirchhoffLovePlate per-structure skeleton keys" begin
  @test G.skeleton_keys(:Γη) == (:dΛη, :n_Λ_η, :h_η)
  @test G.skeleton_keys(:Γ_a) == (:dΛ_Γ_a, :n_Λ_Γ_a, :h_Γ_a)

  order = 2
  fe = PH.FESpaceConfig(order = order, vector_type = Vector{Float64}, γ = 6.0)
  plate(sym) = P.KirchhoffLovePlate(E = 1e9, ν = 0.3, hb = 0.2, fe = fe, space_domain_symbol = sym)
  K(p, dom, V) = assemble_matrix((η, v) -> P.stiffness(p, dom, Dict(:η => η), Dict(:η => v)),
                                 TrialFESpace(V), V)
  lag = ReferenceFE(lagrangian, Float64, order)

  @testset "single plate: own keys reproduce :Γη" begin
    sa = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [2.0, 1.0, 1.0], domain_symbol = :Γ_a)
    tank = G.TankDomain(L = 8.0, W = 4.0, H = 1.0, nx = 8, ny = 8, nz = 2, structure_domains = [sa])
    tr = G.build_triangulations(tank, G.build_model(tank))
    dom = G.get_integration_domains(tr; degree = Dict(:Γη => 4, :Γ_a => 4, :default => 4))
    for k in (:dΛ_Γ_a, :n_Λ_Γ_a, :h_Γ_a)
      @test haskey(dom, k)
    end
    @test dom[:h_Γ_a] ≈ dom[:h_η]
    Kη = K(plate(:Γη), dom, TestFESpace(tr[:Γη], lag))
    Ka = K(plate(:Γ_a), dom, TestFESpace(tr[:Γ_a], lag))
    @test norm(Ka - Kη) <= 1e-12 * norm(Kη)
  end

  @testset "two plates: each uses only its own facets" begin
    # Adjacent plates sharing the edge x = 4: the shared edge belongs to
    # neither per-structure skeleton (it is an interface, not interior).
    sa = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [2.0, 1.0, 1.0], domain_symbol = :Γ_a)
    sb = G.StructureDomain(L = 2.0, W = 2.0, x₀ = [4.0, 1.0, 1.0], domain_symbol = :Γ_b)
    tank = G.TankDomain(L = 8.0, W = 4.0, H = 1.0, nx = 8, ny = 8, nz = 2, structure_domains = [sa, sb])
    tr = G.build_triangulations(tank, G.build_model(tank))
    dom = G.get_integration_domains(tr; degree = Dict(:Γ_a => 4, :Γ_b => 4, :default => 4))
    # 2 × 4 cells each: interior facets x = 3 (length 2) and y = 1.5, 2, 2.5 (3 × 2)
    @test sum(∫(1.0)dom[:dΛ_Γ_a]) ≈ 8.0
    @test sum(∫(1.0)dom[:dΛ_Γ_b]) ≈ 8.0
    @test sum(∫(1.0)dom[:dΛη]) ≈ 8.0 + 8.0 + 2.0   # the global skeleton also has the shared edge

    pa, pb = plate(:Γ_a), plate(:Γ_b)
    Va = TestFESpace(tr[:Γ_a], lag)
    Ka = K(pa, dom, Va)
    # hand-written C/DG form on the plate-a skeleton only
    C, D_ρ = pa.C, pa.C[1, 1, 1, 1]
    n, h, dΛ = dom[:n_Λ_Γ_a], dom[:h_Γ_a], dom[:dΛ_Γ_a]
    a(η, v) = ∫(∇∇(v) ⊙ (C ⊙ ∇∇(η)) + pa.g * v * η)dom[:dΓ_a] +
      ∫(-jump(∇(v)) ⊙ (mean(C ⊙ ∇∇(η)) ⋅ n.⁺) - (mean(C ⊙ ∇∇(v)) ⋅ n.⁺) ⊙ jump(∇(η)) +
        D_ρ * 6.0 / h * jump(∇(v)) ⊙ jump(∇(η)))dΛ
    Kref = assemble_matrix(a, TrialFESpace(Va), Va)
    @test norm(Ka - Kref) <= 1e-12 * norm(Kref)
    # plate b sees a translated copy of the same problem: identical matrix
    Kb = K(pb, dom, TestFESpace(tr[:Γ_b], lag))
    @test size(Kb) == size(Ka)
    @test norm(sort(abs.(Kb.nzval)) - sort(abs.(Ka.nzval))) <= 1e-10 * norm(Ka.nzval)
  end
end
