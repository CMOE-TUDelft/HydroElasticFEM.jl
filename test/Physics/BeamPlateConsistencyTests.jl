using Test
using Gridap

import HydroElasticFEM.Physics as P
import HydroElasticFEM.Geometry as D
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.ParameterHandler as FES

@testset "Beam-plate 1D consistency" begin
  ρw = 1025.0
  E = 11.9e9
  ν = 0.13
  hb = 2.0
  L = 20.0

  EIρ = E * hb^3 / (12 * ρw)

  beam = P.EulerBernoulliBeam(
    L = L,
    mᵨ = 1.0,
    EIᵨ = EIρ,
  )
  plate = P.KirchhoffLovePlate(
    E = E,
    ν = ν,
    hb = hb,
    ρ = ρw,
    ambient_dim = 2,
    manifold_dim = 1,
  )

  @test P.equivalent_beam_rigidity(plate; width = 1.0) ≈ beam.EIᵨ atol = 1e-10
  @test plate.C[1, 1, 1, 1] ≈ beam.EIᵨ atol = 1e-10
  @test plate.C[1, 1, 2, 2] ≈ 0.0 atol = 1e-12
end

@testset "KirchhoffLovePlate 2D tensor — major symmetry" begin
  # 2D plane-stress tensor: ambient_dim=2, manifold_dim=2
  E  = 11.9e9
  ν  = 0.13
  hb = 2.0
  ρw = 1025.0
  C2 = P.build_kl_tensor(2, 2, E, ν, hb, ρw)

  @test P.check_major_symmetry(C2; atol = 1e-10, dim = 2)

  # Minor symmetry
  for i in 1:2, j in 1:2, k in 1:2, l in 1:2
    @test C2[i, j, k, l] ≈ C2[j, i, k, l]  atol = 1e-12
    @test C2[i, j, k, l] ≈ C2[i, j, l, k]  atol = 1e-12
  end

  # Diagonal bending rigidity C[1,1,1,1] = D/ρ
  D_rho = E * hb^3 / (12 * (1 - ν^2) * ρw)
  @test C2[1, 1, 1, 1] ≈ D_rho  atol = 1e-8
  @test C2[2, 2, 2, 2] ≈ D_rho  atol = 1e-8

  # Cross-term: C[1,1,2,2] = ν·D/ρ
  @test C2[1, 1, 2, 2] ≈ ν * D_rho  atol = 1e-8
end

@testset "TimoshenkoBeam → EulerBernoulli thin limit" begin
  # Timoshenko beam converges to Euler-Bernoulli as h/L → 0.
  # For a simply-supported beam with uniform load q:
  #   w_TB = 5qL^4/(384EI) + qL^2/(8κGA)
  #   w_EB = 5qL^4/(384EI)
  # The shear term scales as (h/L)^2; at h/L = 0.005 it is < 0.01 %.
  E, ν, κ = 210e9, 0.3, 5 / 6
  b, L, q = 0.05, 1.0, 1e4
  h = 0.005 * L     # h/L = 0.005

  I_val = b * h^3 / 12
  EI    = E * I_val

  w_EB = 5 * q * L^4 / (384 * EI)

  # Build and solve on a 1-D mesh (ρ_w=1.0, no normalisation)
  model = CartesianDiscreteModel((0.0, L), (20,))
  Ω  = Triangulation(model)
  dΩ = Measure(Ω, 8)

  dom = D.IntegrationDomains(dΓη = dΩ)

  beam = P.TimoshenkoBeam(
    E=E, ν=ν, h_beam=h, b_beam=b, ρ_s=1.0, ρ_w=1.0, g=0.0, κ=κ,
    tangent = VectorValue(1.0),
    fe_w = FES.FESpaceConfig(order=2, vector_type=Vector{Float64}),
    fe_θ = FES.FESpaceConfig(order=1, vector_type=Vector{Float64}),
  )

  reffe_w = ReferenceFE(lagrangian, Float64, 2)
  V_w = TestFESpace(model, reffe_w;
    conformity=:H1, dirichlet_tags="boundary", vector_type=Vector{Float64})
  U_w = TrialFESpace(V_w, 0.0)
  reffe_θ = ReferenceFE(lagrangian, Float64, 1)
  V_θ = TestFESpace(model, reffe_θ;
    conformity=:H1, vector_type=Vector{Float64})
  U_θ = TrialFESpace(V_θ)

  X = MultiFieldFESpace([U_w, U_θ])
  Y = MultiFieldFESpace([V_w, V_θ])

  fmap = Dict(:w => 1, :θ => 2)
  src_w(x) = q
  fmap_rhs = Dict(:w => 1)

  a((w, θ), (v_w, v_θ)) = P.stiffness(beam, dom,
    FO.FieldMap((w, θ), fmap),
    FO.FieldMap((v_w, v_θ), fmap))
  l((v_w, v_θ)) = P.rhs(beam, dom,
    FO.FieldMap((src_w,), fmap_rhs),
    FO.FieldMap((v_w, v_θ), fmap))

  op = AffineFEOperator(a, l, X, Y)
  uh = solve(LUSolver(), op)
  w_TB = uh[1](Point(L / 2))

  @test abs(w_TB - w_EB) / abs(w_EB) < 0.01   # <1 % from EB at h/L=0.005
end
