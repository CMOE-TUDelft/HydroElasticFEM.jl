using Test
using Gridap
using Gridap.CellData
using Gridap.FESpaces

import HydroElasticFEM.Physics as P
import HydroElasticFEM.Geometry as D
import HydroElasticFEM.Simulation.FEOperators as FO

# =========================================================================
# StructuralDampingZone: dashpot and Kelvin-Voigt damping zones on
# AbstractHydroelasticStructures, assembled structure-only (no fluid).
# =========================================================================

function _sdz_problem(; L=10.0, n=8, order=2, dims=1,
                      T=Float64)
  domain = dims == 1 ? (0.0, L) : (0.0, L, 0.0, L)
  model = CartesianDiscreteModel(domain, ntuple(_ -> n, dims))
  Ω = Triangulation(model)
  Λ = Skeleton(Ω)
  dom = D.IntegrationDomains(
    dΓη = Measure(Ω, 2 * order + 2),
    dΛη = Measure(Λ, 2 * order + 2),
    n_Λ_η = get_normal_vector(Λ),
    h_η = L / n,
  )
  reffe = ReferenceFE(lagrangian, Float64, order)
  V = TestFESpace(model, reffe; conformity=:H1, vector_type=Vector{T})
  U = TrialFESpace(V)
  return (; Ω, dom, X=MultiFieldFESpace([U]), Y=MultiFieldFESpace([V]))
end

# Matrix of `form(s, dom, x, y)` (or of the closure `form(x, y)`).
function _sdz_matrix(form, s, prob)
  fmap = Dict(P.variable_symbol(s) => 1)
  a((u,), (v,)) = form(s, prob.dom, FO.FieldMap((u,), fmap),
                       FO.FieldMap((v,), fmap))
  Matrix(assemble_matrix(a, prob.X, prob.Y))
end

_sdz_operator(weight) = (s, dom, x, y) -> P.stiffness_operator(s, dom, x, y, weight)

const _SDZ_L = 10.0
const _SDZ_ramp = P.ramp_profile(x_start=6.0, x_end=_SDZ_L)
const _SDZ_one = x -> 1.0

# Structures (1D unless noted) built with the given keyword arguments
_sdz_membrane(; kw...) = P.Membrane(L=_SDZ_L, mᵨ=0.5, Tᵨ=9.8; kw...)
_sdz_beam(; kw...) = P.EulerBernoulliBeam(L=_SDZ_L, mᵨ=0.5, EIᵨ=4.9; kw...)
_sdz_tbeam(; kw...) = P.TensionedEulerBernoulliBeam(
  L=_SDZ_L, mᵨ=0.5, EIᵨ=4.9, Tᵨ=9.8; kw...)
_sdz_plate(; kw...) = P.KirchhoffLovePlate(
  E=1.0e6, ν=0.3, hb=0.1, ρ=1.0, ρb=5.0, ambient_dim=2, manifold_dim=2;
  kw...)

const _SDZ_1D = (("Membrane", _sdz_membrane), ("EulerBernoulliBeam", _sdz_beam),
                 ("TensionedEulerBernoulliBeam", _sdz_tbeam))

@testset "ramp_profile" begin
  σ = P.ramp_profile(x_start=6.0, x_end=10.0)
  @test σ(VectorValue(3.0)) == 0.0
  @test σ(VectorValue(8.0)) ≈ 0.25
  @test σ(VectorValue(12.0)) == 1.0          # full strength beyond x_end
  σl = P.ramp_profile(x_start=4.0, x_end=0.0, power=1)   # lower-end zone
  @test σl(VectorValue(1.0)) ≈ 0.75
  @test σl(VectorValue(5.0)) == 0.0
  σd = P.ramp_profile(x_start=0.0, x_end=2.0, power=1,
                      direction=(0.0, 3.0))  # normalised internally
  @test σd(VectorValue(5.0, 1.0)) ≈ 0.5
  @test_throws ArgumentError P.ramp_profile(x_start=1.0, x_end=1.0)
end

@testset "No enabled zones: forms unchanged" begin
  prob = _sdz_problem()
  off = [P.StructuralDampingZone(profile=_SDZ_ramp, cᵨ=1.0, τ=0.1,
                                 enabled=false)]
  for (name, build) in _SDZ_1D
    @testset "$name" begin
      ref = build(τ=0.02)
      dis = build(τ=0.02, damping_zones=off)
      @test isempty(P.damping_zones(ref))
      @test _sdz_matrix(P.damping, dis, prob) == _sdz_matrix(P.damping, ref, prob)
      @test _sdz_matrix(_sdz_operator(nothing), ref, prob) ==
            _sdz_matrix(P.stiffness_operator, ref, prob)
      # Zones never touch the stiffness form
      on = build(τ=0.02, damping_zones=[P.StructuralDampingZone(
        profile=_SDZ_ramp, cᵨ=1.0, τ=0.1)])
      @test _sdz_matrix(P.stiffness, on, prob) == _sdz_matrix(P.stiffness, ref, prob)
    end
  end

  plate = _sdz_plate()
  @test !P.has_damping_form(plate)
  @test !P.has_damping_form(_sdz_plate(damping_zones=off))
  @test P.has_damping_form(_sdz_plate(damping_zones=[
    P.StructuralDampingZone(profile=_SDZ_ramp, cᵨ=1.0)]))
end

@testset "Uniform zones reproduce scalar damping" begin
  prob = _sdz_problem()
  τ0, c0 = 0.03, 0.7
  for (name, build) in _SDZ_1D
    @testset "$name" begin
      # σ ≡ 1 Kelvin-Voigt zone == scalar Rayleigh τ
      Cτ = _sdz_matrix(P.damping, build(τ=τ0), prob)
      Cz = _sdz_matrix(P.damping, build(damping_zones=[
        P.StructuralDampingZone(profile=_SDZ_one, τ=τ0)]), prob)
      @test Cz ≈ Cτ rtol=1e-12
      # σ ≡ 1 dashpot == (c/mᵨ) · mass matrix
      s = build(damping_zones=[P.StructuralDampingZone(profile=_SDZ_one, cᵨ=c0)])
      M = _sdz_matrix(P.mass, s, prob)
      @test _sdz_matrix(P.damping, s, prob) ≈ (c0 / P.mass_density(s)) * M rtol=1e-12
      # Zones add up: two half-strength zones == one full zone
      half = P.StructuralDampingZone(profile=_SDZ_ramp, cᵨ=c0 / 2, τ=τ0 / 2)
      full = P.StructuralDampingZone(profile=_SDZ_ramp, cᵨ=c0, τ=τ0)
      @test _sdz_matrix(P.damping, build(damping_zones=[half, half]), prob) ≈
            _sdz_matrix(P.damping, build(damping_zones=[full]), prob) rtol=1e-12
    end
  end

  # Kirchhoff-Love plate (2D manifold): dashpot and constant-weight operator
  prob2 = _sdz_problem(n=4, dims=2)
  plate = _sdz_plate(damping_zones=[P.StructuralDampingZone(
    profile=_SDZ_one, cᵨ=c0, τ=τ0)])
  M = _sdz_matrix(P.mass, plate, prob2)
  K = _sdz_matrix(P.stiffness_operator, plate, prob2)
  @test _sdz_matrix(P.damping, plate, prob2) ≈
        (c0 / P.mass_density(plate)) * M + τ0 * K rtol=1e-10
end

@testset "Ramp weight enters the material coefficients" begin
  prob = _sdz_problem()
  τ0 = 0.05
  zone = [P.StructuralDampingZone(profile=_SDZ_ramp, τ=τ0)]

  # Bending (bulk and C/DG skeleton): EI(x)·σ(x) as a variable rigidity
  Cz = _sdz_matrix(P.damping, _sdz_beam(damping_zones=zone), prob)
  Kσ = _sdz_matrix(P.stiffness_operator,
                   _sdz_beam(EIᵨ=x -> 4.9 * _SDZ_ramp(x)), prob)
  @test Cz ≈ τ0 * Kσ rtol=1e-12

  # Tension: ∫ τ σ T ∇v⋅∇ηₜ
  Cm = _sdz_matrix(P.damping, _sdz_membrane(damping_zones=zone), prob)
  dΩ = prob.dom[:dΓη]
  σ = CellField(_SDZ_ramp, prob.Ω)
  a((u,), (v,)) = ∫(τ0 * 9.8 * σ * ∇(v) ⋅ ∇(u))dΩ
  @test Cm ≈ Matrix(assemble_matrix(a, prob.X, prob.Y)) rtol=1e-12

  # Tensioned beam = bending part + tension part
  Ct = _sdz_matrix(P.damping, _sdz_tbeam(damping_zones=zone), prob)
  @test Ct ≈ Cz + Cm rtol=1e-12
end

@testset "Frequency domain: τ-zone == complex rigidity" begin
  # The L1 sponge EI(x) = EI(1 − iωτσ(x)) is exactly a Kelvin-Voigt zone.
  prob = _sdz_problem(T=ComplexF64)
  ω, τ0 = 2.5, 0.08
  complex_beam = _sdz_beam(EIᵨ=x -> 4.9 * (1 - im * ω * τ0 * _SDZ_ramp(x)))
  zone_beam = _sdz_beam(damping_zones=[
    P.StructuralDampingZone(profile=_SDZ_ramp, τ=τ0)])
  weak(s, dom, x, y) = P.weakform(s, dom, ω, x, y)
  A_complex = _sdz_matrix(weak, complex_beam, prob)
  A_zone = _sdz_matrix(weak, zone_beam, prob)
  @test A_zone ≈ A_complex rtol=1e-12
end
