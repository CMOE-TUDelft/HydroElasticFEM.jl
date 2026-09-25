using Test

using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.FESpaces
using LinearAlgebra: norm

import HydroElasticFEM.Physics as P
import HydroElasticFEM.Geometry as D
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.ParameterHandler as FES

# =========================================================================
# TensionedEulerBernoulliBeam weak form integration tests
#
#   1. Membrane limit  (EIᵨ = 0): stiffness/damping matrices match Membrane
#      exactly on a shared mesh.
#   2. Beam limit      (Tᵨ = 0): stiffness/damping matrices match
#      EulerBernoulliBeam exactly on a shared mesh.
#   3. Dry beam frequency shift: adding tension strictly increases the
#      analytical (simply-supported) natural frequency, matching the
#      struct's own `ωn1` field.
#   4. Hydroelastic regression: the existing Khabakhpasheva floating-beam
#      benchmark (which exercises the refactored EulerBernoulliBeam through
#      AbstractHydroelasticStructure) still converges to a sane response.
# =========================================================================

@testset "TensionedEulerBernoulliBeam weak forms" begin

  # -----------------------------------------------------------------------
  # Helper: build 1D beam mesh, FE spaces, and IntegrationDomains
  # (mirrors EulerBernoulliBeamWeakFormTests.jl's build_beam_problem)
  # -----------------------------------------------------------------------

  function build_beam_problem(; L, nel, order)
    model = CartesianDiscreteModel((0, L), (nel,))
    Ω  = Triangulation(model)
    Λ  = Skeleton(Ω)
    Λb = Boundary(model, tags="boundary")
    h  = L / nel

    dom = D.IntegrationDomains(
      dΓη    = Measure(Ω, 2 * order + 2),
      dΛη    = Measure(Λ, 2 * order + 2),
      n_Λ_η  = get_normal_vector(Λ),
      h_η    = h,
      dΛ_sb  = Measure(Λb, 2 * order + 2),
      n_Λ_sb = get_normal_vector(Λb),
    )

    reffe = ReferenceFE(lagrangian, Float64, order)
    V = TestFESpace(model, reffe, dirichlet_tags="boundary",
                    vector_type=Vector{Float64})
    U = TrialFESpace(V, 0.0)
    Y = MultiFieldFESpace([V])
    X = MultiFieldFESpace([U])

    return (; model, Ω, dom, X, Y, h)
  end

  # -----------------------------------------------------------------------
  # Helpers: assemble a structure's stiffness/damping/mass matrix on a
  # shared mesh, so two different AbstractHydroelasticStructure instances
  # can be compared directly.
  # -----------------------------------------------------------------------

  function assemble_form(form::Function, structure, prob)
    sym  = P.variable_symbol(structure)
    fmap = Dict(sym => 1)
    a((u,), (v,)) = form(structure, prob.dom,
                         FO.FieldMap((u,), fmap), FO.FieldMap((v,), fmap))
    return assemble_matrix(a, prob.X, prob.Y)
  end

  assemble_stiffness(structure, prob) = assemble_form(P.stiffness, structure, prob)
  assemble_damping(structure, prob)   = assemble_form(P.damping, structure, prob)
  assemble_mass(structure, prob)      = assemble_form(P.mass, structure, prob)

  # Shared mesh/FE spaces for all exact-limit comparisons below.
  L, order, nel = 1.0, 2, 40
  fe   = FES.FESpaceConfig(order=order, vector_type=Vector{Float64})
  prob = build_beam_problem(L=L, nel=nel, order=order)

  mᵨ, EIᵨ, Tᵨ, τ = 1.0, 100.0, 50.0, 0.02

  # -----------------------------------------------------------------------
  # Test 1 — Membrane limit (EIᵨ = 0)
  # -----------------------------------------------------------------------
  @testset "Membrane limit" begin

    tbeam    = P.TensionedEulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=0.0, Tᵨ=Tᵨ, τ=τ, fe=fe)
    membrane = P.Membrane(L=L, mᵨ=mᵨ, Tᵨ=Tᵨ, τ=τ, fe=fe)

    @test tbeam isa P.AbstractHydroelasticStructure
    @test tbeam isa P.Structure

    K1, K2 = assemble_stiffness(tbeam, prob), assemble_stiffness(membrane, prob)
    @test norm(K1 - K2) / norm(K2) < 1e-12

    C1, C2 = assemble_damping(tbeam, prob), assemble_damping(membrane, prob)
    @test norm(C1 - C2) / norm(C2) < 1e-12

    M1, M2 = assemble_mass(tbeam, prob), assemble_mass(membrane, prob)
    @test norm(M1 - M2) / norm(M2) < 1e-12

  end

  # -----------------------------------------------------------------------
  # Test 2 — Beam limit (Tᵨ = 0)
  # -----------------------------------------------------------------------
  @testset "Beam limit" begin

    tbeam = P.TensionedEulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=EIᵨ, Tᵨ=0.0, τ=τ, fe=fe)
    euler = P.EulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=EIᵨ, τ=τ, fe=fe)

    K1, K2 = assemble_stiffness(tbeam, prob), assemble_stiffness(euler, prob)
    @test norm(K1 - K2) / norm(K2) < 1e-12

    C1, C2 = assemble_damping(tbeam, prob), assemble_damping(euler, prob)
    @test norm(C1 - C2) / norm(C2) < 1e-12

    M1, M2 = assemble_mass(tbeam, prob), assemble_mass(euler, prob)
    @test norm(M1 - M2) / norm(M2) < 1e-12

  end

  # -----------------------------------------------------------------------
  # Test 3 — Dry beam frequency shift
  #
  # Analytical dispersion relation for the dry (in-vacuo), simply supported
  # beam:  ω_n² = (EI/m) k_n⁴ + (T/m) k_n²,  k_n = nπ/L.
  # Adding tension must strictly increase the fundamental frequency.
  # -----------------------------------------------------------------------
  @testset "Dry beam frequency shift" begin

    k1 = π / L
    ω_no_tension = sqrt(EIᵨ / mᵨ * k1^4)
    ω_tension    = sqrt(EIᵨ / mᵨ * k1^4 + Tᵨ / mᵨ * k1^2)

    @test ω_tension > ω_no_tension

    # The struct's own derived ωn1 field must match the analytical formula.
    beam = P.TensionedEulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=EIᵨ, Tᵨ=Tᵨ)
    @test beam.ωn1 ≈ ω_tension

    # Membrane limit (EIᵨ=0) and beam limit (Tᵨ=0) of ωn1 recover the
    # respective single-mechanism structures' own ωn1 formulas.
    beam_no_EI = P.TensionedEulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=0.0, Tᵨ=Tᵨ)
    membrane   = P.Membrane(L=L, mᵨ=mᵨ, Tᵨ=Tᵨ)
    @test beam_no_EI.ωn1 ≈ membrane.ωn1

  end

  # -----------------------------------------------------------------------
  # Test 4 — Hydroelastic regression
  #
  # Reuses the existing Khabakhpasheva floating-beam benchmark (which
  # exercises EulerBernoulliBeam's mass/damping/stiffness/rhs — now
  # inherited from AbstractHydroelasticStructure — inside the full
  # monolithic fluid-structure solve) at the test suite's standard coarse
  # resolution, to confirm the refactor has not regressed the existing
  # hydroelastic pathway.
  # -----------------------------------------------------------------------
  @testset "Hydroelastic regression (Khabakhpasheva benchmark)" begin

    isdefined(@__MODULE__, :KhabakhpashevaBeamJointExample) ||
      include(joinpath(@__DIR__, "..", "..", "examples", "KhabakhpashevaBeamJointExample.jl"))
    using .KhabakhpashevaBeamJointExample

    p0 = KhabakhpashevaCaseParams(
      name       = "test_tensioned_beam_refactor_regression",
      nx         = 2,
      ny         = 1,
      order      = 2,
      ξ          = 0.0,
      vtk_output = false,
      make_plot  = false,
    )

    xs, eta, meta = run_khabakhpasheva_case(p0)

    @test all(isfinite, eta)
    @test length(xs) == length(eta)
    @test all(eta .>= 0.0)
    @test maximum(eta) > 0.1
    @test maximum(eta) < 20.0

  end

end
