using Gridap
using Gridap.Geometry
using Gridap.CellData
using SparseArrays

import HydroElasticFEM.PhysicsCore.Entities as Entities
import HydroElasticFEM.Geometry as Geometry
import HydroElasticFEM.Simulation.FEOperators as WF

# =========================================================================
# Mini mesh setup: 50m × 10m tank, membrane from x=15 to x=35
# =========================================================================

@testset "FEOperators" begin

  order = 1
  nx, ny = 20, 4
  domain = (0, 50, -10, 0)
  model = CartesianDiscreteModel(domain, (nx, ny))

  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω, "surface", [3, 4, 6])
  add_tag_from_tags!(labels_Ω, "bottom",  [1, 2, 5])
  add_tag_from_tags!(labels_Ω, "inlet",   [7])
  add_tag_from_tags!(labels_Ω, "outlet",  [8])
  add_tag_from_tags!(labels_Ω, "water",   [9])

  Ω = Interior(model)
  Γ = Boundary(model, tags="surface")
  Γin = Boundary(model, tags="inlet")
  Γot = Boundary(model, tags="outlet")

  # Surface masking
  xm₀, xm₁ = 15.0, 35.0
  is_mem(xs) = let n=length(xs); x=(1/n)*sum(xs); (xm₀ <= x[1] <= xm₁) * (x[2] ≈ 0.0) end

  xΓ = get_cell_coordinates(Γ)
  Γm_mask = lazy_map(is_mem, xΓ)
  Γm  = Triangulation(Γ, findall(Γm_mask))
  Γfs = Triangulation(Γ, findall(!, Γm_mask))
  Γη  = Triangulation(Γ, findall(Γm_mask))
  Γκ  = Triangulation(Γ, findall(!, Γm_mask))

  degree = 2 * order
  dΩ   = Measure(Ω, degree)
  dΓm  = Measure(Γm, degree)
  dΓfs = Measure(Γfs, degree)
  dΓin = Measure(Γin, degree)
  dΓot = Measure(Γot, degree)

  # FE spaces (3-field: phi, kappa, eta)
  reffe = ReferenceFE(lagrangian, Float64, order)
  V_Ω  = TestFESpace(Ω,  reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
  U_Ω  = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  U_Γη = TrialFESpace(V_Γη)
  X = MultiFieldFESpace([U_Ω, U_Γκ, U_Γη])
  Y = MultiFieldFESpace([V_Ω, V_Γκ, V_Γη])

  # Parameters
  ω   = 2.0
  ρw  = 1025.0
  g_  = 9.81
  βₕ  = 0.5

  # Physics entities
  fluid = Entities.PotentialFlow(ρw=ρw, g=g_)
  fsurf = Entities.FreeSurface(ρw=ρw, g=g_, βₕ=βₕ)
  mem   = Entities.Membrane2D(L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1, g=g_)

  # Field mapping
  fmap = Dict(:ϕ => 1, :κ => 2, :η_m => 3)

  dom = Geometry.IntegrationDomains(dΩ=dΩ, dΓ_fs=dΓfs, dΓ_s=dΓm,
                                dΓ_in=dΓin, dΓ_ot=dΓot)

  # =========================================================================
  # Test IntegrationDomains + FieldMap construction
  # =========================================================================

  @testset "IntegrationDomains construction" begin
    d = Geometry.IntegrationDomains(dΩ=dΩ, dΓ_s=dΓm)
    @test d[:dΩ] === dΩ
    @test haskey(d, :dΓ_s)
    @test !haskey(d, :dΓ_fs)
  end

  @testset "FieldMap construction" begin
    fd = WF.FieldMap((1, 2, 3), fmap)
    @test fd[:ϕ] == 1
    @test fd[:κ] == 2
    @test fd[:η_m] == 3
    @test haskey(fd, :ϕ)
    @test !haskey(fd, :η_b)
  end

  @testset "variable_symbol dispatch" begin
    @test Entities.variable_symbol(fluid) == :ϕ
    @test Entities.variable_symbol(fsurf) == :κ
    @test Entities.variable_symbol(mem) == :η_m
    beam = Entities.EulerBernoulliBeam(L=20.0, mᵨ=922.5/ρw, EIᵨ=1e9*1e-4/ρw)
    @test Entities.variable_symbol(beam) == :η_b
  end

  # =========================================================================
  # Test single-variable Membrane2D weak forms (η_m only)
  # =========================================================================

  @testset "Membrane2D single-variable forms" begin
    l((w,u,v)) = ∫(0.0 * w)dΓin

    a_mass((ϕ,κ,η),(w,u,v)) = begin
      xd = WF.FieldMap((ϕ,κ,η), fmap)
      yd = WF.FieldMap((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + Entities.mass(mem, dom, xd, yd)
    end
    op_m = AffineFEOperator(a_mass, l, X, Y)
    @test nnz(get_matrix(op_m)) > 0

    a_stiff((ϕ,κ,η),(w,u,v)) = begin
      xd = WF.FieldMap((ϕ,κ,η), fmap)
      yd = WF.FieldMap((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + Entities.stiffness(mem, dom, xd, yd)
    end
    op_k = AffineFEOperator(a_stiff, l, X, Y)
    @test nnz(get_matrix(op_k)) > 0
  end

  # =========================================================================
  # Test PotentialFlow single-variable (Laplacian)
  # =========================================================================

  @testset "PotentialFlow single-variable stiffness" begin
    l((w,u,v)) = ∫(0.0 * w)dΓin

    a_fluid((ϕ,κ,η),(w,u,v)) = begin
      xd = WF.FieldMap((ϕ,κ,η), fmap)
      yd = WF.FieldMap((w,u,v), fmap)
      Entities.stiffness(fluid, dom, xd, yd)
    end
    op = AffineFEOperator(a_fluid, l, X, Y)
    @test nnz(get_matrix(op)) > 0
  end

  # =========================================================================
  # Test Fluid ↔ Structure coupling (damping only)
  # =========================================================================

  @testset "Fluid-Structure coupling damping" begin
    l((w,u,v)) = ∫(0.0 * w)dΓin

    a_coupling((ϕ,κ,η),(w,u,v)) = begin
      xd = WF.FieldMap((ϕ,κ,η), fmap)
      yd = WF.FieldMap((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + Entities.damping(fluid, mem, dom, xd, yd)
    end
    op = AffineFEOperator(a_coupling, l, X, Y)
    @test nnz(get_matrix(op)) > 0
  end

  # =========================================================================
  # Test composed weakform: fluid + structure + coupling
  # =========================================================================

  @testset "Composed weakform (fluid + free surface + membrane + coupling)" begin
    l((w,u,v)) = ∫(0.0 * w)dΓin

    a((ϕ,κ,η),(w,u,v)) = begin
      xd = WF.FieldMap((ϕ,κ,η), fmap)
      yd = WF.FieldMap((w,u,v), fmap)
      # single-variable
      Entities.weakform(fluid, dom, ω, xd, yd) +
      Entities.weakform(fsurf, dom, ω, xd, yd) +
      Entities.weakform(mem, dom, ω, xd, yd) +
      # coupling
      Entities.weakform(fluid, fsurf, dom, ω, xd, yd) +
      Entities.weakform(fluid, mem, dom, ω, xd, yd)
    end

    op = AffineFEOperator(a, l, X, Y)
    A  = get_matrix(op)
    @test size(A, 1) > 0
    @test nnz(A) > 0
  end

  # =========================================================================
  # Test assemble_weakform with fmap (single-variable terms)
  # =========================================================================

  @testset "assemble_weakform composition" begin
    terms = (fluid, fsurf, mem)

    a((ϕ,κ,η),(w,u,v)) =
      WF.assemble_weakform(terms, dom, ω, fmap, (ϕ,κ,η), (w,u,v)) +
      Entities.weakform(fluid, fsurf, dom, ω,
               WF.FieldMap((ϕ,κ,η), fmap), WF.FieldMap((w,u,v), fmap))

    l((w,u,v)) = ∫(0.0 * w)dΓin

    op = AffineFEOperator(a, l, X, Y)
    A  = get_matrix(op)
    @test size(A, 1) > 0
    @test nnz(A) > 0
  end

end
