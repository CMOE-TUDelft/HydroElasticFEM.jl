using Gridap
using Gridap.Geometry
using Gridap.CellData

# =========================================================================
# Mini mesh setup: 50m × 10m tank, membrane from x=15 to x=35
# =========================================================================

@testset "WeakFormAssembly" begin

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
  αₕ  = -im * ω / g_ * (1 - βₕ) / βₕ

  # Physics entities
  fluid = PotentialFlow(ρw=ρw, g=g_)
  mem   = Membrane2D(L=20.0, m=922.5, T=98.1*ρw, ρw=ρw, g=g_)

  # Field mapping
  fmap = Dict(:ϕ => 1, :κ => 2, :η_m => 3)

  dom = WeakFormDomains(dΩ=dΩ, dΓ_fs=dΓfs, dΓ_s=dΓm,
                        dΓ_in=dΓin, dΓ_ot=dΓot)

  # =========================================================================
  # Test WeakFormDomains + FieldDict construction
  # =========================================================================

  @testset "WeakFormDomains construction" begin
    d = WeakFormDomains(dΩ=dΩ, dΓ_s=dΓm)
    @test d[:dΩ] === dΩ
    @test haskey(d, :dΓ_s)
    @test !haskey(d, :dΓ_fs)
  end

  @testset "FieldDict construction" begin
    fd = FieldDict((1, 2, 3), fmap)
    @test fd[:ϕ] == 1
    @test fd[:κ] == 2
    @test fd[:η_m] == 3
    @test haskey(fd, :ϕ)
    @test !haskey(fd, :η_b)
  end

  @testset "variable_symbol dispatch" begin
    @test variable_symbol(fluid) == :ϕ
    @test variable_symbol(mem) == :η_m
    beam = EulerBernoulliBeam(L=20.0, m=922.5, E=1e9, I=1e-4)
    @test variable_symbol(beam) == :η_b
  end

  # =========================================================================
  # Test single-variable Membrane2D weak forms (η_m only)
  # =========================================================================

  @testset "Membrane2D single-variable forms" begin
    l((w,u,v)) = ∫(0.0 * w)dΓin

    a_mass((ϕ,κ,η),(w,u,v)) = begin
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + mass(mem, dom, xd, yd)
    end
    op_m = AffineFEOperator(a_mass, l, X, Y)
    @test nnz(get_matrix(op_m)) > 0

    a_stiff((ϕ,κ,η),(w,u,v)) = begin
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + stiffness(mem, dom, xd, yd)
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
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      stiffness(fluid, dom, xd, yd)
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
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + damping(fluid, mem, dom, xd, yd)
    end
    op = AffineFEOperator(a_coupling, l, X, Y)
    @test nnz(get_matrix(op)) > 0
  end

  # =========================================================================
  # Test composed weakform: fluid + structure + coupling
  # =========================================================================

  @testset "Composed weakform (fluid + membrane + coupling)" begin
    l((w,u,v)) = ∫(0.0 * w)dΓin

    a((ϕ,κ,η),(w,u,v)) = begin
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      # free surface BC (external, not entity-managed)
      ∫(βₕ * (u + αₕ * w) * (g_ * κ - im * ω * ϕ) + im * ω * w * κ)dΓfs +
      # single-variable
      weakform(fluid, dom, ω, xd, yd) +
      weakform(mem, dom, ω, xd, yd) +
      # coupling
      weakform(fluid, mem, dom, ω, xd, yd)
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
    terms = (fluid, mem)

    a((ϕ,κ,η),(w,u,v)) =
      ∫(βₕ * (u + αₕ * w) * (g_ * κ - im * ω * ϕ) + im * ω * w * κ)dΓfs +
      assemble_weakform(terms, dom, ω, fmap, (ϕ,κ,η), (w,u,v))

    l((w,u,v)) = ∫(0.0 * w)dΓin

    op = AffineFEOperator(a, l, X, Y)
    A  = get_matrix(op)
    @test size(A, 1) > 0
    @test nnz(A) > 0
  end

end
