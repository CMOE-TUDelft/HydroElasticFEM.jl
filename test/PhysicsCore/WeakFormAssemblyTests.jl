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
  k   = ω^2 / g_
  βₕ  = 0.5
  αₕ  = -im * ω / g_ * (1 - βₕ) / βₕ

  # Field mapping: symbol -> positional index in the multi-field tuple
  fmap = Dict(:ϕ => 1, :κ => 2, :η_m => 3)

  # =========================================================================
  # Test WeakFormDomains construction (Dict-based)
  # =========================================================================

  @testset "WeakFormDomains construction" begin
    dom = WeakFormDomains(dΩ=dΩ, dΓ_s=dΓm)
    @test dom[:dΩ] === dΩ
    @test dom[:dΓ_s] === dΓm
    @test haskey(dom, :dΩ)
    @test !haskey(dom, :dΓ_fs)

    dom2 = WeakFormDomains(dΩ=dΩ, dΓ_fs=dΓfs, dΓ_s=dΓm,
                           dΓ_in=dΓin, dΓ_ot=dΓot)
    @test dom2[:dΓ_in] === dΓin
    @test dom2[:dΓ_ot] === dΓot
    @test haskey(dom2, :dΓ_fs)
  end

  # =========================================================================
  # Test FieldDict construction
  # =========================================================================

  @testset "FieldDict construction" begin
    fd = FieldDict((1, 2, 3), fmap)
    @test fd[:ϕ] == 1
    @test fd[:κ] == 2
    @test fd[:η_m] == 3
    @test haskey(fd, :ϕ)
    @test !haskey(fd, :η_b)
  end

  dom = WeakFormDomains(dΩ=dΩ, dΓ_fs=dΓfs, dΓ_s=dΓm,
                        dΓ_in=dΓin, dΓ_ot=dΓot)

  # =========================================================================
  # Test Membrane2D linear weak forms via FieldDict
  # =========================================================================

  @testset "Membrane2D linear forms (FieldDict)" begin
    mem = Membrane2D(L=20.0, m=922.5, T=98.1*1025.0, ρw=ρw, g=g_,
                     bndType=FreeBoundary())

    l((w,u,v)) = ∫(0.0 * w)dΓin

    # mass form via FieldDict
    a_mass((ϕ,κ,η),(w,u,v)) = begin
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + mass(mem, dom, xd, yd)
    end
    op_m = AffineFEOperator(a_mass, l, X, Y)
    @test nnz(get_matrix(op_m)) > 0

    # damping form via FieldDict
    a_damp((ϕ,κ,η),(w,u,v)) = begin
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + damping(mem, dom, xd, yd)
    end
    op_c = AffineFEOperator(a_damp, l, X, Y)
    @test nnz(get_matrix(op_c)) > 0

    # stiffness form via FieldDict
    a_stiff((ϕ,κ,η),(w,u,v)) = begin
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + stiffness(mem, dom, xd, yd)
    end
    op_k = AffineFEOperator(a_stiff, l, X, Y)
    @test nnz(get_matrix(op_k)) > 0
  end

  # =========================================================================
  # Test Membrane2D weakform (derived from linear forms via ω)
  # =========================================================================

  @testset "Membrane2D weakform (derived)" begin
    mem_free = Membrane2D(L=20.0, m=922.5, T=98.1*1025.0, ρw=ρw, g=g_,
                          bndType=FreeBoundary())

    l((w,u,v)) = ∫(0.0 * w)dΓin

    a_wf((ϕ,κ,η),(w,u,v)) = begin
      xd = FieldDict((ϕ,κ,η), fmap)
      yd = FieldDict((w,u,v), fmap)
      ∫(∇(w) ⋅ ∇(ϕ))dΩ +
      ∫(βₕ * (u + αₕ * w) * (g_ * κ - im * ω * ϕ) + im * ω * w * κ)dΓfs +
      weakform(mem_free, dom, ω, xd, yd)
    end

    op = AffineFEOperator(a_wf, l, X, Y)
    A  = get_matrix(op)
    @test nnz(A) > 0
  end

  # =========================================================================
  # Test assemble_weakform (generic loop with fmap wrapping)
  # =========================================================================

  @testset "assemble_weakform composition" begin
    mem = Membrane2D(L=20.0, m=922.5, T=98.1*1025.0, ρw=ρw, g=g_)
    terms = (mem,)

    a((ϕ,κ,η),(w,u,v)) =
      ∫(∇(w) ⋅ ∇(ϕ))dΩ +
      ∫(βₕ * (u + αₕ * w) * (g_ * κ - im * ω * ϕ) + im * ω * w * κ)dΓfs +
      assemble_weakform(terms, dom, ω, fmap, (ϕ,κ,η), (w,u,v))

    l((w,u,v)) = ∫(0.0 * w)dΓin

    op = AffineFEOperator(a, l, X, Y)
    A  = get_matrix(op)
    @test size(A, 1) > 0
    @test nnz(A) > 0
  end

  # =========================================================================
  # Test assemble_mass / assemble_damping / assemble_stiffness
  # =========================================================================

  @testset "assemble linear forms" begin
    mem = Membrane2D(L=20.0, m=922.5, T=98.1*1025.0, ρw=ρw, g=g_)
    terms = (mem,)

    l((w,u,v)) = ∫(0.0 * w)dΓin

    a_m((ϕ,κ,η),(w,u,v)) =
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + assemble_mass(terms, dom, fmap, (ϕ,κ,η), (w,u,v))
    op_m = AffineFEOperator(a_m, l, X, Y)
    @test nnz(get_matrix(op_m)) > 0

    a_c((ϕ,κ,η),(w,u,v)) =
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + assemble_damping(terms, dom, fmap, (ϕ,κ,η), (w,u,v))
    op_c = AffineFEOperator(a_c, l, X, Y)
    @test nnz(get_matrix(op_c)) > 0

    a_k((ϕ,κ,η),(w,u,v)) =
      ∫(∇(w) ⋅ ∇(ϕ))dΩ + assemble_stiffness(terms, dom, fmap, (ϕ,κ,η), (w,u,v))
    op_k = AffineFEOperator(a_k, l, X, Y)
    @test nnz(get_matrix(op_k)) > 0
  end

end
