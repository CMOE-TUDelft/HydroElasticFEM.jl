using Gridap
using Gridap.Geometry
using Gridap.CellData
using SparseArrays

import HydroElasticFEM.PhysicsCore.FESpaceAssembly as FEA
import HydroElasticFEM.PhysicsCore.FESpaces as FES
import HydroElasticFEM.PhysicsCore.Entities as E
import HydroElasticFEM.PhysicsCore.Domains as D

# =========================================================================
# FESpaceAssembly tests
# =========================================================================

@testset "FESpaceAssembly" begin

  # -----------------------------------------------------------------------
  # FESpaceConfig construction
  # -----------------------------------------------------------------------

  @testset "FESpaceConfig defaults and overrides" begin
    fe = FES.FESpaceConfig()
    @test fe.reffe_type == lagrangian
    @test fe.space_type == Float64
    @test fe.order == 1
    @test fe.conformity == :H1
    @test fe.vector_type == Vector{ComplexF64}
    @test fe.γ == 10.0
    @test fe.dirichlet_tags === nothing
    @test fe.dirichlet_value === nothing

    fe2 = FES.FESpaceConfig(order=2)
    @test fe2.order == 2
    @test fe2.γ == 40.0   # 10.0 * 2^2

    fe3 = FES.FESpaceConfig(order=3, γ=100.0, dirichlet_tags="boundary")
    @test fe3.order == 3
    @test fe3.γ == 100.0
    @test fe3.dirichlet_tags == "boundary"
  end

  # -----------------------------------------------------------------------
  # Entity fe field backward compatibility
  # -----------------------------------------------------------------------

  @testset "Entity fe field defaults" begin
    fluid = E.PotentialFlow(ρw=1025.0, g=9.81)
    @test fluid.fe.order == 1

    fsurf = E.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5)
    @test fsurf.fe.order == 1

    mem = E.Membrane2D(L=20.0, mᵨ=1.0, Tᵨ=100.0)
    @test mem.fe.order == 1

    beam = E.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=100.0)
    @test beam.fe.order == 1
    @test beam.fe.γ == 10.0
  end

  @testset "Entity fe field overrides" begin
    beam = E.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=100.0,
                                fe=FES.FESpaceConfig(order=2))
    @test beam.fe.order == 2
    @test beam.fe.γ == 40.0
  end

  # -----------------------------------------------------------------------
  # build_fe_spaces: 3-field coupled problem
  # -----------------------------------------------------------------------

  @testset "build_fe_spaces — 3-field coupled" begin
    order = 1
    nx, ny = 10, 4
    domain = (0, 50, -10, 0)
    model = CartesianDiscreteModel(domain, (nx, ny))

    labels_Ω = get_face_labeling(model)
    add_tag_from_tags!(labels_Ω, "surface", [3, 4, 6])
    add_tag_from_tags!(labels_Ω, "bottom",  [1, 2, 5])

    Ω = Interior(model)
    Γ = Boundary(model, tags="surface")

    # Surface masking
    xm₀, xm₁ = 15.0, 35.0
    is_mem(xs) = let n=length(xs); x=(1/n)*sum(xs); (xm₀ <= x[1] <= xm₁) * (x[2] ≈ 0.0) end
    xΓ = get_cell_coordinates(Γ)
    Γm_mask = lazy_map(is_mem, xΓ)
    Γη  = Triangulation(Γ, findall(Γm_mask))
    Γκ  = Triangulation(Γ, findall(!, Γm_mask))

    fluid = E.PotentialFlow(ρw=1025.0, g=9.81, fe=FES.FESpaceConfig(order=order))
    fsurf = E.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5, fe=FES.FESpaceConfig(order=order))
    mem   = E.Membrane2D(L=20.0, mᵨ=0.9, Tᵨ=98.1, fe=FES.FESpaceConfig(order=order))

    X, Y, fmap = FEA.build_fe_spaces(
        fluid => Ω,
        fsurf => Γκ,
        mem   => Γη,
    )

    # Check fmap
    @test fmap[:ϕ] == 1
    @test fmap[:κ] == 2
    @test fmap[:η_m] == 3
    @test length(fmap) == 3

    # Check spaces are valid MultiFieldFESpaces
    @test length(X) == 3
    @test length(Y) == 3
  end

  # -----------------------------------------------------------------------
  # build_fe_spaces: single beam with Dirichlet BCs
  # -----------------------------------------------------------------------

  @testset "build_fe_spaces — beam with Dirichlet" begin
    model = CartesianDiscreteModel((0, 1.0), (20,))

    beam = E.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=100.0, g=0.0,
                                fe=FES.FESpaceConfig(order=2,
                                                     vector_type=Vector{Float64},
                                                     dirichlet_tags="boundary",
                                                     dirichlet_value=0.0))

    X, Y, fmap = FEA.build_fe_spaces(beam => model)

    @test fmap[:η_b] == 1
    @test length(fmap) == 1
    @test length(X) == 1
    @test length(Y) == 1
  end

  # -----------------------------------------------------------------------
  # build_fe_spaces integrates with weak forms
  # -----------------------------------------------------------------------

  @testset "build_fe_spaces — integration with weak forms" begin
    order = 1
    nx, ny = 10, 4
    domain = (0, 50, -10, 0)
    model = CartesianDiscreteModel(domain, (nx, ny))

    labels_Ω = get_face_labeling(model)
    add_tag_from_tags!(labels_Ω, "surface", [3, 4, 6])
    add_tag_from_tags!(labels_Ω, "bottom",  [1, 2, 5])
    add_tag_from_tags!(labels_Ω, "inlet",   [7])

    Ω = Interior(model)
    Γ = Boundary(model, tags="surface")
    Γin = Boundary(model, tags="inlet")

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

    fluid = E.PotentialFlow(ρw=1025.0, g=9.81, fe=FES.FESpaceConfig(order=order))
    fsurf = E.FreeSurface(ρw=1025.0, g=9.81, βₕ=0.5, fe=FES.FESpaceConfig(order=order))
    mem   = E.Membrane2D(L=20.0, mᵨ=0.9, Tᵨ=98.1, fe=FES.FESpaceConfig(order=order))

    X, Y, fmap = FEA.build_fe_spaces(
        fluid => Ω,
        fsurf => Γκ,
        mem   => Γη,
    )

    dΩ   = Measure(Ω, degree)
    dΓm  = Measure(Γm, degree)
    dΓfs = Measure(Γfs, degree)
    dΓin = Measure(Γin, degree)

    dom = D.WeakFormDomains(dΩ=dΩ, dΓ_fs=dΓfs, dΓ_s=dΓm, dΓ_in=dΓin)

    ω = 2.0

    a((ϕ,κ,η),(w,u,v)) = begin
      xd = D.FieldDict((ϕ,κ,η), fmap)
      yd = D.FieldDict((w,u,v), fmap)
      E.weakform(fluid, dom, ω, xd, yd) +
      E.weakform(fsurf, dom, ω, xd, yd) +
      E.weakform(mem, dom, ω, xd, yd) +
      E.weakform(fluid, fsurf, dom, ω, xd, yd) +
      E.weakform(fluid, mem, dom, ω, xd, yd)
    end

    l((w,u,v)) = ∫(0.0 * w)dΓin

    op = AffineFEOperator(a, l, X, Y)
    A  = get_matrix(op)
    @test size(A, 1) > 0
    @test nnz(A) > 0
  end

end
