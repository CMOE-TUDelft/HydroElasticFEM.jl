using Test

# =========================================================================
# Integration tests for the Khabakhpasheva beam-joint frequency-domain
# example (two stiffness cases: ξ = 0 and ξ = 625).
#
# Reference: Khabakhpasheva (2002) — flexible plate on the free surface
# of a fluid with a hinge joint.
#
# Both cases use the same wave conditions (α·Lb incidence).
# ξ = 0     → no rotational spring at the joint  (kᵣ = 0)
# ξ = 625   → finite rotational spring stiffness at the joint
#
# Tests run at minimal resolution (nx=2, ny=1, order=2) to remain fast.
# They check structural properties of the solution rather than
# high-accuracy values.
# =========================================================================

include(joinpath(@__DIR__, "..", "..", "examples", "KhabakhpashevaBeamJointExample.jl"))
using .KhabakhpashevaBeamJointExample

@testset "Khabakhpasheva beam-joint (two stiffness cases)" begin

    nx_test  = 2
    ny_test  = 1
    ord_test = 2

    # -----------------------------------------------------------------------
    # Case ξ = 0: no rotational spring
    # -----------------------------------------------------------------------
    @testset "ξ = 0 (no spring)" begin
        p0 = KhabakhpashevaCaseParams(
            name       = "test_xi_0",
            nx         = nx_test,
            ny         = ny_test,
            order      = ord_test,
            ξ          = 0.0,
            vtk_output = false,
            make_plot  = false,
        )

        xs0, eta0, meta0 = run_khabakhpasheva_case(p0)

        # spring constant must be zero
        @test meta0.kᵣ == 0.0

        # output arrays have at least 400 probe points (the example appends a
        # near-joint cluster on top of the 400-point uniform grid)
        @test length(xs0) >= 400
        @test length(eta0) >= 400
        @test length(xs0) == length(eta0)

        # normalised coordinates span [0, 1]
        @test xs0[1]   ≈ 0.0  atol=1e-12
        @test xs0[end] ≈ 1.0  atol=1e-12

        # |η|/η₀ must be real and non-negative everywhere
        @test all(eta0 .>= 0.0)

        # physical sanity: maximum response should be of order 1
        # (neither negligibly small nor blowing up)
        @test maximum(eta0) > 0.1
        @test maximum(eta0) < 20.0
    end

    # -----------------------------------------------------------------------
    # Case ξ = 625: finite rotational spring at the joint
    # -----------------------------------------------------------------------
    @testset "ξ = 625 (finite spring)" begin
        p625 = KhabakhpashevaCaseParams(
            name       = "test_xi_625",
            nx         = nx_test,
            ny         = ny_test,
            order      = ord_test,
            ξ          = 625.0,
            vtk_output = false,
            make_plot  = false,
        )

        xs625, eta625, meta625 = run_khabakhpasheva_case(p625)

        # spring constant must be strictly positive
        @test meta625.kᵣ > 0.0

        # same array-length invariants
        @test length(xs625) >= 400
        @test length(eta625) >= 400
        @test length(xs625) == length(eta625)

        # coordinates span [0, 1]
        @test xs625[1]   ≈ 0.0  atol=1e-12
        @test xs625[end] ≈ 1.0  atol=1e-12

        # |η|/η₀ must be real and non-negative
        @test all(eta625 .>= 0.0)

        # physical sanity bounds
        @test maximum(eta625) > 0.1
        @test maximum(eta625) < 20.0
    end

    # -----------------------------------------------------------------------
    # Consistency: both cases share the same wave parameters (ω, k, EIᵨ)
    # -----------------------------------------------------------------------
    @testset "Both cases share the same wave parameters" begin
        p0 = KhabakhpashevaCaseParams(
            name="test_shared_ω_ξ0", nx=nx_test, ny=ny_test, order=ord_test,
            ξ=0.0, vtk_output=false, make_plot=false)
        p625 = KhabakhpashevaCaseParams(
            name="test_shared_ω_ξ625", nx=nx_test, ny=ny_test, order=ord_test,
            ξ=625.0, vtk_output=false, make_plot=false)

        _, _, m0   = run_khabakhpasheva_case(p0)
        _, _, m625 = run_khabakhpasheva_case(p625)

        @test m0.ω   ≈ m625.ω   rtol=1e-12
        @test m0.k   ≈ m625.k   rtol=1e-12
        @test m0.EIᵨ ≈ m625.EIᵨ rtol=1e-12

        # spring stiffness for ξ=625 must exceed that for ξ=0
        @test m625.kᵣ > m0.kᵣ
    end

    # -----------------------------------------------------------------------
    # run_khabakhpasheva_two_cases convenience wrapper
    # -----------------------------------------------------------------------
    @testset "run_khabakhpasheva_two_cases wrapper" begin
        results = run_khabakhpasheva_two_cases(
            nx = nx_test, ny = ny_test, order = ord_test,
            vtk_output = false, make_plot = false)

        @test haskey(results, :with_joint)
        @test haskey(results, :without_joint)

        @test length(results.with_joint.xs)    >= 400
        @test length(results.without_joint.xs) >= 400

        # with_joint aliases the ξ=625 (stiff spring) case, so kᵣ must be positive
        @test results.with_joint.meta.kᵣ > 0.0

        # without_joint aliases the ξ=0 (hinged) case, so kᵣ must be zero
        @test results.without_joint.meta.kᵣ == 0.0

        # Both produce non-negative deflection fields
        @test all(results.with_joint.η_rel    .>= 0.0)
        @test all(results.without_joint.η_rel .>= 0.0)
    end

end
