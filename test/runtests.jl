using Test
using HydroElasticFEM
using Gridap

@testset "HydroElasticFEM" begin

  @testset "Membrane construction" begin
    mem = HydroElasticFEM.Membrane.Membrane2D(
      20.0, 922.5, 98.1 * 1025.0, 0.0,
      HydroElasticFEM.Membrane.Free())
    @test mem.L == 20.0
    @test mem.m == 922.5
    @test mem.T == 98.1 * 1025.0
    @test mem.τ == 0.0
    @test mem.MTotal ≈ 922.5 * 20.0
    @test mem.ωn1 > 0
    @test mem.ωn1 ≈ (π / 20.0) * sqrt(98.1 * 1025.0 / 922.5)
  end

  @testset "Beam construction" begin
    beam = HydroElasticFEM.BeamNoJoints.Beam2D(
      20.0, 192.956, 500e6, 6.667e-4, 0.0,
      HydroElasticFEM.BeamNoJoints.Free())
    @test beam.L == 20.0
    @test beam.EI ≈ 500e6 * 6.667e-4
    @test beam.τEI ≈ 0.0
    @test beam.MTotal ≈ 192.956 * 20.0
    @test beam.ωn1 ≈ 22.3733 * sqrt(beam.EI / (192.956 * 20.0^4))
  end

  @testset "Resonator construction" begin
    resn = HydroElasticFEM.Resonator.Single(
      1e3, 5.9e3, 0.0, VectorValue(10.0, 0.0))
    @test resn.M == 1e3
    @test resn.K == 5.9e3
    @test resn.C == 0.0
    @test resn.ωn1 ≈ sqrt(5.9e3 / 1e3)
  end

  @testset "Resonator Array1D" begin
    xz = [VectorValue(5.0, 0.0), VectorValue(10.0, 0.0), VectorValue(15.0, 0.0)]
    arr = HydroElasticFEM.Resonator.Array1D(3, 100.0, 500.0, 10.0, xz)
    @test length(arr) == 3
    @test arr[1].M == 100.0
    @test arr[2].XZ == VectorValue(10.0, 0.0)
    @test arr[3].ωn1 ≈ sqrt(500.0 / 100.0)
  end

  @testset "WaveInput AiryWaveXZ" begin
    WI = HydroElasticFEM.WaveInput_FrequencyDomain
    wave = WI.AiryWaveXZ(10.0, 2.0, 0.1)
    @test wave.k > 0
    @test wave.kh > 0
    @test wave.h == 10.0
    @test wave.ω == 2.0
    @test wave.η0 == 0.1
    @test wave.λ ≈ 2π / wave.k

    eta = WI.surface_elevation(wave, VectorValue(0.0, 0.0))
    @test abs(eta) ≈ 0.1 atol = 1e-10
  end

  # ==========================================================================
  # Phase 1: PhysicsCore/PhysicalEntities.jl -- module-level tests
  # ==========================================================================

  @testset "PhysicalEntities - module interface" begin
    PE = HydroElasticFEM.PhysicalEntities

    # Abstract base type exists
    @test PE.PhysicsParameters isa Type

    # Default print_parameters throws for unknown subtypes
    struct _TestParams <: PE.PhysicsParameters end
    @test_throws ErrorException PE.print_parameters(_TestParams())

    # print_parameters works for concrete types
    mem = PE.Membrane2D(20.0, 922.5, 98.1 * 1025.0, 0.0, PE.FreeBoundary())
    @test_nowarn PE.print_parameters(mem)
    beam = PE.Beam2D(20.0, 192.956, 500e6, 6.667e-4, 0.0, PE.FreeBoundary())
    @test_nowarn PE.print_parameters(beam)
    resn = PE.ResonatorSingle(1e3, 5.9e3, 0.0, VectorValue(10.0, 0.0))
    @test_nowarn PE.print_parameters(resn)
  end

  # ==========================================================================
  # Phase 1: PhysicsCore/PhysicalEntities.jl -- via top-level exports
  # ==========================================================================

  @testset "PhysicalEntities - Membrane2D" begin
    mem = HydroElasticFEM.Membrane2D(
      20.0, 922.5, 98.1 * 1025.0, 0.0, FreeBoundary())
    @test mem.L == 20.0
    @test mem.m == 922.5
    @test mem.T == 98.1 * 1025.0
    @test mem.τ == 0.0
    @test mem.MTotal ≈ 922.5 * 20.0
    @test mem.ωn1 ≈ (π / 20.0) * sqrt(98.1 * 1025.0 / 922.5)
    @test mem isa HydroElasticFEM.AbstractStructure

    # Fixed boundary
    mem_fix = HydroElasticFEM.Membrane2D(
      10.0, 500.0, 1000.0, 0.1, FixedBoundary())
    @test mem_fix.bndType isa FixedBoundary
  end

  @testset "PhysicalEntities - Beam2D" begin
    beam = HydroElasticFEM.Beam2D(
      20.0, 192.956, 500e6, 6.667e-4, 0.0, FreeBoundary())
    @test beam.EI ≈ 500e6 * 6.667e-4
    @test beam.τEI ≈ 0.0
    @test beam.ωn1 ≈ 22.3733 * sqrt(beam.EI / (192.956 * 20.0^4))
    @test beam isa HydroElasticFEM.AbstractStructure
  end

  @testset "PhysicalEntities - ResonatorSingle" begin
    resn = ResonatorSingle(1e3, 5.9e3, 0.0, VectorValue(10.0, 0.0))
    @test resn.M == 1e3
    @test resn.ωn1 ≈ sqrt(5.9e3 / 1e3)
  end

  @testset "PhysicalEntities - resonator_array" begin
    xz = [VectorValue(5.0, 0.0), VectorValue(10.0, 0.0), VectorValue(15.0, 0.0)]
    arr = resonator_array(3, 100.0, 500.0, 10.0, xz)
    @test length(arr) == 3
    @test arr[1] isa ResonatorSingle
    @test arr[2].XZ == VectorValue(10.0, 0.0)
    @test arr[3].ωn1 ≈ sqrt(500.0 / 100.0)

    # Vector variant
    arr2 = resonator_array(2, [100.0, 200.0], [500.0, 600.0], [10.0, 20.0],
      [VectorValue(5.0, 0.0), VectorValue(15.0, 0.0)])
    @test arr2[1].M == 100.0
    @test arr2[2].M == 200.0
  end

  @testset "PhysicalEntities - backward compat" begin
    # Old API still produces same physics
    mem_old = HydroElasticFEM.Membrane.Membrane2D(
      20.0, 922.5, 98.1 * 1025.0, 0.0,
      HydroElasticFEM.Membrane.Free())
    mem_new = HydroElasticFEM.Membrane2D(
      20.0, 922.5, 98.1 * 1025.0, 0.0, FreeBoundary())
    @test mem_old.ωn1 ≈ mem_new.ωn1
    @test mem_old.MTotal ≈ mem_new.MTotal

    beam_old = HydroElasticFEM.BeamNoJoints.Beam2D(
      20.0, 192.956, 500e6, 6.667e-4, 0.0,
      HydroElasticFEM.BeamNoJoints.Free())
    beam_new = HydroElasticFEM.Beam2D(
      20.0, 192.956, 500e6, 6.667e-4, 0.0, FreeBoundary())
    @test beam_old.ωn1 ≈ beam_new.ωn1
    @test beam_old.EI ≈ beam_new.EI
  end

  # ==========================================================================
  # Utilities
  # ==========================================================================

  @testset "Utilities" begin
    # Uniform mesh (r ≈ 1)
    y_uniform = HydroElasticFEM.map_vertical_GP_for_const_dep(-5.0, 1.0, 10, 10.0)
    @test y_uniform ≈ -5.0

    # GP mesh (r > 1)
    y_gp = HydroElasticFEM.map_vertical_GP_for_const_dep(-5.0, 1.2, 10, 10.0)
    @test y_gp < 0
    @test y_gp >= -10.0

    # Boundary values
    y_zero = HydroElasticFEM.map_vertical_GP_for_const_dep(0.0, 1.2, 10, 10.0)
    @test y_zero ≈ 0.0
  end

  # ==========================================================================
  # Time-domain solvers
  # ==========================================================================

  include("test_mem_time_lrmm.jl")

end
