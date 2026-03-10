using Test
using Gridap
using HydroElasticFEM

@testset "Legacy Membrane construction" begin
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

@testset "Legacy Beam construction" begin
  beam = HydroElasticFEM.BeamNoJoints.Beam2D(
    20.0, 192.956, 500e6, 6.667e-4, 0.0,
    HydroElasticFEM.BeamNoJoints.Free())
  @test beam.L == 20.0
  @test beam.EI ≈ 500e6 * 6.667e-4
  @test beam.τEI ≈ 0.0
  @test beam.MTotal ≈ 192.956 * 20.0
  @test beam.ωn1 ≈ 22.3733 * sqrt(beam.EI / (192.956 * 20.0^4))
end

@testset "Legacy Resonator construction" begin
  resn = HydroElasticFEM.Resonator.Single(
    1e3, 5.9e3, 0.0, VectorValue(10.0, 0.0))
  @test resn.M == 1e3
  @test resn.K == 5.9e3
  @test resn.C == 0.0
  @test resn.ωn1 ≈ sqrt(5.9e3 / 1e3)
end

@testset "Legacy Resonator Array1D" begin
  xz = [VectorValue(5.0, 0.0), VectorValue(10.0, 0.0), VectorValue(15.0, 0.0)]
  arr = HydroElasticFEM.Resonator.Array1D(3, 100.0, 500.0, 10.0, xz)
  @test length(arr) == 3
  @test arr[1].M == 100.0
  @test arr[2].XZ == VectorValue(10.0, 0.0)
  @test arr[3].ωn1 ≈ sqrt(500.0 / 100.0)
end

@testset "Backward compat - old vs new API" begin
  ρw = 1025.0
  # Old API still produces same physics
  mem_old = HydroElasticFEM.Membrane.Membrane2D(
    20.0, 922.5, 98.1 * ρw, 0.0,
    HydroElasticFEM.Membrane.Free())
  mem_new = HydroElasticFEM.PhysicsCore.Entities.Membrane2D(
    L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1, τ=0.0)
  @test mem_old.ωn1 ≈ mem_new.ωn1

  beam_old = HydroElasticFEM.BeamNoJoints.Beam2D(
    20.0, 192.956, 500e6, 6.667e-4, 0.0,
    HydroElasticFEM.BeamNoJoints.Free())
  beam_new = HydroElasticFEM.PhysicsCore.Entities.EulerBernoulliBeam(
    L=20.0, mᵨ=192.956/ρw, EIᵨ=500e6*6.667e-4/ρw, τ=0.0)
  @test beam_old.ωn1 ≈ beam_new.ωn1
end
