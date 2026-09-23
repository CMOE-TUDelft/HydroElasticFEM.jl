using Test
import HydroElasticFEM.Physics as P

@testset "TensionedEulerBernoulliBeam struct" begin
  ρw = 1025.0
  EIᵨ = 500e6 * 6.667e-4 / ρw
  Tᵨ  = 98.1

  beam = P.TensionedEulerBernoulliBeam(
    L=20.0, mᵨ=192.956/ρw, EIᵨ=EIᵨ, Tᵨ=Tᵨ, τ=0.0)
  @test beam.EIᵨ ≈ EIᵨ
  @test beam.Tᵨ == Tᵨ
  @test beam.τ == 0.0
  @test beam.ωn1 ≈ sqrt(
    EIᵨ / (192.956/ρw) * (π/20.0)^4 + Tᵨ / (192.956/ρw) * (π/20.0)^2)
  @test beam isa P.AbstractHydroelasticStructure
  @test beam isa P.Structure

  # Defaults: τ, g, joints default; derived ωn1 auto-computed
  beam_def = P.TensionedEulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=EIᵨ, Tᵨ=Tᵨ)
  @test beam_def.τ == 0.0
  @test beam_def.g == 9.81
  @test isempty(beam_def.joints)
  @test beam_def.ωn1 ≈ sqrt(
    EIᵨ / (192.956/ρw) * (π/20.0)^4 + Tᵨ / (192.956/ρw) * (π/20.0)^2)

  # Nonzero τ
  beam_d = P.TensionedEulerBernoulliBeam(
    L=10.0, mᵨ=100.0/ρw, EIᵨ=1e9*1e-3/ρw, Tᵨ=50.0, τ=0.05)
  @test beam_d.τ == 0.05

  # EIᵨ as a Function ⟹ ωn1 is not analytically available
  beam_var = P.TensionedEulerBernoulliBeam(
    L=20.0, mᵨ=192.956/ρw, EIᵨ=(x -> EIᵨ * (1 + 0.1x[1])), Tᵨ=Tᵨ)
  @test beam_var.ωn1 === nothing
end

@testset "TensionedEulerBernoulliBeam limiting cases (struct level)" begin
  ρw = 1025.0
  mᵨ  = 192.956 / ρw
  EIᵨ = 500e6 * 6.667e-4 / ρw
  Tᵨ  = 98.1
  L   = 20.0

  # EIᵨ = 0 ⟹ ωn1 matches Membrane's own ωn1 formula exactly.
  beam_no_EI = P.TensionedEulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=0.0, Tᵨ=Tᵨ)
  membrane   = P.Membrane(L=L, mᵨ=mᵨ, Tᵨ=Tᵨ)
  @test beam_no_EI.ωn1 ≈ membrane.ωn1

  # Tᵨ = 0 ⟹ mass_density/damping_parameter/mᵨ/EIᵨ match a plain EulerBernoulliBeam.
  beam_no_T = P.TensionedEulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=EIᵨ, Tᵨ=0.0, τ=0.03)
  euler     = P.EulerBernoulliBeam(L=L, mᵨ=mᵨ, EIᵨ=EIᵨ, τ=0.03)
  @test P.mass_density(beam_no_T) == P.mass_density(euler)
  @test P.damping_parameter(beam_no_T) == P.damping_parameter(euler)
end

@testset "TensionedEulerBernoulliBeam joints" begin
  ρw = 1025.0
  EIᵨ = 500e6 * 6.667e-4 / ρw
  Tᵨ  = 98.1

  # Reuses the same JointRotationalSpring type as EulerBernoulliBeam.
  joint = P.JointRotationalSpring(:dΛj_1, :n_Λ_j_1, 1.25e4)
  beam_one = P.TensionedEulerBernoulliBeam(
    L=20.0, mᵨ=192.956 / ρw, EIᵨ=EIᵨ, Tᵨ=Tᵨ,
    joints=[joint],
  )
  @test length(beam_one.joints) == 1
  @test beam_one.joints[1] == joint

  # Backward-compatible default remains no joints.
  beam_default = P.TensionedEulerBernoulliBeam(L=20.0, mᵨ=192.956 / ρw, EIᵨ=EIᵨ, Tᵨ=Tᵨ)
  @test isempty(beam_default.joints)
end

@testset "AbstractHydroelasticStructure interface" begin
  ρw = 1025.0
  beam = P.TensionedEulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=1.0, Tᵨ=1.0, τ=0.02)

  @test P.mass_density(beam) == beam.mᵨ
  @test P.damping_parameter(beam) == beam.τ
  @test P.gravitational_acceleration(beam) == beam.g

  # Membrane and EulerBernoulliBeam are also AbstractHydroelasticStructure now.
  membrane = P.Membrane(L=20.0, mᵨ=192.956/ρw, Tᵨ=98.1)
  euler    = P.EulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=1.0)
  @test membrane isa P.AbstractHydroelasticStructure
  @test euler isa P.AbstractHydroelasticStructure
  @test P.mass_density(membrane) == membrane.mᵨ
  @test P.mass_density(euler) == euler.mᵨ
end
