using Test
import HydroElasticFEM.Physics as P

@testset "EulerBernoulliBeam struct" begin
  ρw = 1025.0
  EIᵨ = 500e6 * 6.667e-4 / ρw
  beam = P.EulerBernoulliBeam(
    L=20.0, mᵨ=192.956/ρw, EIᵨ=EIᵨ, τ=0.0)
  @test beam.EIᵨ ≈ EIᵨ
  @test beam.ωn1 ≈ 22.3733 * sqrt(EIᵨ / ((192.956/ρw) * 20.0^4))
  @test beam isa P.Structure

  # Defaults: τ and bndType default, derived fields auto-computed
  beam_def = P.EulerBernoulliBeam(L=20.0, mᵨ=192.956/ρw, EIᵨ=EIᵨ)
  @test beam_def.τ == 0.0
  @test beam_def.EIᵨ ≈ EIᵨ
  @test beam_def.ωn1 ≈ 22.3733 * sqrt(beam_def.EIᵨ / ((192.956/ρw) * 20.0^4))

  # Nonzero τ
  beam_d = P.EulerBernoulliBeam(L=10.0, mᵨ=100.0/ρw, EIᵨ=1e9*1e-3/ρw, τ=0.05)
  @test beam_d.τ == 0.05
end

@testset "EulerBernoulliBeam joints" begin
  ρw = 1025.0
  EIᵨ = 500e6 * 6.667e-4 / ρw

  # Explicit rotational spring joint definition
  joint = P.JointRotationalSpring(:dΛj_1, :n_Λ_j_1, 1.25e4)
  @test joint.domain_symbol == :dΛj_1
  @test joint.normal_symbol == :n_Λ_j_1
  @test joint.kᵣ == 1.25e4

  # Beam accepts and stores one or many joints
  beam_one = P.EulerBernoulliBeam(
    L=20.0, mᵨ=192.956 / ρw, EIᵨ=EIᵨ,
    joints=[joint],
  )
  @test length(beam_one.joints) == 1
  @test beam_one.joints[1] == joint

  joint_2 = P.JointRotationalSpring(:dΛj_2, :n_Λ_j_2, 2.50e4)
  beam_two = P.EulerBernoulliBeam(
    L=20.0, mᵨ=192.956 / ρw, EIᵨ=EIᵨ,
    joints=[joint, joint_2],
  )
  @test length(beam_two.joints) == 2
  @test beam_two.joints[2].domain_symbol == :dΛj_2
  @test beam_two.joints[2].normal_symbol == :n_Λ_j_2
  @test beam_two.joints[2].kᵣ == 2.50e4

  # Backward-compatible default remains no joints
  beam_default = P.EulerBernoulliBeam(L=20.0, mᵨ=192.956 / ρw, EIᵨ=EIᵨ)
  @test isempty(beam_default.joints)
end
