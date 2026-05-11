using Test
import HydroElasticFEM.Physics as P

@testset "Membrane struct" begin
  ρw = 1025.0
  mem = P.Membrane(
    L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1, τ=0.0)
  @test mem.L == 20.0
  @test mem.mᵨ ≈ 922.5 / ρw
  @test mem.Tᵨ == 98.1
  @test mem.τ == 0.0
  @test mem.ωn1 ≈ (π / 20.0) * sqrt(98.1 / (922.5 / ρw))
  @test mem isa P.Structure

  # Defaults: τ and bndType default, derived fields auto-computed
  mem_def = P.Membrane(L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1)
  @test mem_def.τ == 0.0
  @test mem_def.ωn1 ≈ (π / 20.0) * sqrt(98.1 / (922.5 / ρw))

  # Backward-compatible alias
  mem_alias = P.Membrane(L=20.0, mᵨ=922.5/ρw, Tᵨ=98.1)
  @test mem_alias isa P.Membrane
end
