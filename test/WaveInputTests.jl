using WaveSpec

@testset "WaveInput AiryWaveXZ" begin
  WI = HydroElasticFEM.WaveInput_FrequencyDomain
  wave = WI.AiryWaveXZ(10.0, 2.0, 0.1)

  @test wave.k ≈ WaveSpec.AiryWaves.solve_wavenumber(2.0, 10.0)
  @test wave.k > 0
  @test wave.kh > 0
  @test wave.h == 10.0
  @test wave.ω == 2.0
  @test wave.η0 == 0.1
  @test wave.λ ≈ 2π / wave.k

  eta = WI.surface_elevation(wave, VectorValue(0.0, 0.0))
  @test abs(eta) ≈ 0.1 atol = 1e-10
end
