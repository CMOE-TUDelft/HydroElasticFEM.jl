@testset "ResonatorSingle" begin
  resn = ResonatorSingle(M=1e3, K=5.9e3, C=0.0, XZ=VectorValue(10.0, 0.0))
  @test resn.M == 1e3
  @test resn.ωn1 ≈ sqrt(5.9e3 / 1e3)

  # Defaults: C and XZ default, derived ωn1
  resn_def = ResonatorSingle(M=1e3, K=5.9e3)
  @test resn_def.C == 0.0
  @test resn_def.XZ == VectorValue(0.0, 0.0)
  @test resn_def.ωn1 ≈ sqrt(5.9e3 / 1e3)
end

@testset "resonator_array" begin
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
