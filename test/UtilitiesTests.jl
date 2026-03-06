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
