using Test
import HydroElasticFEM.Physics as P

@testset "FreeSurface struct" begin
  fs = P.FreeSurface()
  @test fs isa P.PhysicsParameters
  @test fs.ρw == 1025.0
  @test fs.g  == 9.81
  @test fs.βₕ == 0.5

  # Custom values
  fs2 = P.FreeSurface(ρw=1000.0, g=9.80, βₕ=0.3)
  @test fs2.ρw == 1000.0
  @test fs2.g  == 9.80
  @test fs2.βₕ == 0.3

  # variable_symbol
  @test P.variable_symbol(fs) == :κ

  # print_parameters should not throw
  @test_nowarn P.print_parameters(fs)
end
