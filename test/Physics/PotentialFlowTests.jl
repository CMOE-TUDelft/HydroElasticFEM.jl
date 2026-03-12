using Test
import HydroElasticFEM.Physics as P

@testset "PotentialFlow struct" begin
  pf = P.PotentialFlow()
  @test pf isa P.PhysicsParameters
  @test pf.ρw == 1025.0
  @test pf.g  == 9.81

  # Custom values
  pf2 = P.PotentialFlow(ρw=1000.0, g=9.80)
  @test pf2.ρw == 1000.0
  @test pf2.g  == 9.80

  # variable_symbol
  @test P.variable_symbol(pf) == :ϕ

  # trait queries
  @test P.has_mass_form(pf)    == false
  @test P.has_damping_form(pf) == false
  @test P.has_stiffness_form(pf) == true
  @test P.has_rhs_form(pf)      == true

  # print_parameters should not throw
  @test_nowarn P.print_parameters(pf)
end
