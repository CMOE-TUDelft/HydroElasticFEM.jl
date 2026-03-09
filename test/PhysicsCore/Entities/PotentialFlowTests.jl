using Test
import HydroElasticFEM.PhysicsCore.Entities as E

@testset "PotentialFlow struct" begin
  pf = E.PotentialFlow()
  @test pf isa E.PhysicsParameters
  @test pf.ρw == 1025.0
  @test pf.g  == 9.81

  # Custom values
  pf2 = E.PotentialFlow(ρw=1000.0, g=9.80)
  @test pf2.ρw == 1000.0
  @test pf2.g  == 9.80

  # variable_symbol
  @test E.variable_symbol(pf) == :ϕ

  # trait queries
  @test E.has_mass_form(pf)    == false
  @test E.has_damping_form(pf) == false
  @test E.has_stiffness_form(pf) == true
  @test E.has_rhs_form(pf)      == true

  # print_parameters should not throw
  @test_nowarn E.print_parameters(pf)
end
