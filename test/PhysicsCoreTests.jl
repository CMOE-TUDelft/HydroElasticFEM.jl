import HydroElasticFEM.Physics as Physics
using Test

@testset "Physics Core Tests" begin
    # Test abstract type existence
    @test Physics.PhysicsParameters isa Type

    # Test that default print_parameters throws error
    struct TestParams <: Physics.PhysicsParameters end
    params = TestParams()
    
    @test_throws ErrorException Physics.print_parameters(params)
    
    # Test checking the error message specifically
    try
        Physics.print_parameters(params)
    catch e
        @test e isa ErrorException
        @test startswith(e.msg, "print_parameters not implemented for")
    end
end
