import HydroElasticFEM.Physics as Physics
using Test

@testset "EulerBernoulliBeam Tests" begin
    # Test default values
    beam_params = Physics.EBBeamParameters()
    @test beam_params isa Physics.EBBeamParameters
    @test beam_params isa Physics.PhysicsParameters
    @test beam_params.E == 210e9
    @test beam_params.I == 1e-6

    # Test custom values
    custom_E = 70e9
    custom_I = 2e-5
    custom_beam = Physics.EBBeamParameters(E=custom_E, I=custom_I)
    @test custom_beam.E == custom_E
    @test custom_beam.I == custom_I

    # Test kwarg constructor
    kwarg_beam = Physics.EBBeamParameters(I=5e-6)
    @test kwarg_beam.E == 210e9
    @test kwarg_beam.I == 5e-6
    
    # Test print_parameters
    # Should not throw error
    @test_nowarn Physics.print_parameters(beam_params)
end
