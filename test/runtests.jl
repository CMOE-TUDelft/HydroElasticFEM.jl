using HydroElasticFEM
using Test

@testset "HydroElasticFEM.jl" begin
    include("PhysicsCoreTests.jl")
    include("EulerBernoulliBeamTests.jl")
end
