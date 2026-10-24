using Test

@testset "PostProcessing - wave decomposition" include("WaveDecompositionTests.jl")
@testset "PostProcessing - R/T/A coefficients" include("CoefficientsTests.jl")
@testset "PostProcessing - rigid-body coefficients" include("RigidBodyCoefficientsTests.jl")
