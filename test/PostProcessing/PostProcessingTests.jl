using Test

@testset "PostProcessing - wave decomposition" include("WaveDecompositionTests.jl")
@testset "PostProcessing - R/T/A coefficients" include("CoefficientsTests.jl")
