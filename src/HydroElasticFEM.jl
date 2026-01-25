module HydroElasticFEM

const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

include(joinpath(PKG_ROOT, "src", "Resonator", "Resonator.jl"))
include(joinpath(PKG_ROOT, "src", "WaveInput_FrequencyDomain.jl"))

using .Resonator
using .WaveInput_FrequencyDomain

greet() = print("Hello World!")

end # module HydroElasticFEM
