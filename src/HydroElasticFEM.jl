module HydroElasticFEM

const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

include(joinpath(PKG_ROOT, "src", "StructuralComponents", "Resonator.jl"))
include(joinpath(PKG_ROOT, "src", "StructuralComponents", "Membrane.jl"))
include(joinpath(PKG_ROOT, "src", "Utilities.jl"))
include(joinpath(PKG_ROOT, "src", "WaveInput_FrequencyDomain.jl"))

using .MeshModifier
using .Resonator
using .Membrane
using .WaveInput_FrequencyDomain

greet() = print("Hello World!")

end # module HydroElasticFEM
