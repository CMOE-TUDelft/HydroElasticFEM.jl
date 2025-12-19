module HydroElasticFEM

const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

include(joinpath(PKG_ROOT, "src", "Resonator", "Resonator.jl"))

using .Resonator


greet() = print("Hello World!")

end # module HydroElasticFEM
