module HydroElasticFEM

const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

include(joinpath(PKG_ROOT, "src", "StructuralComponents", "BeamNoJoints.jl"))
include(joinpath(PKG_ROOT, "src", "StructuralComponents", "Resonator.jl"))
include(joinpath(PKG_ROOT, "src", "StructuralComponents", "Membrane.jl"))
include(joinpath(PKG_ROOT, "src", "Utilities.jl"))
include(joinpath(PKG_ROOT, "src", "WaveInput_FrequencyDomain.jl"))

using .MeshModifier
using .BeamNoJoints
using .Resonator
using .Membrane
using .WaveInput_FrequencyDomain

greet() = print("Hello World!")

## print properties functions
# --------------------Start--------------------

function print_properties(ele::BeamNoJoints.Beam2D)
    BeamNoJoints.print_properties(ele)
end

function print_properties(ele::Membrane.Membrane2D)
    print_properties(ele)
end

function print_properties(ele::Resonator.Single)
    print_properties(ele)
end
# ----------------------End---------------------

end # module HydroElasticFEM
