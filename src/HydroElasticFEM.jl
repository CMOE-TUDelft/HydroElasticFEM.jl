module HydroElasticFEM

  const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

  include(joinpath(PKG_ROOT, "src", "StructuralComponents", "BeamNoJoints.jl"))
  include(joinpath(PKG_ROOT, "src", "StructuralComponents", "Resonator.jl"))
  include(joinpath(PKG_ROOT, "src", "StructuralComponents", "Membrane.jl"))
  include(joinpath(PKG_ROOT, "src", "WaveInput_FrequencyDomain.jl"))

  include(joinpath(PKG_ROOT, "src", "Utilities.jl"))

  using .BeamNoJoints
  using .Resonator
  using .Membrane
  using .WaveInput_FrequencyDomain

  ## Utilities.jl functions included here
  # print_properties()
  # map_vertical_GP_for_const_dep()

  greet() = print("Hello World!")

end # module HydroElasticFEM
