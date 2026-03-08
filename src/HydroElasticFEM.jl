module HydroElasticFEM

  const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

  # PhysicsCore (new canonical types)
  include(joinpath(@__DIR__, "PhysicsCore", "PhysicalEntities.jl"))
  using .PhysicalEntities

  include(joinpath(@__DIR__, "PhysicsCore", "WeakFormAssembly.jl"))
  using .WeakFormAssembly

  # Backward-compatible shims (kept during transition)
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

  # Re-export PhysicalEntities public API
  export PhysicsParameters, print_parameters
  export BoundaryCondition, FreeBoundary, FixedBoundary
  export AbstractStructure, PotentialFlow, Membrane2D, Beam2D
  export ResonatorSingle, resonator_array
  export WeakFormDomains
  export η_symbol
  export weakform, mass, damping, stiffness, rhs
  export residual, jacobian, jacobian_t, jacobian_tt

  # Re-export WeakFormAssembly public API
  export FieldDict
  export assemble_weakform
  export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
  export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module HydroElasticFEM
