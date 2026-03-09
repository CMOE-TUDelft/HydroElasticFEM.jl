module HydroElasticFEM

  const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

  # PhysicsCore (new canonical types)
  include(joinpath(@__DIR__, "PhysicsCore", "PhysicsCore.jl"))
  using .PhysicsCore

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

  # Re-export PhysicsCore public API
  export PhysicsParameters, print_parameters
  export BoundaryCondition, FreeBoundary, FixedBoundary
  export AbstractStructure, PotentialFlow, FreeSurface, Membrane2D, EulerBernoulliBeam
  export ResonatorSingle, resonator_array
  export WeakFormDomains, FESpaceConfig
  export variable_symbol
  export weakform, mass, damping, stiffness, rhs
  export residual, jacobian, jacobian_t, jacobian_tt
  export has_mass_form, has_damping_form, has_stiffness_form, has_rhs_form
  export FieldDict
  export assemble_weakform
  export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
  export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt
  export build_fe_spaces, build_test_fe_space, build_trial_fe_space

end # module HydroElasticFEM
