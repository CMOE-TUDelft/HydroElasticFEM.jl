module HydroElasticFEM

  const PKG_ROOT = normpath(joinpath(@__DIR__, ".."))  # because @__DIR__ here is src/

  # IO / ParameterHandler (config structs, loaded early)
  include(joinpath(@__DIR__, "IO", "ParameterHandler.jl"))
  using .ParameterHandler

  # Geometry
  include(joinpath(@__DIR__, "Geometry", "Geometry.jl"))
  using .Geometry

  # Shared immutable assembly contexts
  include(joinpath(@__DIR__, "AssemblyContexts.jl"))
  using .AssemblyContexts

  # Physics (new canonical types)
  include(joinpath(@__DIR__, "Physics", "Physics.jl"))
  using .Physics

  # Simulation (simulation orchestrator)
  include(joinpath(@__DIR__, "Simulation", "Simulation.jl"))
  using .Simulation

  include(joinpath(PKG_ROOT, "src", "Utilities.jl"))

  ## Utilities.jl functions included here
  # print_properties()
  # map_vertical_GP_for_const_dep()

  # Re-export Physics public API
  export PhysicsParameters, print_parameters
  export BoundaryCondition, FreeBoundary, FixedBoundary
  export PotentialFlow, FreeSurface, Membrane2D, EulerBernoulliBeam
  export ResonatorSingle, resonator_array
  export IntegrationDomains
  export AbstractAssemblyContext, FrequencyAssemblyContext, TimeAssemblyContext
  export domains, frequency, current_time, stabilization_parameter, has_stabilization, with_time

  # Re-export ParameterHandler public API
  export FESpaceConfig, SimConfig, TimeConfig
  export variable_symbol
  export weakform, mass, damping, stiffness, rhs
  export residual, jacobian, jacobian_t, jacobian_tt
  export has_mass_form, has_damping_form, has_stiffness_form, has_rhs_form

  # Re-export Simulation public API (includes FESpaceAssembly + FEOperators)
  # Re-export Simulation public API (includes FESpaceAssembly + FEOperators)
  export build_fe_spaces, build_test_fe_space, build_trial_fe_space
  export SimResult
  export simulate, detect_couplings, build_fe_operator
  export build_frequency_context, build_time_context
  export build_frequency_fe_operator, build_time_fe_operator
  export get_assembly_context
  export FieldMap
  export assemble_weakform
  export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
  export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module HydroElasticFEM
