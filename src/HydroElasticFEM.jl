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

  # ── Geometry public API ────────────────────────────────────────────────────
  # Domain constructors — users build one of these to describe the geometry
  export AbstractDomain, STANDARD_TAGS
  export TankDomain
  export StructureDomain, DampingZone, JointDomain
  export CartesianDomain
  export GmshDomain
  # Domain query helpers
  export boundary_tags, ambient_dimension
  export get_boundary
  export get_plate_triangulation
  export validate_gmsh_tags
  # Geometry pipeline (build_model → build_triangulations → get_integration_domains)
  export build_model, build_triangulations, get_integration_domains
  export TankTriangulations, IntegrationDomains

  # ── Physics public API ──────────────────────────────────────────────────────
  # Entity types — users construct instances of these
  export PhysicsParameters, print_parameters
  export PotentialFlow, FreeSurface, Membrane, EulerBernoulliBeam
  export KirchhoffLovePlate
  export TimoshenkoBeam
  export ResonatorSingle, resonator_array
  # Physics helper functions
  export build_kl_tensor, build_KL_tensor
  export equivalent_beam_rigidity
  # Extension interface — needed by users implementing new physics models
  export variable_symbol
  export has_mass_form, has_damping_form, has_stiffness_form, has_rhs_form
  # Assembly context types — needed for custom physics dispatch
  export AbstractAssemblyContext, FrequencyAssemblyContext, TimeAssemblyContext

  # ── ParameterHandler public API ────────────────────────────────────────────
  # Configuration types — users construct these to set up simulations
  export FESpaceConfig, TimeConfig
  export FreqDomainConfig, TimeDomainConfig

  # ── Simulation public API ───────────────────────────────────────────────────
  # Primary simulation entry points — users call build_problem then simulate
  export build_problem
  export SimResult
  export simulate

end # module HydroElasticFEM
