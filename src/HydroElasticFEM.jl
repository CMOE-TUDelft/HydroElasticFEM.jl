"""
    HydroElasticFEM

Finite element analysis of hydroelastic wave–structure interaction
for very large floating structures.

Uses a continuous/discontinuous Galerkin (C/DG) formulation coupling linearised
potential flow fluid models with structural models (membranes, Euler–Bernoulli beams,
Kirchhoff–Love plates, Timoshenko beams, resonators), implemented using
[Gridap.jl](https://github.com/gridap/Gridap.jl) as a base finite element library.

# Quick start

```julia
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S

# Build 2-D tank domain
domain = G.TankDomain(L=10.0, H=1.0, nx=60, ny=8)

# Define physics entities
fluid   = P.PotentialFlow()
surface = P.FreeSurface()

# Run frequency-domain simulation at ω = 1 rad/s
config  = PH.FreqDomainConfig(ω=1.0)
problem = S.build_problem(domain, P.PhysicsParameters[fluid, surface], config)
result  = S.simulate(problem)

# Unpack solution fields: velocity potential ϕ, free-surface elevation κ
phi_h, kappa_h = result.solution
```

See the `examples/` directory for complete simulation scripts, and the online
documentation for the theory, API reference, and validation benchmarks.

# Reference
Colomés, O., Verdugo, F., & Akkerman, I. (2023). A monolithic finite element
formulation for the hydroelastic analysis of very large floating structures.
*Int. J. Numer. Methods Eng.*, 124(3), 714–751. DOI: 10.1002/nme.7140
"""
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

  # Re-export Geometry public API
  export AbstractDomain, STANDARD_TAGS
  export TankDomain
  export StructureDomain, DampingZone, JointDomain
  export CartesianDomain
  export GmshDomain
  export triangulation, boundary_tags, ambient_dimension
  export manifold_dimension, get_boundary
  export get_plate_triangulation
  export validate_gmsh_tags
  export build_model, build_triangulations, get_integration_domains
  export TankTriangulations, IntegrationDomains

  # Re-export Physics public API
  export PhysicsParameters, print_parameters
  export PotentialFlow, FreeSurface, Membrane, EulerBernoulliBeam
  export KirchhoffLovePlate
  export TimoshenkoBeam
  export build_kl_tensor, build_KL_tensor
  export equivalent_beam_rigidity
  export ResonatorSingle, resonator_array
  export IntegrationDomains
  export AbstractAssemblyContext, FrequencyAssemblyContext, TimeAssemblyContext
  export domains, frequency, current_time, stabilization_parameter, has_stabilization, with_time

  # Re-export ParameterHandler public API
  export FESpaceConfig, TimeConfig
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
