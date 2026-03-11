"""
    module Simulation

High-level simulation orchestrator for HydroElasticFEM.

Given physics entities, triangulations and integration domains, builds
FE spaces, assembles the FE operator (frequency or time-domain), solves,
and returns a `SimResult`.

# Main entry point
- `simulate(config, entities_trians...; dom, ...)` — frequency-domain
- `simulate(config, tconfig, entities_trians...; dom, ...)` — time-domain
"""
module Simulation

using Parameters
using Gridap
using Gridap.ODEs

import ..Geometry as G
import ..PhysicsCore.Entities as E
import ..PhysicsCore.FESpaceAssembly as FA

# FEOperators (FieldMap, assemble_*, detect_couplings)
include("FEOperators.jl")
using .FEOperators

include("SimConfig.jl")
include("TimeConfig.jl")
include("SimResult.jl")
include("simulate.jl")

export SimConfig, TimeConfig, SimResult
export simulate
export FEOperators
export FieldMap, detect_couplings
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module Simulation
