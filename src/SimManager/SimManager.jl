"""
    module SimManager

High-level simulation orchestrator for HydroElasticFEM.

Given physics entities, triangulations and integration domains, builds
FE spaces, assembles the FE operator (frequency or time-domain), solves,
and returns a `SimResult`.

# Main entry point
- `simulate(config, entities_trians...; dom, ...)` — frequency-domain
- `simulate(config, tconfig, entities_trians...; dom, ...)` — time-domain
"""
module SimManager

using Parameters
using Gridap
using Gridap.ODEs

import ..Geometry as G
import ..PhysicsCore.Entities as E
import ..PhysicsCore.FESpaceAssembly as FA

include("SimConfig.jl")
include("TimeConfig.jl")
include("SimResult.jl")
include("simulate.jl")

export SimConfig, TimeConfig, SimResult
export simulate, detect_couplings

end # module SimManager
