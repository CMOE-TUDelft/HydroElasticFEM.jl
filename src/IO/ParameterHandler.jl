"""
    module ParameterHandler

Configuration structs for HydroElasticFEM.

Provides lightweight, `@with_kw`-constructed parameter containers
consumed by other modules:

- **`FESpaceConfig`** — numerical FE discretisation parameters (order, conformity,
  vector type, Dirichlet tags/values) stored inside each physics entity.
- **`SimConfig`** — simulation run settings (domain type, frequency, solver).
- **`TimeConfig`** — time-domain integration parameters (time step, final time,
  initial conditions, spectral radius).

Loaded early in the module dependency chain so that both `PhysicsCore`
(entities) and `Simulation` can depend on these types.
"""
module ParameterHandler

using Parameters
using Gridap

include("FESpaceConfig.jl")
include("SimConfig.jl")
include("TimeConfig.jl")

export FESpaceConfig
export SimConfig, TimeConfig

end # module ParameterHandler
