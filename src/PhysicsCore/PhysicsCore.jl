"""
    module PhysicsCore

Wrapper module for physics-related types.

Loads submodules in dependency order:
1. **Geometry** (parent) — `IntegrationDomains`
2. **ParameterHandler** (parent) — `FESpaceConfig`, `SimConfig`, `TimeConfig`
3. **Entities** — physics types, traits, composed weak forms (uses Geometry + ParameterHandler)

Note: `FESpaceAssembly` and `FEOperators` live in the `Simulation` sibling module.
"""
module PhysicsCore

# 1. IntegrationDomains comes from the Geometry sibling module
using ..Geometry

# 2. Config structs from ParameterHandler sibling module
using ..ParameterHandler

# 3. Entities (depends on Geometry + ParameterHandler)
include("Entities/Entities.jl")
using .Entities

export Entities

end # module PhysicsCore
