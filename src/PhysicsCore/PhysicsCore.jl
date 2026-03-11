"""
    module PhysicsCore

Wrapper module for physics-related types.

Loads submodules in dependency order:
1. **Geometry** (parent) — `IntegrationDomains`
2. **FESpaces** — `FESpaceConfig` (standalone)
3. **Entities** — physics types, traits, composed weak forms (uses Geometry + FESpaces)

Note: `FESpaceAssembly` and `FEOperators` live in the `Simulation` sibling module.
"""
module PhysicsCore

# 1. IntegrationDomains comes from the Geometry sibling module
using ..Geometry

# 2. FESpaces config (standalone)
include("FESpaces/FESpaces.jl")
using .FESpaces

# 3. Entities (depends on Geometry + FESpaces)
include("Entities/Entities.jl")
using .Entities

export Entities
export FESpaces

end # module PhysicsCore
