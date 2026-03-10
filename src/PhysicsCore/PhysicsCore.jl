"""
    module PhysicsCore

Wrapper module for all physics-related types and assembly routines.

Loads submodules in dependency order:
1. **Geometry** (parent) — `IntegrationDomains`, `FieldDict`
2. **FESpaces** — `FESpaceConfig` (standalone)
3. **Entities** — physics types, traits, composed weak forms (uses Geometry + FESpaces)
4. **WeakFormAssembly** — `assemble_*` helpers (uses Entities + Geometry)
5. **FESpaceAssembly** — `build_fe_spaces` (uses Entities)
"""
module PhysicsCore

# 1. IntegrationDomains & FieldDict come from the Geometry sibling module
using ..Geometry

# 2. FESpaces config (standalone)
include("FESpaces/FESpaces.jl")
using .FESpaces

# 3. Entities (depends on Geometry + FESpaces)
include("Entities/Entities.jl")
using .Entities

# 4. WeakForm assembly (depends on Entities + Geometry)
include("WeakFormAssembly.jl")
using .WeakFormAssembly

# 5. FESpace assembly (depends on Entities)
include("FESpaceAssembly.jl")
using .FESpaceAssembly

export Entities
export WeakFormAssembly
export FESpaces
export FESpaceAssembly

end # module PhysicsCore
