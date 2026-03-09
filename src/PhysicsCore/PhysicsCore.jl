"""
    module PhysicsCore

Wrapper module for all physics-related types and assembly routines.

Loads submodules in dependency order:
1. **Domains** — `WeakFormDomains`, `FieldDict` (standalone)
2. **FESpaces** — `FESpaceConfig` (standalone)
3. **Entities** — physics types, traits, composed weak forms (uses Domains + FESpaces)
4. **WeakFormAssembly** — `assemble_*` helpers (uses Entities + Domains)
5. **FESpaceAssembly** — `build_fe_spaces` (uses Entities)
"""
module PhysicsCore

# 1. Domains (standalone)
include("Domains/Domains.jl")
using .Domains

# 2. FESpaces config (standalone)
include("FESpaces/FESpaces.jl")
using .FESpaces

# 3. Entities (depends on Domains + FESpaces)
include("Entities/Entities.jl")
using .Entities

# 4. WeakForm assembly (depends on Entities + Domains)
include("WeakFormAssembly.jl")
using .WeakFormAssembly

# 5. FESpace assembly (depends on Entities)
include("FESpaceAssembly.jl")
using .FESpaceAssembly

export Domains
export Entities
export WeakFormAssembly
export FESpaces
export FESpaceAssembly

end # module PhysicsCore
