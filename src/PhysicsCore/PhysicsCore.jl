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
include("Entities/WeakFormAssembly.jl")
using .WeakFormAssembly

# 5. FESpace assembly (depends on Entities)
include("FESpaces/FESpaceAssembly.jl")
using .FESpaceAssembly

# ─────────────────────────────────────────────────────────────
# Re-export Domains
# ─────────────────────────────────────────────────────────────
export WeakFormDomains, FieldDict

# ─────────────────────────────────────────────────────────────
# Re-export FESpaces
# ─────────────────────────────────────────────────────────────
export FESpaceConfig
export build_fe_spaces, build_test_fe_space, build_trial_fe_space

# ─────────────────────────────────────────────────────────────
# Re-export Entities
# ─────────────────────────────────────────────────────────────
export PhysicsParameters, print_parameters
export BoundaryCondition, FreeBoundary, FixedBoundary
export AbstractStructure, PotentialFlow, FreeSurface, Membrane2D, EulerBernoulliBeam
export ResonatorSingle, resonator_array
export variable_symbol
export weakform, mass, damping, stiffness, rhs
export residual, jacobian, jacobian_t, jacobian_tt
export has_mass_form, has_damping_form, has_stiffness_form, has_rhs_form

# ─────────────────────────────────────────────────────────────
# Re-export WeakFormAssembly
# ─────────────────────────────────────────────────────────────
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module PhysicsCore
