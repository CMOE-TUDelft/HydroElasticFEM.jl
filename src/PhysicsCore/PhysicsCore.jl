"""
    module PhysicsCore

Wrapper module for all physics-related types and assembly routines.

Re-exports the public API from `PhysicalEntities` and `WeakFormAssembly`.
"""
module PhysicsCore

include("PhysicalEntities.jl")
using .PhysicalEntities

include("WeakFormAssembly.jl")
using .WeakFormAssembly

include("FESpaceAssembly.jl")
using .FESpaceAssembly

# ─────────────────────────────────────────────────────────────
# Re-export PhysicalEntities public API
# ─────────────────────────────────────────────────────────────
export PhysicsParameters, print_parameters
export BoundaryCondition, FreeBoundary, FixedBoundary
export AbstractStructure, PotentialFlow, FreeSurface, Membrane2D, EulerBernoulliBeam
export ResonatorSingle, resonator_array
export WeakFormDomains, FESpaceConfig
export variable_symbol
export weakform, mass, damping, stiffness, rhs
export residual, jacobian, jacobian_t, jacobian_tt

# ─────────────────────────────────────────────────────────────
# Re-export WeakFormAssembly public API
# ─────────────────────────────────────────────────────────────
export FieldDict
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

# ─────────────────────────────────────────────────────────────
# Re-export FESpaceAssembly public API
# ─────────────────────────────────────────────────────────────
export build_fe_spaces, build_test_fe_space, build_trial_fe_space

end # module PhysicsCore
