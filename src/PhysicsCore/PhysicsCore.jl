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

# ─────────────────────────────────────────────────────────────
# Re-export PhysicalEntities public API
# ─────────────────────────────────────────────────────────────
export PhysicsParameters, print_parameters
export BoundaryCondition, FreeBoundary, FixedBoundary
export AbstractStructure, PotentialFlow, Membrane2D, EulerBernoulliBeam
export ResonatorSingle, resonator_array
export WeakFormDomains
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

end # module PhysicsCore
