"""
    module PhysicalEntities

Unified type hierarchy for all physical entities in HydroElasticFEM.

Follows the PhysicsCore pattern: abstract base type `PhysicsParameters`,
concrete structs for each entity, and `print_parameters` as the generic
display interface.
"""
module PhysicalEntities

using Parameters
using Gridap
using Printf

# ─────────────────────────────────────────────────────────────
# Abstract base
# ─────────────────────────────────────────────────────────────

"""
    abstract type PhysicsParameters

Abstract base type for all physics parameter structures.
"""
abstract type PhysicsParameters end

"""
    print_parameters(params::PhysicsParameters)

Print the parameters of a concrete `PhysicsParameters` subtype.
Must be implemented for every concrete subtype.
"""
function print_parameters(params::PhysicsParameters)
    error("print_parameters not implemented for $(typeof(params))")
end

# ─────────────────────────────────────────────────────────────
# Boundary conditions
# ─────────────────────────────────────────────────────────────

abstract type BoundaryCondition end
struct FreeBoundary  <: BoundaryCondition end
struct FixedBoundary <: BoundaryCondition end

# ─────────────────────────────────────────────────────────────
# Structural entities
# ─────────────────────────────────────────────────────────────

abstract type AbstractStructure <: PhysicsParameters end

include("Membrane2D.jl")
include("Beam2D.jl")
include("Resonator.jl")

# ─────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────

export PhysicsParameters, print_parameters
export BoundaryCondition, FreeBoundary, FixedBoundary
export AbstractStructure, Membrane2D, Beam2D
export ResonatorSingle, resonator_array

end # module PhysicalEntities
