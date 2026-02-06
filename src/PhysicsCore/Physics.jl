"""
    module Physics

This module defines the core interface for defining models parameters in the HydroElasticFEM package.
"""
module Physics

using Parameters
using Printf

export PhysicsParameters, print_parameters
export EBBeamParameters

"""
    abstract type PhysicsParameters

Abstract base type for all physics parameter structures. 
Any specific physics implementation (e.g., Fluid, Structural) should define a struct that subtypes `PhysicsParameters`.
"""
abstract type PhysicsParameters end

"""
    print_parameters(params::PhysicsParameters)

Print the parameters associated with a specific physics implementation.
This function must be implemented for any concrete subtype of `PhysicsParameters`.

# Arguments
- `params::PhysicsParameters`: The physics parameters to be printed.

# Throws
- `ErrorException`: If the method is not implemented for the specific subtype of `PhysicsParameters`.
"""
function print_parameters(params::PhysicsParameters)
    error("print_parameters not implemented for $(typeof(params))")
end

include("EulerBernoulliBeam.jl")

end