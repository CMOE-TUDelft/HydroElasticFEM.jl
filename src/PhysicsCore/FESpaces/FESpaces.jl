"""
    module FESpaces

FE space configuration for HydroElasticFEM.

Provides `FESpaceConfig`, a struct holding numerical FE discretisation
parameters (polynomial order, conformity, vector type, Nitsche penalty,
Dirichlet tags/values) stored inside each physics entity.
"""
module FESpaces

using Parameters
using Gridap

include("FESpaceConfig.jl")

export FESpaceConfig

end # module FESpaces
