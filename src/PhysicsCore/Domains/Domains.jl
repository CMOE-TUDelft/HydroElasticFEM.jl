"""
    module Domains

Domain containers for HydroElasticFEM weak form assembly.

Provides `WeakFormDomains` (dict-based container for Gridap measures,
normals, and metadata) and `FieldDict` (symbol-indexed wrapper around
Gridap multi-field tuples).
"""
module Domains

include("WeakFormDomains.jl")
include("FieldDict.jl")

export WeakFormDomains
export FieldDict

end # module Domains
