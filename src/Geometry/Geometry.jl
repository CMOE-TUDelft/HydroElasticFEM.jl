"""
    module Geometry

Geometry, mesh construction, integration domains and field helpers.

Provides:
- `TankDomain2D`, `StructureDomain1D`, `DampingZone1D` — geometry specs
- `build_model`, `build_triangulations` — mesh construction
- `IntegrationDomains` — dict-based container for Gridap measures/normals
- `FieldDict` — symbol-indexed wrapper around Gridap field tuples
- `get_integration_domains` — build `IntegrationDomains` from triangulations
"""
module Geometry
using Parameters
using Gridap

include("IntegrationDomains.jl")
include("FieldDict.jl")
include("CartesianGeometry.jl")

export IntegrationDomains
export FieldDict

end # module Geometry