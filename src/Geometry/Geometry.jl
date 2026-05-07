"""
    module Geometry

Geometry, mesh construction, integration domains and field helpers for
HydroElasticFEM.jl.

## File map

| File                  | Purpose                                              |
|-----------------------|------------------------------------------------------|
| `AbstractDomain.jl`   | `AbstractDomain` interface, `STANDARD_TAGS`          |
| `Triangulations.jl`   | `TankTriangulations` container, `_tank_triangulation_dict` |
| `CartesianDomain.jl`  | `CartesianDomain{D}`, centroid helpers, `build_model` / `build_triangulations` for plain Cartesian meshes, `map_fn` / `f_z` |
| `TankDomain.jl`       | `TankDomain{D}`, `StructureDomain`, `DampingZone`, `JointDomain`, surface masking, `build_triangulations` for structured tanks, `get_plate_triangulation` |
| `GmshDomain.jl`       | `GmshDomain`, `build_triangulations` for Gmsh meshes, `validate_gmsh_tags` |
| `IntegrationDomains.jl` | `IntegrationDomains` container, `get_integration_domains` |

## Dependency order (top → bottom)

```
AbstractDomain
    └── Triangulations
            └── CartesianDomain
                    └── TankDomain
                            └── GmshDomain
                                    └── IntegrationDomains
```

## Quick-start

```julia
import HydroElasticFEM.Geometry as G

# Structured Cartesian tank
tank = G.TankDomain(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains = [G.StructureDomain(L=1.0, x₀=[1.5, 1.0])],
    damping_zones     = [G.DampingZone(L=0.5, x₀=[0.0, 1.0])],
)
model  = G.build_model(tank)
trians = G.build_triangulations(tank, model)
dom    = G.get_integration_domains(trians; degree=4)

# Unstructured Gmsh mesh
gmsh_dom = G.GmshDomain("tank.msh"; dim=2)
model    = G.build_model(gmsh_dom)
trians   = G.build_triangulations(gmsh_dom, model)
dom      = G.get_integration_domains(trians)
```
"""
module Geometry
using Parameters
using Gridap

# Load order: each file may use symbols from files above it.
include("AbstractDomain.jl")       # AbstractDomain, STANDARD_TAGS
include("Triangulations.jl")       # TankTriangulations, _tank_triangulation_dict
include("CartesianDomain.jl")      # CartesianDomain{D}, centroid helpers, map_fn/f_z
include("TankDomain.jl")           # TankDomain{D}, StructureDomain, DampingZone, JointDomain
include("GmshDomain.jl")           # GmshDomain, validate_gmsh_tags
include("IntegrationDomains.jl")   # IntegrationDomains, get_integration_domains

export TankTriangulations
export IntegrationDomains
export TankDomain
export StructureDomain, DampingZone
export JointDomain
export AbstractDomain, STANDARD_TAGS
export GmshDomain
export validate_gmsh_tags
export triangulation, boundary_tags, ambient_dimension
export manifold_dimension, get_boundary
export CartesianDomain
export f_z, map_fn
export get_plate_triangulation
export build_model, build_triangulations, get_integration_domains
export surface_mask, surface_masks, joint_mask

end # module Geometry
