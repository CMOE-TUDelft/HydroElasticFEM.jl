"""
    module Geometry

Geometry, mesh construction, integration domains, and field helpers for
HydroElasticFEM.jl.

## Two domain families

All geometry in HydroElasticFEM is accessed through the [`AbstractDomain`](@ref)
interface.  Two concrete implementations cover the two main use cases:

**Structured Cartesian** (`TankDomain` / `CartesianDomain`) — use when the
fluid region is a simple rectangular box and structure/damping sub-domains can
be described by axis-aligned coordinates.  Boundaries are identified by
coordinate masks computed at runtime; no external mesh file is needed.

**Unstructured Gmsh** (`GmshDomain`) — use when the geometry is irregular or
comes from a CAD tool.  Boundaries are identified exclusively by Gmsh
*physical-group names* baked into the `.msh` file.  Coordinate-based filtering
is never used, because it is fragile for unstructured meshes and couples
geometry decisions to Julia source code.

Both families expose the same `AbstractDomain` API and feed the same
three-step simulation pipeline described below.

## Simulation pipeline

Every simulation follows this three-step pattern, regardless of domain type:

```
1. model  = build_model(domain)
2. trians = build_triangulations(domain, model)
3. dom    = get_integration_domains(trians; degree=4)
```

**Step 1** constructs (or retrieves) a Gridap `DiscreteModel`.
**Step 2** partitions the model into named sub-triangulations (`TankTriangulations`).
**Step 3** wraps each triangulation in a Gridap `Measure` and outward normal
(`IntegrationDomains`), ready for weak-form assembly.

All three steps must share the same `model` instance.  Building triangulations
from separate model objects will trigger a Gridap assertion error at assembly
time.

## File map

| File                    | Purpose                                              |
|-------------------------|------------------------------------------------------|
| `AbstractDomain.jl`     | `AbstractDomain` interface, `STANDARD_TAGS`          |
| `Triangulations.jl`     | `TankTriangulations` container                       |
| `CartesianDomain.jl`    | `CartesianDomain{D}`, `build_model` / `build_triangulations` (plain box), `map_fn` / `f_z` |
| `TankDomain.jl`         | `TankDomain{D}`, `StructureDomain`, `DampingZone`, `JointDomain`, surface-mask partition, `get_plate_triangulation` |
| `GmshDomain.jl`         | `GmshDomain`, tag-based triangulations, `validate_gmsh_tags` |
| `IntegrationDomains.jl` | `IntegrationDomains` container, `get_integration_domains` |

## Quick-start

```julia
import HydroElasticFEM.Geometry as G

# Structured Cartesian tank with an embedded beam and a damping zone
tank = G.TankDomain(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains = [G.StructureDomain(L=1.0, x₀=[1.5, 1.0])],
    damping_zones     = [G.DampingZone(L=0.5, x₀=[0.0, 1.0])],
)
model  = G.build_model(tank)
trians = G.build_triangulations(tank, model)
dom    = G.get_integration_domains(trians; degree=4)

# Unstructured Gmsh mesh — boundaries come from physical-group names
gmsh_tank = G.GmshDomain("tank.msh"; dim=2)
model     = G.build_model(gmsh_tank)
trians    = G.build_triangulations(gmsh_tank, model)
dom       = G.get_integration_domains(trians)
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
