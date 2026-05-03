# ─────────────────────────────────────────────────────────────────────────────
# AbstractDomain — unified geometry interface for HydroElasticFEM
# ─────────────────────────────────────────────────────────────────────────────

"""
    const STANDARD_TAGS

Canonical physical-group names every domain must expose.

| Tag name         | Meaning                              | Dim (2D) |
|------------------|--------------------------------------|----------|
| `"fluid"`        | Bulk fluid volume                    | 2        |
| `"free_surface"` | Top open-water boundary              | 1        |
| `"seabed"`       | Bottom boundary                      | 1        |
| `"inlet"`        | Wave-generation / left wall          | 1        |
| `"outlet"`       | Absorbing / right wall               | 1        |
| `"structure"`    | Fluid–structure interface            | 1        |
"""
const STANDARD_TAGS = [
  "fluid",
  "free_surface",
  "seabed",
  "inlet",
  "outlet",
  "structure",
]

# ─────────────────────────────────────────────────────────────────────────────
# Abstract type
# ─────────────────────────────────────────────────────────────────────────────

"""
    abstract type AbstractDomain

Unified interface for geometry domains used in HydroElasticFEM simulations.

Every concrete subtype must implement the five required methods below.
Two concrete implementations are provided:

- [`TankDomain2D`](@ref) — structured Cartesian rectangular tank
- [`GmshDomain`](@ref)   — unstructured mesh loaded from a `.msh` file

## Required interface

| Method                              | Returns                        |
|-------------------------------------|--------------------------------|
| `triangulation(d)`                  | Gridap `Triangulation` (bulk)  |
| `boundary_tags(d)`                  | `Dict{String,Any}`             |
| `ambient_dimension(d)`              | `Int` (2 or 3)                 |
| `manifold_dimension(d)`             | `Int` (same as ambient for vol)|
| `get_boundary(d, name)`             | Gridap `BoundaryTriangulation` |

All six [`STANDARD_TAGS`](@ref) must be accessible via `get_boundary`.

## Example

```julia
# Works with any AbstractDomain subtype
function describe_domain(d::AbstractDomain)
    println("dim = ", ambient_dimension(d))
    println("tags = ", keys(boundary_tags(d)))
    Γfs = get_boundary(d, "free_surface")
    println("free-surface cells = ", num_cells(Γfs))
end
```
"""
abstract type AbstractDomain end

# ─────────────────────────────────────────────────────────────────────────────
# Interface stubs (throw informative errors when not implemented)
# ─────────────────────────────────────────────────────────────────────────────

"""
    triangulation(d::AbstractDomain) -> Triangulation

Return the Gridap bulk (interior) triangulation for domain `d`.

The returned object corresponds to the full fluid-volume mesh.  For 2D
problems this is the 2-D interior; for 3D the 3-D interior.
"""
function triangulation(d::AbstractDomain)
  error(
    "triangulation not implemented for $(typeof(d)). " *
    "Subtypes of AbstractDomain must implement this method.",
  )
end

"""
    boundary_tags(d::AbstractDomain) -> Dict{String, Any}

Return a dictionary mapping tag name strings to Gridap tag objects (strings,
integers, or triangulations depending on the domain type).

All six [`STANDARD_TAGS`](@ref) must be present as keys.  Additional
domain-specific tags may also appear.
"""
function boundary_tags(d::AbstractDomain)
  error(
    "boundary_tags not implemented for $(typeof(d)). " *
    "Subtypes of AbstractDomain must implement this method.",
  )
end

"""
    ambient_dimension(d::AbstractDomain) -> Int

Return the ambient spatial dimension of the domain (2 or 3).
"""
function ambient_dimension(d::AbstractDomain)
  error(
    "ambient_dimension not implemented for $(typeof(d)). " *
    "Subtypes of AbstractDomain must implement this method.",
  )
end

"""
    manifold_dimension(d::AbstractDomain) -> Int

Return the manifold dimension of the domain.  For volume meshes this equals
[`ambient_dimension`](@ref).
"""
function manifold_dimension(d::AbstractDomain)
  error(
    "manifold_dimension not implemented for $(typeof(d)). " *
    "Subtypes of AbstractDomain must implement this method.",
  )
end

"""
    get_boundary(d::AbstractDomain, name::String) -> BoundaryTriangulation

Return the Gridap boundary triangulation for the physical region named `name`.

`name` must be one of the keys returned by [`boundary_tags`](@ref).  All six
[`STANDARD_TAGS`](@ref) must be supported.

## Example

```julia
Γfs = get_boundary(domain, "free_surface")
Γin = get_boundary(domain, "inlet")
```
"""
function get_boundary(d::AbstractDomain, name::String)
  error(
    "get_boundary not implemented for $(typeof(d)). " *
    "Subtypes of AbstractDomain must implement this method.",
  )
end

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers available to all subtypes
# ─────────────────────────────────────────────────────────────────────────────

"""
    _assert_standard_tags(tags::Dict, domain_description::String)

Check that all [`STANDARD_TAGS`](@ref) are present as keys in `tags`.
Raises an `ArgumentError` listing every missing tag.
"""
function _assert_standard_tags(tags::Dict, domain_description::String)
  missing = filter(t -> !haskey(tags, t), STANDARD_TAGS)
  if !isempty(missing)
    error(
      "$(domain_description) is missing required standard tags: " *
      join(missing, ", ") *
      ". All of $(STANDARD_TAGS) must be defined.",
    )
  end
end
