using GridapGmsh

# ─────────────────────────────────────────────────────────────────────────────
# GmshDomain — unstructured mesh wrapper
# ─────────────────────────────────────────────────────────────────────────────

"""
    GmshDomain <: AbstractDomain

Wraps a [`GridapGmsh`](https://github.com/gridap/GridapGmsh.jl) discrete model
loaded from a `.msh` file.

All boundary conditions are identified exclusively by their Gmsh physical-group
name.  The constructor enforces that every [`STANDARD_TAGS`](@ref) name is
present in the mesh before the object is created.

## Constructor

    GmshDomain(msh_file; dim=2, verbose=false)

## Fields (read-only)
- `model`      — the underlying `UnstructuredDiscreteModel`
- `msh_file`   — absolute path of the source `.msh` file
- `tags`       — `Dict{String,String}` mapping every physical-group name to
                 itself (standard and extra tags)
- `dim`        — ambient spatial dimension (2 or 3)

## Example

```julia
dom = GmshDomain("tank.msh"; dim=2)
Ω   = triangulation(dom)
Γfs = get_boundary(dom, "free_surface")
```
"""
struct GmshDomain <: AbstractDomain
  model::DiscreteModel
  msh_file::String
  tags::Dict{String, String}
  dim::Int
end

"""
    GmshDomain(msh_file::String; dim::Int=2, verbose::Bool=false) -> GmshDomain

Load a Gmsh `.msh` file and build a `GmshDomain`.

# Arguments
- `msh_file` — path to the `.msh` file (v2 or v4 format)
- `dim`      — ambient spatial dimension; default 2
- `verbose`  — if `true`, print the list of detected physical groups

# Errors
Raises an `ArgumentError` listing all missing required tags if any of the six
[`STANDARD_TAGS`](@ref) (`"fluid"`, `"free_surface"`, `"seabed"`,
`"inlet"`, `"outlet"`, `"structure"`) are absent from the mesh.

# Example

```julia
dom = GmshDomain("path/to/tank.msh")
dom = GmshDomain("path/to/tank.msh"; dim=3, verbose=true)
```
"""
function GmshDomain(
  msh_file::String;
  dim::Int   = 2,
  verbose::Bool = false,
)
  if !isfile(msh_file)
    error("Mesh file not found: \"$msh_file\"")
  end

  model = GridapGmsh.GmshDiscreteModel(msh_file)
  tags  = _extract_gmsh_tags(model, msh_file)

  if verbose
    println("GmshDomain: loaded \"$msh_file\"")
    println("  Physical groups detected: ", join(sort(collect(keys(tags))), ", "))
  end

  _assert_standard_tags(tags, "GmshDomain(\"$msh_file\")")
  GmshDomain(model, abspath(msh_file), tags, dim)
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _extract_gmsh_tags(model, msh_file) -> Dict{String, String}

Read physical-group names from the Gridap face labeling of `model` and
return an identity-valued dictionary `name => name`.

Tags starting with `"tag_"` are internal Gridap labels and are skipped.
"""
function _extract_gmsh_tags(model, msh_file::String)
  labels   = get_face_labeling(model)
  all_names = labels.tag_to_name
  tags = Dict{String, String}()
  for name in all_names
    # Skip internal Gridap tags (e.g. "tag_1", "interior", "boundary")
    startswith(name, "tag_")   && continue
    name == "interior"         && continue
    name == "boundary"         && continue
    tags[name] = name
  end
  tags
end

# ─────────────────────────────────────────────────────────────────────────────
# AbstractDomain interface
# ─────────────────────────────────────────────────────────────────────────────

"""
    ambient_dimension(d::GmshDomain) -> Int

Return the ambient spatial dimension passed at construction (2 or 3).
"""
ambient_dimension(d::GmshDomain) = d.dim

"""
    manifold_dimension(d::GmshDomain) -> Int

Return the manifold dimension of the fluid bulk (equal to `ambient_dimension`
for volume meshes).
"""
manifold_dimension(d::GmshDomain) = d.dim

"""
    boundary_tags(d::GmshDomain) -> Dict{String, String}

Return all physical-group names extracted from the `.msh` file as an
identity-valued dictionary `name => name`.

All six [`STANDARD_TAGS`](@ref) are guaranteed to be present.
"""
boundary_tags(d::GmshDomain) = d.tags

"""
    triangulation(d::GmshDomain) -> Triangulation

Return the Gridap interior (bulk-fluid) triangulation from the Gmsh model.
"""
function triangulation(d::GmshDomain)
  Interior(d.model)
end

"""
    get_boundary(d::GmshDomain, name::String) -> BoundaryTriangulation

Return the boundary triangulation for the physical group named `name`.

`name` must be one of the keys in [`boundary_tags`](@ref).  All six
[`STANDARD_TAGS`](@ref) are always valid.
"""
function get_boundary(d::GmshDomain, name::String)
  if !haskey(d.tags, name)
    error(
      "Boundary tag \"$name\" not found in GmshDomain loaded from " *
      "\"$(d.msh_file)\". Available tags: " *
      join(sort(collect(keys(d.tags))), ", ") * ".",
    )
  end
  Boundary(d.model, tags=name)
end

# ─────────────────────────────────────────────────────────────────────────────
# Simulation pipeline bridge
# ─────────────────────────────────────────────────────────────────────────────

"""
    build_model(d::GmshDomain) -> UnstructuredDiscreteModel

Return the pre-loaded Gridap discrete model.  This is a zero-cost accessor
(the model is stored in `d.model`).
"""
build_model(d::GmshDomain) = d.model

"""
    build_triangulations(d::GmshDomain, model) -> TankTriangulations

Build a [`TankTriangulations`](@ref) from the Gmsh model using the standard
physical-group names.

The returned container has the following entries:

| Key           | Source                              |
|---------------|-------------------------------------|
| `:Ω`          | `Interior(model)`                   |
| `:Γ`          | `Boundary(model, tags="free_surface")` |
| `:Γfs`        | `Boundary(model, tags="free_surface")` |
| `:Γκ`         | `Boundary(model, tags="free_surface")` |
| `:Γη`         | `Boundary(model, tags="structure")` |
| `:Γbot`       | `Boundary(model, tags="seabed")`    |
| `:Γin`        | `Boundary(model, tags="inlet")`     |
| `:Γout`       | `Boundary(model, tags="outlet")`    |
| `:Γ_structures` | `[Boundary(model, tags="structure")]` |
| `:Γ_dampings` | damping-zone triangulations (see below) |
| `:Λη`         | `Skeleton(Γη)`                      |
| `:Λ_joints`   | `[]`                                |
| `:joint_domains` | `[]`                             |

**Damping zone support**: if the mesh contains physical groups named
`"damping_in"` and/or `"damping_out"`, their boundary triangulations are
collected (in that order) into `:Γ_dampings`.  `get_integration_domains`
then exposes them as `:dΓd_1` / `:nΓd_1` (for `"damping_in"`) and
`:dΓd_2` / `:nΓd_2` (for `"damping_out"`), which matches the domain
symbols expected by `DampingZoneBC(domain=:dΓd_1, ...)`.  The `:Γκ`
triangulation is automatically expanded to include these zones.

Joint domains are not supported for `GmshDomain`; extend by sub-classing
or by operating on the model directly.
"""
function build_triangulations(d::GmshDomain, model)
  Ω    = Interior(model)
  Γfs  = Boundary(model, tags="free_surface")
  Γη   = Boundary(model, tags="structure")
  Γbot = Boundary(model, tags="seabed")
  Γin  = Boundary(model, tags="inlet")
  Γout = Boundary(model, tags="outlet")
  Λη   = Skeleton(Γη)

  # ── Optional damping zone boundaries ────────────────────────────────────
  # Detect "damping_in" and "damping_out" physical groups (in that order so
  # they map to :dΓd_1 and :dΓd_2 in get_integration_domains).
  DAMP_TAGS = ["damping_in", "damping_out"]
  present_damp = filter(t -> haskey(d.tags, t), DAMP_TAGS)

  Γ_dampings = Any[Boundary(model, tags=t) for t in present_damp]

  # :Γκ = free surface ∪ all damping zones (everything except structure)
  κ_tags = vcat(["free_surface"], present_damp)
  Γκ = isempty(present_damp) ? Γfs : Boundary(model, tags=κ_tags)

  TankTriangulations(Dict{Symbol, Any}(
    :Ω            => Ω,
    :Γ            => Γfs,
    :Γfs          => Γfs,
    :Γκ           => Γκ,
    :Γη           => Γη,
    :Γbot         => Γbot,
    :Γin          => Γin,
    :Γout         => Γout,
    :Γ_structures => Any[Γη],
    :Γ_dampings   => Γ_dampings,
    :Λη           => Λη,
    :Λ_joints     => Any[],
    :joint_domains => JointDomain[],
  ))
end
