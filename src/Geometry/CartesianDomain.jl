# ─────────────────────────────────────────────────────────────────────────────
# CartesianDomain — generic axis-aligned structured domain
# ─────────────────────────────────────────────────────────────────────────────

"""
    CartesianDomain{D} <: AbstractDomain

Generic axis-aligned Cartesian domain for 2D or 3D fluid problems.

`CartesianDomain{2}` uses coordinates `(x, z)` (horizontal, vertical).
`CartesianDomain{3}` uses coordinates `(x, y, z)`.

This type is functional and dispatch-driven: dimension-specific behavior is
selected with `Val(D)` dispatch so the same API works for both 2D and 3D.

## Fields
- `mins::NTuple{D,Float64}` — lower bounds for each axis [m]
- `maxs::NTuple{D,Float64}` — upper bounds for each axis [m]
- `parts::NTuple{D,Int}`    — number of cells along each axis
- `map::Function`           — optional coordinate map (default: identity)

## Constructors

Use the keyword constructor [`CartesianDomain(; ...)`](@ref) instead of the
inner constructor.

## Example

```julia
# 2D box [0,4] × [0,1], 80 × 10 cells
d = CartesianDomain(L=4.0, H=1.0, nx=80, ny=10)

# 3D box [0,8] × [-2,2] × [0,2], graded in z
d = CartesianDomain(
    mins=(0.0,-2.0,0.0), maxs=(8.0,2.0,2.0), parts=(40,20,10),
    map = x -> G.map_fn(x, 2.0, 10))
```
"""
struct CartesianDomain{D, F <: Function} <: AbstractDomain
  mins::NTuple{D, Float64}
  maxs::NTuple{D, Float64}
  parts::NTuple{D, Int}
  map::F
end

"""
    CartesianDomain(; L, H, nx, ny, W=nothing, nz=nothing,
                     mins=nothing, maxs=nothing, parts=nothing,
                     map=x->x) -> CartesianDomain

Build a Cartesian domain.  Three calling forms are supported:

**Shorthand 2D** — axis-aligned box `[0,L]×[0,H]`:
```julia
CartesianDomain(L=4.0, H=1.0, nx=80, ny=10)
```

**Shorthand 3D** — axis-aligned box `[0,L]×[0,W]×[0,H]`:
```julia
CartesianDomain(L=8.0, W=4.0, H=2.0, nx=40, ny=20, nz=10)
```

**Explicit bounds** — arbitrary origin, shifted, or non-unit aspect ratios:
```julia
CartesianDomain(
    mins=(0.0, -2.0, 0.0),
    maxs=(8.0,  2.0, 2.0),
    parts=(40, 20, 10),
    map = x -> map_fn(x, 2.0, 10),
)
```

The `map` keyword is a `Function` applied to every node coordinate by
`CartesianDiscreteModel`.  Use [`map_fn`](@ref) for vertical grading.
"""
function CartesianDomain(;
  L = nothing, H = nothing, nx = nothing, ny = nothing,
  W = nothing, nz = nothing,
  mins = nothing, maxs = nothing, parts = nothing,
  map = x -> x,
)
  if !isnothing(mins) || !isnothing(maxs) || !isnothing(parts)
    isnothing(mins) && error("`mins` must be provided for explicit CartesianDomain construction.")
    isnothing(maxs) && error("`maxs` must be provided for explicit CartesianDomain construction.")
    isnothing(parts) && error("`parts` must be provided for explicit CartesianDomain construction.")
    D = length(mins)
    length(maxs) == D || error("length(maxs) must match length(mins).")
    length(parts) == D || error("length(parts) must match length(mins).")
    return CartesianDomain{D, typeof(map)}(
      Tuple(Float64.(mins)),
      Tuple(Float64.(maxs)),
      Tuple(Int.(parts)),
      map,
    )
  end

  if isnothing(W) && isnothing(nz)
    return CartesianDomain{2, typeof(map)}(
      (0.0, 0.0),
      (Float64(L), Float64(H)),
      (Int(nx), Int(ny)),
      map,
    )
  end

  if !isnothing(W) && !isnothing(nz)
    return CartesianDomain{3, typeof(map)}(
      (0.0, 0.0, 0.0),
      (Float64(L), Float64(W), Float64(H)),
      (Int(nx), Int(ny), Int(nz)),
      map,
    )
  end

  error(
    "CartesianDomain expects 2D args (L,H,nx,ny), " *
    "3D args (L,W,H,nx,ny,nz), or explicit (mins,maxs,parts).",
  )
end

# ─────────────────────────────────────────────────────────────────────────────
# Low-level cell-coordinate helpers (used here and in TankDomain)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _centroid(xs) -> VectorValue

Compute the centroid of a cell given its `Vector{VectorValue}` of node
coordinates.
"""
function _centroid(xs)
  n = length(xs)
  (1 / n) * sum(xs)
end

"""
    _always_false_mask() -> Function

Return a cell predicate that always returns `false`.
Used as a placeholder for empty structural regions.
"""
_always_false_mask() = xs -> false

"""
    _centroid_axis_eq_mask(axis, value; tol=1e-10) -> Function

Return a cell predicate that is `true` when the centroid's coordinate
along `axis` equals `value` within `tol`.
"""
function _centroid_axis_eq_mask(axis::Int, value; tol = 1.0e-10)
  xs -> abs(_centroid(xs)[axis] - value) < tol
end

"""
    _centroid_axis_either_eq_mask(axis, value_a, value_b; tol=1e-10) -> Function

Return a cell predicate that is `true` when the centroid's `axis` coordinate
equals either `value_a` or `value_b` within `tol`.
"""
function _centroid_axis_either_eq_mask(axis::Int, value_a, value_b; tol = 1.0e-10)
  function (xs)
    c = _centroid(xs)
    abs(c[axis] - value_a) < tol || abs(c[axis] - value_b) < tol
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# AbstractDomain interface
# ─────────────────────────────────────────────────────────────────────────────

"""
    ambient_dimension(d::CartesianDomain{D}) -> Int

Return the ambient spatial dimension `D`.
"""
ambient_dimension(::CartesianDomain{D}) where {D} = D

"""
    manifold_dimension(d::CartesianDomain{D}) -> Int

Return the manifold dimension `D` (equals ambient for volume meshes).
"""
manifold_dimension(::CartesianDomain{D}) where {D} = D

"""
    boundary_tags(d::CartesianDomain{D}) -> Dict{String, String}

Return a dictionary of valid boundary tag names for a Cartesian domain of
dimension `D`.  All six [`STANDARD_TAGS`](@ref) are present; 3D adds
`"lateral_walls"`.
"""
boundary_tags(::CartesianDomain{D}) where {D} = _cartesian_boundary_tags(Val(D))

function _cartesian_boundary_tags(::Val{2})
  Dict{String, String}(
    "fluid" => "fluid",
    "free_surface" => "free_surface",
    "seabed" => "seabed",
    "inlet" => "inlet",
    "outlet" => "outlet",
    "structure" => "structure",
  )
end

function _cartesian_boundary_tags(::Val{3})
  Dict{String, String}(
    "fluid" => "fluid",
    "free_surface" => "free_surface",
    "seabed" => "seabed",
    "inlet" => "inlet",
    "outlet" => "outlet",
    "structure" => "structure",
    "lateral_walls" => "lateral_walls",
  )
end

# ─────────────────────────────────────────────────────────────────────────────
# Model and triangulation construction
# ─────────────────────────────────────────────────────────────────────────────

_flatten_bounds(mins::NTuple{2, Float64}, maxs::NTuple{2, Float64}) =
  (mins[1], maxs[1], mins[2], maxs[2])

_flatten_bounds(mins::NTuple{3, Float64}, maxs::NTuple{3, Float64}) =
  (mins[1], maxs[1], mins[2], maxs[2], mins[3], maxs[3])

"""
    build_model(d::CartesianDomain{D}) -> CartesianDiscreteModel

Build a `CartesianDiscreteModel` from `d`.

Each call constructs a new model object.  For performance-critical code that
needs multiple triangulations from the same mesh, call `build_model` once
and pass the model to [`build_triangulations`](@ref) and
[`get_boundary`](@ref).
"""
function build_model(d::CartesianDomain{D}) where {D}
  bounds = _flatten_bounds(d.mins, d.maxs)
  CartesianDiscreteModel(bounds, d.parts; map = d.map)
end

"""
    triangulation(d::CartesianDomain) -> Triangulation

Build the model and return the bulk-fluid interior triangulation.

Note: builds a fresh model each call.  Use [`build_model`](@ref) /
[`build_triangulations`](@ref) when multiple triangulations are needed.
"""
triangulation(d::CartesianDomain) = Interior(build_model(d))

# ─────────────────────────────────────────────────────────────────────────────
# Boundary extraction helpers
# ─────────────────────────────────────────────────────────────────────────────

_vertical_axis(::Val{2}) = 2
_vertical_axis(::Val{3}) = 3

_cartesian_valid_tags(::Val{2}) = STANDARD_TAGS
_cartesian_valid_tags(::Val{3}) = vcat(STANDARD_TAGS, ["lateral_walls"])

# "surface" is a 3D alias for "free_surface" (backwards-compatibility shim)
_normalize_cartesian_tag(::Val{2}, name::String) = name
_normalize_cartesian_tag(::Val{3}, name::String) =
  name == "surface" ? "free_surface" : name

"""
    _face_mask(::Val{D}, mins, maxs, name) -> Function

Return a centroid-based predicate selecting boundary faces named `name` on a
`CartesianDomain{D}` with bounds `mins`/`maxs`.
"""
function _face_mask(
  ::Val{D},
  mins::NTuple{D, Float64},
  maxs::NTuple{D, Float64},
  name::String,
) where {D}
  tol = 1.0e-10
  iz = _vertical_axis(Val(D))
  if name == "free_surface"
    return _centroid_axis_eq_mask(iz, maxs[iz]; tol = tol)
  elseif name == "seabed"
    return _centroid_axis_eq_mask(iz, mins[iz]; tol = tol)
  elseif name == "inlet"
    return _centroid_axis_eq_mask(1, mins[1]; tol = tol)
  elseif name == "outlet"
    return _centroid_axis_eq_mask(1, maxs[1]; tol = tol)
  elseif name == "structure"
    return _always_false_mask()
  elseif D == 3 && name == "lateral_walls"
    return _centroid_axis_either_eq_mask(2, mins[2], maxs[2]; tol = tol)
  else
    error("Unknown boundary tag \"$name\" for CartesianDomain{$D}")
  end
end

function _face_mask(::Val{D}, d::CartesianDomain{D}, name::String) where {D}
  _face_mask(Val(D), d.mins, d.maxs, name)
end

"""
    _empty_boundary(Γ) -> Triangulation

Return an empty sub-triangulation of `Γ` (zero cells).
Used as a placeholder when no cells match a predicate.
"""
function _empty_boundary(Γ)
  Triangulation(Γ, Int[])
end

"""
    _cartesian_boundary_from_model(::Val{D}, model, mins, maxs, name) -> Triangulation

Extract a single named boundary triangulation from a pre-built `model`.

`name` must be a valid tag for `CartesianDomain{D}`.  Returns
`Interior(model)` for `"fluid"`, and a coordinate-filtered
`Triangulation` for all surface tags.
"""
function _cartesian_boundary_from_model(
  ::Val{D},
  model,
  mins::NTuple{D, Float64},
  maxs::NTuple{D, Float64},
  name::String,
) where {D}
  normalized_name = _normalize_cartesian_tag(Val(D), name)
  valid = _cartesian_valid_tags(Val(D))
  normalized_name in valid || error(
    "Unknown boundary tag \"$name\" for CartesianDomain{$D}. " *
    "Valid tags: " * join(valid, ", ") *
    (D == 3 ? ", surface" : ""),
  )

  if normalized_name == "fluid"
    return Interior(model)
  end

  Γ = Boundary(model)
  xΓ = get_cell_coordinates(Γ)
  bits = lazy_map(_face_mask(Val(D), mins, maxs, normalized_name), xΓ)
  Triangulation(Γ, findall(bits))
end

"""
    get_boundary(d::CartesianDomain{D}, name::String) -> Triangulation

Return the boundary triangulation for `name` on `d`.

Supported names: all six [`STANDARD_TAGS`](@ref) plus `"lateral_walls"`
for 3D.  For `"fluid"` the interior triangulation is returned.  For
`"structure"` an empty triangulation is returned (no embedded structures on
a plain `CartesianDomain`; see [`TankDomain`](@ref)).

Note: builds a fresh model each call.
"""
function get_boundary(d::CartesianDomain{D}, name::String) where {D}
  _cartesian_boundary_from_model(Val(D), build_model(d), d.mins, d.maxs, name)
end

# ─────────────────────────────────────────────────────────────────────────────
# Full triangulation set
# ─────────────────────────────────────────────────────────────────────────────

"""
    build_triangulations(d::CartesianDomain{D}, model) -> TankTriangulations

Build the full set of triangulations for a plain Cartesian domain from a
pre-built `model`.

No structure or damping sub-domains are created; `:Γ_structures`,
`:Γ_dampings`, `:Λ_joints`, and `:joint_domains` are empty.  For
`CartesianDomain{3}` an additional `:Γlateral` key is populated.

This is the generic Cartesian fallback.  Use [`TankDomain`](@ref) when
embedded structure or damping sub-domains are needed.
"""
function build_triangulations(d::CartesianDomain{D}, model) where {D}
  _cartesian_boundary_triangulations(Val(D), model, d.mins, d.maxs)
end

function _cartesian_boundary_triangulations(
  ::Val{D},
  model,
  mins::NTuple{D, Float64},
  maxs::NTuple{D, Float64},
) where {D}
  Ω   = Interior(model)
  Γ   = Boundary(model)
  Γfs  = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "free_surface")
  Γbot = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "seabed")
  Γin  = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "inlet")
  Γout = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "outlet")
  Γη   = _empty_boundary(Γ)

  data = _tank_triangulation_dict(
    Ω            = Ω,
    Γ            = Γ,
    Γbot         = Γbot,
    Γin          = Γin,
    Γout         = Γout,
    Γ_structures = Any[],
    Γ_dampings   = Any[],
    Γfs          = Γfs,
    Γκ           = Γfs,
    Γη           = Γη,
    Λη           = nothing,
    Λ_joints     = Any[],
    joint_domains = Any[],
  )

  if D == 3
    data[:Γlateral] = _cartesian_boundary_from_model(
      Val(D), model, mins, maxs, "lateral_walls",
    )
  end

  TankTriangulations(data)
end

# ─────────────────────────────────────────────────────────────────────────────
# 3D coordinate grading helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    f_z(x, H, nz, grading_base=2.5) -> Float64

Vertical grading map: maps a uniform coordinate `x ∈ [0, H]` to a
logarithmically graded coordinate with smaller cells near `z = H` (the
free surface) and larger cells near `z = 0` (seabed).

Used internally by [`map_fn`](@ref).

# Arguments
- `x`            — input coordinate in `[0, H]` [m]
- `H`            — total water depth [m]
- `nz`           — number of layers in the vertical direction
- `grading_base` — geometric ratio controlling cell size ratio (default 2.5)
"""
function f_z(x, H, nz, grading_base = 2.5)
  x == H && return H
  i = x / (H / nz)
  return H - H / (grading_base^i)
end

"""
    map_fn(x, H, nz; grading_base=2.5) -> VectorValue

Coordinate map for graded 3D Cartesian meshes.  Passes `x` and `y`
through unchanged and applies [`f_z`](@ref) to the vertical `z` coordinate.

Designed to be passed to `CartesianDomain(map=...)`.

# Example

```julia
domain = CartesianDomain(
    mins=(0.0, -1.0, 0.0), maxs=(10.0, 1.0, 2.0), parts=(50, 10, 8),
    map = x -> map_fn(x, 2.0, 8; grading_base=3.0),
)
```
"""
map_fn(x, H, nz; grading_base = 2.5) =
  VectorValue(x[1], x[2], f_z(x[3], H, nz, grading_base))
