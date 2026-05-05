"""
  CartesianDomain{D} <: AbstractDomain

Generic axis-aligned Cartesian domain for 2D/3D fluid problems.

`CartesianDomain{2}` uses coordinates `(x, z)`.
`CartesianDomain{3}` uses coordinates `(x, y, z)`.

This type is intentionally functional and dispatch-driven:
- behavior by dimension is selected with `Val(D)` dispatch
- the same API is shared by both 2D and 3D variants

# Fields
- `mins::NTuple{D,Float64}`: lower bounds for each axis [m]
- `maxs::NTuple{D,Float64}`: upper bounds for each axis [m]
- `parts::NTuple{D,Int}`: cell partitions by axis
- `map::Function`: optional coordinate map (default identity)
"""
struct CartesianDomain{D, F <: Function} <: AbstractDomain
  mins::NTuple{D, Float64}
  maxs::NTuple{D, Float64}
  parts::NTuple{D, Int}
  map::F
end

"""
  CartesianDomain(; L, H, nx, ny, W=nothing, nz=nothing, map=x->x)

Build a generic Cartesian domain:
- 2D when `W` and `nz` are not provided
- 3D when both `W` and `nz` are provided
"""
function CartesianDomain(; L, H, nx, ny, W = nothing, nz = nothing, map = x -> x)
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
    "CartesianDomain expects either 2D args (L,H,nx,ny) or " *
    "3D args (L,W,H,nx,ny,nz).",
  )
end

_flatten_bounds(mins::NTuple{2, Float64}, maxs::NTuple{2, Float64}) =
  (mins[1], maxs[1], mins[2], maxs[2])

_flatten_bounds(mins::NTuple{3, Float64}, maxs::NTuple{3, Float64}) =
  (mins[1], maxs[1], mins[2], maxs[2], mins[3], maxs[3])

function build_model(d::CartesianDomain{D}) where {D}
  bounds = _flatten_bounds(d.mins, d.maxs)
  CartesianDiscreteModel(bounds, d.parts; map = d.map)
end

ambient_dimension(::CartesianDomain{D}) where {D} = D
manifold_dimension(::CartesianDomain{D}) where {D} = D

function boundary_tags(::CartesianDomain{2})
  Dict{String, String}(
    "fluid" => "fluid",
    "free_surface" => "free_surface",
    "seabed" => "seabed",
    "inlet" => "inlet",
    "outlet" => "outlet",
    "structure" => "structure",
  )
end

function boundary_tags(::CartesianDomain{3})
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

triangulation(d::CartesianDomain) = Interior(build_model(d))

_vertical_axis(::Val{2}) = 2
_vertical_axis(::Val{3}) = 3

function _face_mask(
  ::Val{D},
  d::CartesianDomain{D},
  name::String,
) where {D}
  tol = 1.0e-10
  iz = _vertical_axis(Val(D))
  if name == "free_surface"
    return xs -> abs(_centroid(xs)[iz] - d.maxs[iz]) < tol
  elseif name == "seabed"
    return xs -> abs(_centroid(xs)[iz] - d.mins[iz]) < tol
  elseif name == "inlet"
    return xs -> abs(_centroid(xs)[1] - d.mins[1]) < tol
  elseif name == "outlet"
    return xs -> abs(_centroid(xs)[1] - d.maxs[1]) < tol
  elseif name == "structure"
    return xs -> false
  elseif D == 3 && name == "lateral_walls"
    return xs -> begin
      c = _centroid(xs)
      abs(c[2] - d.mins[2]) < tol || abs(c[2] - d.maxs[2]) < tol
    end
  else
    error("Unknown boundary tag \"$name\" for CartesianDomain{$D}")
  end
end

function _valid_tags(::Val{2})
  STANDARD_TAGS
end

function _valid_tags(::Val{3})
  vcat(STANDARD_TAGS, ["lateral_walls"])
end

function get_boundary(d::CartesianDomain{D}, name::String) where {D}
  valid = _valid_tags(Val(D))
  name in valid || error(
    "Unknown boundary tag \"$name\" for CartesianDomain{$D}. " *
    "Valid tags: " * join(valid, ", "),
  )

  model = build_model(d)
  if name == "fluid"
    return Interior(model)
  end

  Γ = Boundary(model)
  xΓ = get_cell_coordinates(Γ)
  bits = lazy_map(_face_mask(Val(D), d, name), xΓ)
  Triangulation(Γ, findall(bits))
end

"""
    f_z(x, H, nz, grading_base=2.5)

Vertical grading map used by the Yago 3D benchmark.
Maps a uniform coordinate `x ∈ [0, H]` to a graded coordinate with
smaller cells near `z = H` and larger cells near `z = 0`.
"""
function f_z(x, H, nz, grading_base=2.5)
  x == H && return H
  i = x / (H / nz)
  return H - H / (grading_base^i)
end

"""
    map_fn(x, H, nz; grading_base=2.5)

Coordinate map for graded 3D Cartesian meshes.
"""
map_fn(x, H, nz; grading_base=2.5) =
  VectorValue(x[1], x[2], f_z(x[3], H, nz, grading_base))

"""
    verify_cartesian_3d_tags(model)

Print all raw Gridap Cartesian face tags and verify the expected tag ids
for the top surface and inlet used by the Yago benchmark setup.

Returns `(surface_tag, inlet_tag)`.
"""
function verify_cartesian_3d_tags(model)
  labels = get_face_labeling(model)
  println("[CartesianDomain3D] Raw Gridap face tags:")
  for (i, name) in enumerate(labels.tag_to_name)
    println("  tag $i -> $name")
  end

  @assert "interior" in labels.tag_to_name

  surface_tag = 22
  inlet_tag = 25
  @assert surface_tag <= length(labels.tag_to_name)
  @assert inlet_tag <= length(labels.tag_to_name)

  return surface_tag, inlet_tag
end

"""
    CartesianDomain3D <: AbstractDomain

3D Cartesian fluid domain with geometrically graded vertical mesh.

Physical labels are created for at least `"surface"` and `"inlet"` using
Gridap internal Cartesian tags. Additional standard boundaries are extracted
by coordinate masks in `get_boundary`.
"""
struct CartesianDomain3D <: AbstractDomain
  LΩ::Float64
  BΩ::Float64
  H::Float64
  nx_total::Int
  ny_total::Int
  nz::Int
  grading_base::Float64
  model
  tags::Dict{String, Any}
end

"""
    CartesianDomain3D(; LΩ, BΩ, H, nx_total, ny_total, nz, grading_base=2.5)

Build a graded 3D Cartesian fluid domain for the Yago benchmark.
"""
function CartesianDomain3D(; LΩ, BΩ, H, nx_total, ny_total, nz,
                           grading_base=2.5)
  dom = (0.0, LΩ, -BΩ / 2, BΩ / 2, 0.0, H)
  part = (nx_total, ny_total, nz)
  map = x -> map_fn(x, H, nz; grading_base=grading_base)
  model = CartesianDiscreteModel(dom, part; map=map)

  surface_tag, inlet_tag = verify_cartesian_3d_tags(model)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "surface", [surface_tag])
  add_tag_from_tags!(labels, "inlet", [inlet_tag])

  tags = Dict{String, Any}(
    "fluid" => "fluid",
    "surface" => "surface",
    "free_surface" => "surface",
    "inlet" => "inlet",
    "seabed" => "seabed",
    "outlet" => "outlet",
    "structure" => "structure",
  )

  return CartesianDomain3D(
    Float64(LΩ),
    Float64(BΩ),
    Float64(H),
    Int(nx_total),
    Int(ny_total),
    Int(nz),
    Float64(grading_base),
    model,
    tags,
  )
end

ambient_dimension(::CartesianDomain3D) = 3
manifold_dimension(::CartesianDomain3D) = 3
triangulation(d::CartesianDomain3D) = Interior(d.model)
boundary_tags(d::CartesianDomain3D) = d.tags

function _surface_mask_3d(d::CartesianDomain3D, name::String)
  tol = 1.0e-10
  if name == "surface" || name == "free_surface"
    return xs -> begin
      c = _centroid(xs)
      abs(c[3] - d.H) < tol
    end
  elseif name == "seabed"
    return xs -> begin
      c = _centroid(xs)
      abs(c[3] - 0.0) < tol
    end
  elseif name == "inlet"
    return xs -> begin
      c = _centroid(xs)
      abs(c[1] - 0.0) < tol
    end
  elseif name == "outlet"
    return xs -> begin
      c = _centroid(xs)
      abs(c[1] - d.LΩ) < tol
    end
  elseif name == "structure"
    return xs -> false
  else
    error("Unknown boundary name \"$name\" for CartesianDomain3D")
  end
end

"""
    get_boundary(d::CartesianDomain3D, name::String)

Return interior or boundary triangulation for the requested tag name.
"""
function get_boundary(d::CartesianDomain3D, name::String)
  if name == "fluid"
    return Interior(d.model)
  end

  valid = ["surface", "free_surface", "inlet", "outlet", "seabed", "structure"]
  name in valid || error(
    "Unknown boundary tag \"$name\" for CartesianDomain3D. Valid tags: " *
    join(valid, ", "),
  )

  Γ = Boundary(d.model)
  xΓ = get_cell_coordinates(Γ)
  bits = lazy_map(_surface_mask_3d(d, name), xΓ)
  return Triangulation(Γ, findall(bits))
end
