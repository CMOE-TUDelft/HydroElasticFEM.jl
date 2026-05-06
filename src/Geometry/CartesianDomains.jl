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

boundary_tags(::CartesianDomain{D}) where {D} = _cartesian_boundary_tags(Val(D))

triangulation(d::CartesianDomain) = Interior(build_model(d))

_vertical_axis(::Val{2}) = 2
_vertical_axis(::Val{3}) = 3

_cartesian_valid_tags(::Val{2}) = STANDARD_TAGS
_cartesian_valid_tags(::Val{3}) = vcat(STANDARD_TAGS, ["lateral_walls"])

_normalize_cartesian_tag(::Val{2}, name::String) = name
_normalize_cartesian_tag(::Val{3}, name::String) =
  name == "surface" ? "free_surface" : name

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

function get_boundary(d::CartesianDomain{D}, name::String) where {D}
  _cartesian_boundary_from_model(Val(D), build_model(d), d.mins, d.maxs, name)
end

function _empty_boundary(Γ)
  Triangulation(Γ, Int[])
end

function _cartesian_boundary_triangulations(
  ::Val{D},
  model,
  mins::NTuple{D, Float64},
  maxs::NTuple{D, Float64},
) where {D}
  Ω = Interior(model)
  Γ = Boundary(model)

  Γfs = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "free_surface")
  Γbot = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "seabed")
  Γin = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "inlet")
  Γout = _cartesian_boundary_from_model(Val(D), model, mins, maxs, "outlet")
  Γη = _empty_boundary(Γ)

  data = _tank_triangulation_dict(
    Ω = Ω,
    Γ = Γ,
    Γbot = Γbot,
    Γin = Γin,
    Γout = Γout,
    Γ_structures = Any[],
    Γ_dampings = Any[],
    Γfs = Γfs,
    Γκ = Γfs,
    Γη = Γη,
    Λη = nothing,
    Λ_joints = Any[],
    joint_domains = JointDomain[],
  )

  if D == 3
    data[:Γlateral] = _cartesian_boundary_from_model(
      Val(D),
      model,
      mins,
      maxs,
      "lateral_walls",
    )
  end

  TankTriangulations(data)
end

function build_triangulations(d::CartesianDomain{D}, model) where {D}
  _cartesian_boundary_triangulations(Val(D), model, d.mins, d.maxs)
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
  map = x -> map_fn(x, H, nz; grading_base=grading_base)
  domain = CartesianDomain{3, typeof(map)}(
    (0.0, -Float64(BΩ) / 2, 0.0),
    (Float64(LΩ), Float64(BΩ) / 2, Float64(H)),
    (Int(nx_total), Int(ny_total), Int(nz)),
    map,
  )
  model = build_model(domain)

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
    "lateral_walls" => "lateral_walls",
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
build_model(d::CartesianDomain3D) = d.model
triangulation(d::CartesianDomain3D) = Interior(d.model)
boundary_tags(d::CartesianDomain3D) = d.tags

"""
    get_boundary(d::CartesianDomain3D, name::String)

Return interior or boundary triangulation for the requested tag name.
"""
function get_boundary(d::CartesianDomain3D, name::String)
  mins = (0.0, -d.BΩ / 2, 0.0)
  maxs = (d.LΩ, d.BΩ / 2, d.H)
  _cartesian_boundary_from_model(Val(3), d.model, mins, maxs, name)
end

function build_triangulations(d::CartesianDomain3D, model)
  mins = (0.0, -d.BΩ / 2, 0.0)
  maxs = (d.LΩ, d.BΩ / 2, d.H)
  _cartesian_boundary_triangulations(Val(3), model, mins, maxs)
end
