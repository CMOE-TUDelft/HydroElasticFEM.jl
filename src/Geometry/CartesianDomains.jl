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
  CartesianDomain(; L, H, nx, ny, W=nothing, nz=nothing,
                    mins=nothing, maxs=nothing, parts=nothing, map=x->x)

Build a generic Cartesian domain.

- Shorthand 2D: provide `L, H, nx, ny`
- Shorthand 3D: provide `L, W, H, nx, ny, nz`
- Explicit bounds: provide `mins, maxs, parts`

The explicit form is intended for shifted, centered, or graded structured
meshes that do not fit the box-at-origin shorthand.
"""
function CartesianDomain(; L = nothing, H = nothing, nx = nothing, ny = nothing,
                         W = nothing, nz = nothing, mins = nothing,
                         maxs = nothing, parts = nothing, map = x -> x)
  if !isnothing(mins) || !isnothing(maxs) || !isnothing(parts)
    isnothing(mins) && error("`mins` must be provided with explicit CartesianDomain construction.")
    isnothing(maxs) && error("`maxs` must be provided with explicit CartesianDomain construction.")
    isnothing(parts) && error("`parts` must be provided with explicit CartesianDomain construction.")

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
    "CartesianDomain expects either 2D args (L,H,nx,ny) or " *
    "3D args (L,W,H,nx,ny,nz), or explicit (mins,maxs,parts).",
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

