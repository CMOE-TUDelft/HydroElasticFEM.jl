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

  error("CartesianDomain expects either 2D args (L,H,nx,ny) or " *
        "3D args (L,W,H,nx,ny,nz).")
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
