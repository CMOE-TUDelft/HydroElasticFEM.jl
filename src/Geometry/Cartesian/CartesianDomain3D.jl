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

# Fields
- `LΩ::Float64`      : domain length in x [m]
- `BΩ::Float64`      : domain width in y [m]
- `H::Float64`       : domain depth in z [m]
- `nx_total::Int`    : total cells in x
- `ny_total::Int`    : total cells in y
- `nz::Int`          : total cells in z
- `grading_base`     : geometric grading base in z
- `model`            : Gridap Cartesian discrete model
- `tags`             : mapping for user-facing boundary tags
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
