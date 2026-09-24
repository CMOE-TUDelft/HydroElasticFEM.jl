# ─────────────────────────────────────────────────────────────────────────────
# TankDomain — structured Cartesian tank geometry with embedded sub-domains
# ─────────────────────────────────────────────────────────────────────────────

"""
    AbstractSurfaceZone

Supertype of axis-aligned sub-domain descriptors that partition the tank's
top (free) surface: [`StructureDomain`](@ref) and [`DampingZone`](@ref).

A surface zone is a codimension-1 box on the free surface:

- 2D tank (`ambient_dim = 2`): a segment `x ∈ [x₀[1], x₀[1] + L]` at `y = x₀[2]`.
- 3D tank (`ambient_dim = 3`): a rectangle `x ∈ [x₀[1], x₀[1] + L]`,
  `y ∈ [x₀[2], x₀[2] + W]` at `z = x₀[3]`.

The zone carries geometry only.  Whether it hosts a beam, membrane or plate,
or which damping law acts on it, is decided by the physics layer.
"""
abstract type AbstractSurfaceZone end

"""
  StructureDomain <: AbstractSurfaceZone

Descriptor for an axis-aligned structure sub-domain on the Cartesian free
surface (a floating beam or membrane in 2D; a plate or membrane in 3D).

At `build_triangulations` time, a centroid-based mask selects the surface
cells whose centroid lies inside this descriptor's box.  The selected cells
become the per-structure triangulation (stored under `domain_symbol` in
`TankTriangulations`) and contribute to `Γη`.  The structure's interior
skeleton, its boundary (`∂Γs`: beam end points in 2D, plate edge in 3D) and a
model face tag `boundary_tag` on that boundary are generated automatically.
Use `boundary_tag` as `dirichlet_tags` to clamp or simply support the
structure edge.

## Fields
- `L::Float64`               — extent along x [m]
- `x₀::Vector{Float64}`      — lower corner [m]; its length sets the ambient dimension
- `ambient_dim::Int`         — ambient dimension (default `length(x₀)`)
- `manifold_dim::Int`        — structure manifold dimension (default `ambient_dim - 1`)
- `domain_symbol::Symbol`    — key used in `TankTriangulations` (default `:Γ_s`)
- `W::Union{Nothing,Float64}` — extent along y [m]; required in 3D
- `boundary_tag::String`     — model face tag on the structure boundary
  (default `"\$(domain_symbol)_boundary"`)

## Example
```julia
beam  = StructureDomain(L=1.0, x₀=[1.5, 1.0])                                  # 2D
plate = StructureDomain(L=300.0, W=60.0, x₀=[1350.0, -30.0, 58.5], domain_symbol=:Γ_p)  # 3D
```
"""
@with_kw struct StructureDomain <: AbstractSurfaceZone
  L::Float64 = 1.0
  x₀::Vector{Float64} = [0.0, 1.0]
  ambient_dim::Int = length(x₀)
  manifold_dim::Int = ambient_dim - 1
  domain_symbol::Symbol = :Γ_s
  W::Union{Nothing, Float64} = nothing
  boundary_tag::String = string(domain_symbol, "_boundary")
end

"""
  DampingZone <: AbstractSurfaceZone

Descriptor for an axis-aligned numerical wave absorber on the Cartesian free
surface.

At `build_triangulations` time, a centroid-based mask selects the surface
cells inside this zone.  Damping-zone cells are excluded from the open
free-surface triangulation `Γfs` but included in `Γκ` (the wave-equation
free-surface that includes damping), so the damping operator can apply a
spatially varying coefficient.  The damping law itself (`μ₁`, `μ₂`, …) is
defined by the physics layer (`DampingZoneBC`).

## Fields
- `L::Float64`               — extent along x [m]
- `σ::Float64`               — deprecated, unused by the geometry layer; the
  damping law lives in `DampingZoneBC`
- `x₀::Vector{Float64}`      — lower corner [m]; its length sets the ambient dimension
- `ambient_dim::Int`         — ambient dimension (default `length(x₀)`)
- `manifold_dim::Int`        — manifold dimension (default `ambient_dim - 1`)
- `domain_symbol::Symbol`    — key used in `TankTriangulations` (default `:Γ_d`)
- `W::Union{Nothing,Float64}` — extent along y [m] (3D only); `nothing`
  (default) spans the full tank width
"""
@with_kw struct DampingZone <: AbstractSurfaceZone
  L::Float64 = 0.5
  σ::Float64 = 1.0
  x₀::Vector{Float64} = [0.0, 1.0]
  ambient_dim::Int = length(x₀)
  manifold_dim::Int = ambient_dim - 1
  domain_symbol::Symbol = :Γ_d
  W::Union{Nothing, Float64} = nothing
end

"""
    JointDomain

Declares a rotational-spring joint located at a specific point on a 1D beam
boundary embedded in the 2D tank mesh.

`JointDomain` is purely a geometry-level descriptor.  After you add it to
`TankDomain{2}.joint_domains`, `build_triangulations` automatically builds the
corresponding Gridap skeleton sub-triangulation (one interior facet of
`Skeleton(Γη)` whose centroid is closest to `location`).  Then
`get_integration_domains` populates `IntegrationDomains` with the resulting
measure and outward normal, ready to be consumed by
[`HydroElasticFEM.Physics.EulerBernoulliBeam`](@ref) via the matching
`JointRotationalSpring`.

# Fields
- `location::Vector{Float64}` — 2D coordinates `[x, y]` of the joint point
  on the structure surface (e.g. `[2.0, 1.0]`).
- `domain_symbol::Symbol` — Key under which the Gridap skeleton `Measure` is
  stored in `IntegrationDomains` (e.g. `:dΛj_1`).  Must match
  `JointRotationalSpring.domain_symbol` in the corresponding beam entity.
- `normal_symbol::Symbol` — Key under which the skeleton outward-normal field
  is stored (e.g. `:n_Λ_j_1`).  Must match
  `JointRotationalSpring.normal_symbol` in the corresponding beam entity.
- `tol::Float64` — Spatial tolerance for selecting the skeleton facet whose
  centroid is within `tol` of `location`.  Default `1e-6` m.

# Example

```julia
# Beam spanning x ∈ [1.5, 2.5] at y = 1.0 — joint at midpoint x = 2.0
s1 = StructureDomain(L=1.0, x₀=[1.5, 1.0])
j1 = JointDomain(location=[2.0, 1.0], domain_symbol=:dΛj_1, normal_symbol=:n_Λ_j_1)

tank = TankDomain(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains=[s1],
    joint_domains=[j1])
model  = build_model(tank)
trians = build_triangulations(tank, model)   # skeleton facet extracted here
dom    = get_integration_domains(trians)      # :dΛj_1 and :n_Λ_j_1 populated

# Physics side — link via matching symbols
beam = EulerBernoulliBeam(L=1.0, mᵨ=0.5, EIᵨ=100.0,
    joints=[JointRotationalSpring(:dΛj_1, :n_Λ_j_1, kᵣ)])
```
"""
@with_kw struct JointDomain
  location::Vector{Float64}
  domain_symbol::Symbol
  normal_symbol::Symbol
  tol::Float64 = 1.0e-6
end

"""
    ResonatorDomain

Declares a point interaction used by lumped resonators attached to a tank
sub-domain.

`ResonatorDomain` is a geometry-level descriptor.  After you add it to
`TankDomain.resonator_domains`, `build_triangulations` stores the
descriptor together with the other tank sub-domain metadata.  Then
`get_integration_domains` automatically builds a Gridap `DiracDelta` on
`trian_symbol` at `location` and stores it in `IntegrationDomains` under
`delta_symbol`.  Multiple resonators sharing the same `delta_symbol` are
grouped in input order, matching the DOF order expected by `ResonatorArray`.

# Fields
- `location::Vector{Float64}` — Coordinates of the point interaction (in meters)
  (e.g. `[2.0, 1.0]` on the structure surface of a 2D tank).
- `trian_symbol::Symbol` — Key in `TankTriangulations` of the triangulation
  that supports the point interaction.  Defaults to `:Γη`, the union of
  structure-surface cells.
  For plain tanks without `:Γη`, pass another existing key such as `:Ω` or
  `:Γfs`.
- `delta_symbol::Symbol` — Key under which the vector of `DiracDelta`
  functionals is stored in `IntegrationDomains`.  Defaults to `:δ_p`, the key
  consumed by resonator weak forms.

# Example

```julia
s1 = StructureDomain(L=1.0, x₀=[1.5, 1.0])
r1 = ResonatorDomain(location=[2.0, 1.0])

tank = TankDomain(L=4.0, H=1.0, nx=40, ny=4,
    structure_domains=[s1],
    resonator_domains=[r1])
model  = build_model(tank)
trians = build_triangulations(tank, model)
dom    = get_integration_domains(trians)      # :δ_p populated
```
"""
@with_kw struct ResonatorDomain
  location::Vector{Float64}
  trian_symbol::Symbol = :Γη
  delta_symbol::Symbol = :δ_p
end

"""
    TankDomain{D} <: AbstractDomain

Structured Cartesian tank domain with embedded sub-domain descriptors.

`TankDomain` wraps a [`CartesianDomain`](@ref) and extends it with four
optional lists of sub-domain descriptors:

- `structure_domains` — beam / plate segments on the free surface
- `damping_zones`     — numerical wave absorbers on the free surface
- `joint_domains`     — rotational-spring joints between beam segments
- `resonator_domains` — point interactions for lumped resonators

At `build_triangulations` time, the top-surface `Boundary` is partitioned
into named sub-triangulations using coordinate masks derived from each
descriptor.  This partition is the bridge between the geometric description
and the Gridap `Measure` objects consumed by weak forms.

The dimension `D` is inferred from the keyword arguments:
- 2D: provide `L, H, nx, ny`
- 3D: provide `L, W, H, nx, ny, nz`

`structure_domains`, `damping_zones` and `resonator_domains` are supported
for `D = 2` and `D = 3`; every surface zone must have `ambient_dim == D`.
`joint_domains` (point joints between beam segments) are 2D only.

For a 3D tank that does not start at the origin (e.g. `y ∈ [-W/2, W/2]` or a
free surface at `z = 0`), build the `CartesianDomain` with explicit
`mins`/`maxs` and pass it to `TankDomain(cartesian; ...)`.

## Periodicity

Pass `is_periodic` to enable periodic boundary conditions along one or more
axes.  It must be an `NTuple{D,Bool}` whose length matches the domain
dimension:

```julia
# 2D domain periodic in x (useful for regular wave trains)
tank = TankDomain(L=4.0, H=1.0, nx=40, ny=8, is_periodic=(true, false))
```

## Compatibility properties

`TankDomain` exposes the wrapped Cartesian dimensions as direct properties for
backwards compatibility: `domain.L`, `domain.H`, `domain.nx`, `domain.ny`,
`domain.is_periodic` (2D) and additionally `domain.W`, `domain.nz` (3D).
These are read-only derived values; the canonical data lives in
`domain.cartesian`.
"""
struct TankDomain{D, C, SZ, DZ, JZ, RZ} <: AbstractDomain
  cartesian::C
  structure_domains::SZ
  damping_zones::DZ
  joint_domains::JZ
  resonator_domains::RZ
end

function _validate_tank_domain_inputs(
  ::Val{D},
  structure_domains,
  damping_zones,
  joint_domains,
  resonator_domains,
) where {D}
  for zone in Iterators.flatten((structure_domains, damping_zones))
    zone.ambient_dim == D || error(
      "$(nameof(typeof(zone))) :$(zone.domain_symbol) has ambient_dim=$(zone.ambient_dim) " *
      "(length(x₀)=$(length(zone.x₀))) but the tank is $(D)D.",
    )
    _zone_bounds(zone)  # validates x₀, dimensions and extents
  end
  D == 2 || isempty(joint_domains) ||
    error("TankDomain{$D} does not yet support joint_domains (2D point joints only).")
  nothing
end

"""
    TankDomain(cartesian; structure_domains=[], damping_zones=[], joint_domains=[], resonator_domains=[])

Construct a `TankDomain` from an existing `CartesianDomain`.

This is the low-level constructor used when you already have a `CartesianDomain`
object.  See also the keyword-argument constructor `TankDomain(; L, H, nx, ny, ...)`
for a one-shot variant that builds the `CartesianDomain` internally.

# Arguments
- `cartesian`: a `CartesianDomain{D}` for D = 2 or 3
- `structure_domains::Vector{StructureDomain}`: structural subregions on the free surface
- `damping_zones::Vector{DampingZone}`: sponge-layer regions on the free surface
- `joint_domains::Vector{JointDomain}`: interior skeleton facets for beam joints (2D only)
- `resonator_domains::Vector{ResonatorDomain}`: point interactions for lumped resonators

# Returns
- `TankDomain{D,...}`: configured tank domain

# Example
```julia
cart = G.CartesianDomain(L=10.0, H=1.0, nx=40, ny=8)
domain = G.TankDomain(cart; structure_domains=[G.StructureDomain(L=4.0, x₀=[3.0, 1.0])])
```
"""
function TankDomain(
  cartesian;
  structure_domains = StructureDomain[],
  damping_zones = DampingZone[],
  joint_domains = JointDomain[],
  resonator_domains = ResonatorDomain[],
)
  D = ambient_dimension(cartesian)
  _validate_tank_domain_inputs(
    Val(D),
    structure_domains,
    damping_zones,
    joint_domains,
    resonator_domains,
  )
  TankDomain{
    D,
    typeof(cartesian),
    typeof(structure_domains),
    typeof(damping_zones),
    typeof(joint_domains),
    typeof(resonator_domains),
  }(
    cartesian,
    structure_domains,
    damping_zones,
    joint_domains,
    resonator_domains,
  )
end

"""
    TankDomain(; L, H, nx, ny, W=nothing, nz=nothing, map=identity,
               is_periodic=nothing, structure_domains=[], damping_zones=[],
               joint_domains=[], resonator_domains=[])

Build a `TankDomain` from scratch using keyword arguments.

Creates the underlying `CartesianDomain` internally and wraps it with the
given structural, damping-zone, joint, and resonator sub-regions.  Use `W`
and `nz` to create a 3D domain; omit them for 2D.

# Arguments
- `L::Real`: domain length in x [m]
- `H::Real`: domain height (depth) in z [m]
- `nx::Int`: number of cells in x
- `ny::Int`: number of cells in y (vertical for 2D)
- `W::Union{Real,Nothing}`: domain width in y [m]; `nothing` for 2D
- `nz::Union{Int,Nothing}`: number of cells in z; `nothing` for 2D
- `map`: optional coordinate mapping function, default `identity`
- `is_periodic`: pass-through to `CartesianDomain`; `nothing` for default behaviour
- `structure_domains::Vector{StructureDomain}`: structural subregions
- `damping_zones::Vector{DampingZone}`: sponge-layer regions
- `joint_domains::Vector{JointDomain}`: interior beam-joint facets (2D only)
- `resonator_domains::Vector{ResonatorDomain}`: point interactions for lumped resonators

# Returns
- `TankDomain{D,...}` with D = 2 (no `W`/`nz`) or D = 3

# Example
```julia
domain = G.TankDomain(L=10.0, H=1.0, nx=60, ny=8)
```
"""
function TankDomain(;
  L = 4.0,
  H = 1.0,
  nx = 16,
  ny = 2,
  W = nothing,
  nz = nothing,
  map = x -> x,
  is_periodic = nothing,
  structure_domains = StructureDomain[],
  damping_zones = DampingZone[],
  joint_domains = JointDomain[],
  resonator_domains = ResonatorDomain[],
)
  cartesian = CartesianDomain(
    L = L,
    H = H,
    nx = nx,
    ny = ny,
    W = W,
    nz = nz,
    map = map,
    is_periodic = is_periodic,
  )
  TankDomain(
    cartesian;
    structure_domains = structure_domains,
    damping_zones = damping_zones,
    joint_domains = joint_domains,
    resonator_domains = resonator_domains,
  )
end

# Retrieve the underlying CartesianDomain stored in the :cartesian field.
function _cartesian_domain(domain::TankDomain)
  getfield(domain, :cartesian)
end

# Extract the named dimension/partition properties of a 2D TankDomain for
# use by legacy code and property accessors.
function _tank_legacy_dimensions(domain::TankDomain{2})
  cartesian = _cartesian_domain(domain)
  (
    L = cartesian.maxs[1] - cartesian.mins[1],
    H = cartesian.maxs[2] - cartesian.mins[2],
    nx = cartesian.parts[1],
    ny = cartesian.parts[2],
    map = cartesian.map,
    is_periodic = cartesian.is_periodic,
  )
end

function _tank_legacy_dimensions(domain::TankDomain{3})
  cartesian = _cartesian_domain(domain)
  (
    L = cartesian.maxs[1] - cartesian.mins[1],
    W = cartesian.maxs[2] - cartesian.mins[2],
    H = cartesian.maxs[3] - cartesian.mins[3],
    nx = cartesian.parts[1],
    ny = cartesian.parts[2],
    nz = cartesian.parts[3],
    map = cartesian.map,
    is_periodic = cartesian.is_periodic,
  )
end

_tank_legacy_property_names(::Val{2}) = (:L, :H, :nx, :ny, :map, :is_periodic)
_tank_legacy_property_names(::Val{3}) = (:L, :W, :H, :nx, :ny, :nz, :map, :is_periodic)

function _tank_public_property_names(::Val{D}) where {D}
  (
    :cartesian,
    :structure_domains,
    :damping_zones,
    :joint_domains,
    :resonator_domains,
    _tank_legacy_property_names(Val(D))...,
  )
end

function Base.getproperty(domain::TankDomain{D}, name::Symbol) where {D}
  if name in _tank_legacy_property_names(Val(D))
    return getproperty(_tank_legacy_dimensions(domain), name)
  end
  getfield(domain, name)
end

function Base.propertynames(::TankDomain{D}, private::Bool = false) where {D}
  _tank_public_property_names(Val(D))
end

"""
    build_model(domain::TankDomain{D})

Build a `CartesianDiscreteModel` from the wrapped `CartesianDomain`.

For `D=2`, the model spans `[0,L]×[0,H]` with `nx×ny` elements.
For `D=3`, the model spans `[0,L]×[0,W]×[0,H]` with `nx×ny×nz` elements.
"""
function build_model(domain::TankDomain{D}) where {D}
  build_model(_cartesian_domain(domain))
end


# ─────────────────────────────────────────────────────────────────────────────
# Surface masking
# ─────────────────────────────────────────────────────────────────────────────

# Absolute tolerance [m] for the off-surface (vertical) centroid coordinate.
# Absolute rather than relative so that a free surface at z = 0 is matched.
const _SURFACE_ZONE_ATOL = 1.0e-8

# In-plane extent of `zone` along manifold axis `i`: `L` along x, `W` along y.
# `nothing` means the zone is unbounded along that axis (full tank width).
_zone_extent(zone, i::Int) = i == 1 ? zone.L : zone.W

_missing_width_error(zone::StructureDomain) = error(
  "StructureDomain :$(zone.domain_symbol) in 3D requires the y-extent `W`.",
)
_missing_width_error(::DampingZone) = nothing  # spans the full tank width

"""
    _zone_bounds(zone::AbstractSurfaceZone) -> (lo, hi)

Return per-manifold-axis bounds `lo[i] ≤ c[i] ≤ hi[i]` of a surface zone,
after validating its dimensions.  An unbounded axis (a `DampingZone` with
`W = nothing` in 3D) gets `(-Inf, Inf)`.
"""
function _zone_bounds(zone::AbstractSurfaceZone)
  _check_surface_zone(zone)
  m = zone.manifold_dim
  lo = Vector{Float64}(undef, m)
  hi = Vector{Float64}(undef, m)
  for i in 1:m
    ext = _zone_extent(zone, i)
    if ext === nothing
      _missing_width_error(zone)
      lo[i], hi[i] = -Inf, Inf
    else
      ext > 0 || error("$(nameof(typeof(zone))) :$(zone.domain_symbol) has non-positive extent $ext along axis $i.")
      lo[i], hi[i] = zone.x₀[i], zone.x₀[i] + ext
    end
  end
  return lo, hi
end

function _surface_zone_mask(zone::AbstractSurfaceZone)
  lo, hi = _zone_bounds(zone)
  x0 = zone.x₀
  m = zone.manifold_dim
  D = zone.ambient_dim
  return function (xs)
    c = _centroid(xs)
    in_plane = all(lo[i] <= c[i] <= hi[i] for i in 1:m)
    on_surface = all(abs(c[i] - x0[i]) <= _SURFACE_ZONE_ATOL for i in (m + 1):D)
    in_plane && on_surface
  end
end

"""
    surface_mask(zone::AbstractSurfaceZone) -> Function

Return a closure `(xs) -> Bool` that tests whether the centroid of a surface
cell lies within `zone` ([`StructureDomain`](@ref) or [`DampingZone`](@ref)).

The centroid must satisfy `x₀[1] ≤ c[1] ≤ x₀[1] + L` and, in 3D,
`x₀[2] ≤ c[2] ≤ x₀[2] + W` (no y-constraint for a `DampingZone` with
`W = nothing`).  The off-surface coordinate must equal `x₀[end]` up to an
absolute tolerance of `1e-8`.
"""
surface_mask(zone::AbstractSurfaceZone) = _surface_zone_mask(zone)

function _check_surface_zone(zone)
  length(zone.x₀) == zone.ambient_dim ||
    error("length(x₀)=$(length(zone.x₀)) must equal ambient_dim=$(zone.ambient_dim).")
  zone.ambient_dim in (2, 3) || error("ambient_dim must be 2 or 3.")
  zone.manifold_dim == zone.ambient_dim - 1 ||
    error("manifold_dim must be ambient_dim - 1 for a free-surface sub-domain.")
  nothing
end

"""
    surface_masks(domain::TankDomain)

Build mask closures for every structure domain and damping zone in `domain`.

Returns `(structure_masks, damping_masks)` where each is a `Vector{Function}`,
ordered to match `domain.structure_domains` and `domain.damping_zones`.
"""
function surface_masks(domain::TankDomain)
  smasks = [surface_mask(s) for s in domain.structure_domains]
  dmasks = [surface_mask(d) for d in domain.damping_zones]
  return smasks, dmasks
end

"""
    joint_mask(joint::JointDomain) -> Function

Return a closure `(xs) -> Bool` that selects skeleton cells whose centroid is
at `joint.location` within the joint tolerance.
"""
function joint_mask(joint::JointDomain)
  loc = joint.location
  tol = joint.tol
  return function (xs)
    c = _centroid(xs)
    all(abs(c[i] - loc[i]) <= tol for i in 1:length(loc))
  end
end


# ─────────────────────────────────────────────────────────────────────────────
# Surface partition helpers
# ─────────────────────────────────────────────────────────────────────────────

function _unique_domain_symbols(zones, kind)
  syms = [z.domain_symbol for z in zones]
  seen = Set{Symbol}()
  for sym in syms
    if sym in seen
      error(
        "Duplicate domain_symbol :$sym found in $kind. Each $(kind == "structure_domains" ? "structure" : "damping zone") must have a unique domain_symbol.",
      )
    end
    push!(seen, sym)
  end
  syms
end

function _partition_surface_zones(Γ, xΓ, zones, kind)
  syms = _unique_domain_symbols(zones, kind)
  bits = [lazy_map(surface_mask(z), xΓ) for z in zones]
  trians = Any[Triangulation(Γ, findall(b)) for b in bits]
  return (triangulations = trians, symbols = syms, bits = bits)
end

function _partition_surface_by_zones(Γ, xΓ, structure_zones, damping_zones)
  structures = _partition_surface_zones(Γ, xΓ, structure_zones, "structure_domains")
  dampings = _partition_surface_zones(Γ, xΓ, damping_zones, "damping_zones")

  n = length(xΓ)
  any_structure = _or_bits(structures.bits, n)
  any_damping = _or_bits(dampings.bits, n)
  any_zone = any_structure .| any_damping

  return (
    structures = structures,
    dampings = dampings,
    Γfs = Triangulation(Γ, findall(!, any_zone)),
    Γκ = Triangulation(Γ, findall(!, any_structure)),
    Γη = Triangulation(Γ, findall(any_structure)),
    any_structure = any_structure,
  )
end

function _partition_joint_skeletons(Γη, joints)
  if isempty(joints)
    return (Λη = nothing, Λ_joints = Any[], symbols = Symbol[])
  end

  Λη = Skeleton(Γη)
  xΛη = get_cell_coordinates(Λη)
  joint_bits_all = falses(length(xΛη))
  non_joint_bits = trues(length(xΛη))
  Λ_joints = Any[]
  joint_domain_syms = Symbol[]
  seen_joint_domain_syms = Set{Symbol}()
  seen_joint_normal_syms = Set{Symbol}()

  for joint in joints
    if joint.domain_symbol in seen_joint_domain_syms
      error(
        "Duplicate joint domain_symbol :$(joint.domain_symbol) found in joint_domains. " *
        "Each joint must have a unique domain_symbol.",
      )
    end
    if joint.normal_symbol in seen_joint_normal_syms
      error(
        "Duplicate joint normal_symbol :$(joint.normal_symbol) found in joint_domains. " *
        "Each joint must have a unique normal_symbol.",
      )
    end
    push!(seen_joint_domain_syms, joint.domain_symbol)
    push!(seen_joint_normal_syms, joint.normal_symbol)

    bits = lazy_map(joint_mask(joint), xΛη)
    joint_idxs = findall(bits)
    if isempty(joint_idxs)
      error(
        "Joint at location $(joint.location) did not match any structure skeleton cell. " *
        "Check location and tolerance.",
      )
    end

    overlap_bits = joint_bits_all .& bits
    if any(overlap_bits)
      overlap_idxs = findall(overlap_bits)
      error(
        "Joint domain_symbol :$(joint.domain_symbol) overlaps previously assigned " *
        "joint skeleton facet(s) at index/indices $(overlap_idxs). " *
        "Each skeleton facet may belong to at most one joint.",
      )
    end

    joint_bits_all .= joint_bits_all .| bits
    non_joint_bits .= non_joint_bits .& .!bits
    push!(Λ_joints, Triangulation(Λη, joint_idxs))
    push!(joint_domain_syms, joint.domain_symbol)
  end

  Λη_no_joints = Triangulation(Λη, findall(non_joint_bits))
  if (num_cells(Λη_no_joints) + count(joint_bits_all)) != length(xΛη)
    error("Joint/non-joint skeleton partition is inconsistent for Γη.")
  end

  return (Λη = Λη_no_joints, Λ_joints = Λ_joints, symbols = joint_domain_syms)
end

function _add_triangulations_by_symbol!(trian_dict, symbols, trians, label)
  for (sym, tri) in zip(symbols, trians)
    if haskey(trian_dict, sym)
      error(
        "domain_symbol :$sym for $label triangulation would overwrite an existing " *
        "triangulation key in TankTriangulations. All keys must be unique.",
      )
    end
    trian_dict[sym] = tri
  end
  trian_dict
end

"""
    _or_bits(bits, n) -> BitVector

Element-wise OR of a vector of boolean arrays.
Returns a `falses(n)` vector when `bits` is empty.
"""
function _or_bits(bits, n)
  result = falses(n)
  for b in bits
    result .= result .| b
  end
  return result
end


# ─────────────────────────────────────────────────────────────────────────────
# Model labelling
# ─────────────────────────────────────────────────────────────────────────────

"""
    _tank_label_map(::Val{D}) -> Dict{String, String}

Internal mapping from [`STANDARD_TAGS`](@ref) (plus `"lateral_walls"` in 3D)
to the Gridap face-label strings added by [`_label_tank_model!`](@ref).

`"fluid"` maps to the interior label `"water"`.  `"structure"` is handled
separately via coordinate masks in [`get_boundary`](@ref).
"""
function _tank_label_map(::Val{2})
  Dict{String, String}(
    "fluid"        => "water",
    "free_surface" => "surface",
    "seabed"       => "bottom",
    "inlet"        => "inlet",
    "outlet"       => "outlet",
    "structure"    => "structure",  # virtual; resolved via masks
  )
end

function _tank_label_map(::Val{3})
  lmap = _tank_label_map(Val(2))
  lmap["lateral_walls"] = "lateral_walls"
  lmap
end

# Gridap entity ids of a `CartesianDiscreteModel`, i.e. the face ids of the
# reference n-cube.  Entity ids are shared across dimensions: the vertices and
# edges interior to a boundary face carry that face's id.
#
# 2D (QUAD):  vertices 1–4, edges 5 (y=min), 6 (y=max), 7 (x=min), 8 (x=max),
#             interior 9.
# 3D (HEX):   vertices 1–8, edges 9–20, faces 21 (z=min), 22 (z=max),
#             23 (y=min), 24 (y=max), 25 (x=min), 26 (x=max), interior 27.
#
# "surface"/"bottom" include the closure of the top/bottom face.  Inlet,
# outlet and lateral walls include face interiors only (as in 2D, the shared
# corners belong to the surface/bottom).
_tank_entity_ids(::Val{2}) = (
  surface = [3, 4, 6],
  bottom  = [1, 2, 5],
  inlet   = [7],
  outlet  = [8],
  water   = [9],
)

_tank_entity_ids(::Val{3}) = (
  surface       = [5, 6, 7, 8, 11, 12, 15, 16, 22],
  bottom        = [1, 2, 3, 4, 9, 10, 13, 14, 21],
  inlet         = [25],
  outlet        = [26],
  lateral_walls = [23, 24],
  water         = [27],
)

"""
    _label_tank_model!(model) -> nothing

Apply the standard Cartesian face labels to `model` in-place.

Adds `"surface"`, `"bottom"`, `"inlet"`, `"outlet"`, `"water"` (and
`"lateral_walls"` in 3D) tags using the fixed Gridap entity-id convention of
Cartesian meshes (see `_tank_entity_ids`).  Tags are assigned by entity id,
so they are independent of the tank's coordinates (e.g. a free surface at
`z = 0`).  Calling it twice on the same model is a no-op.
"""
function _label_tank_model!(model)
  labels = get_face_labeling(model)
  "water" in labels.tag_to_name && return nothing
  for (name, ids) in pairs(_tank_entity_ids(Val(num_cell_dims(model))))
    add_tag_from_tags!(labels, String(name), ids)
  end
  nothing
end

# Build a fresh model with the standard tank labels applied.
function _labelled_tank_model(domain::TankDomain)
  model = build_model(domain)
  _label_tank_model!(model)
  model
end


# ─────────────────────────────────────────────────────────────────────────────
# Structure skeletons and boundaries
# ─────────────────────────────────────────────────────────────────────────────

function _check_nonempty_structures(structures, trians)
  for (s, Γs) in zip(structures, trians)
    num_cells(Γs) > 0 || error(
      "StructureDomain :$(s.domain_symbol) does not contain any free-surface cell. " *
      "Check x₀, L (and W in 3D) against the tank dimensions and mesh.",
    )
  end
  nothing
end

# Interior skeleton of one structure, excluding facets owned by joints.
function _structure_skeleton(Γs, joints)
  Λ = Skeleton(Γs)
  isempty(joints) && return Λ
  masks = [joint_mask(j) for j in joints]
  keep = findall(xs -> !any(m(xs) for m in masks), collect(get_cell_coordinates(Λ)))
  Triangulation(Λ, keep)
end

"""
    _structure_boundary_faces(model, Γs) -> Dict{Int, Vector{Int}}

Return the model faces on the boundary `∂Γs` of a free-surface
sub-triangulation `Γs`, grouped by face dimension.

A `(D-2)`-face (edge in 3D, vertex in 2D) is on `∂Γs` when it belongs to
exactly one cell of `Γs`.  The lower-dimensional closure (the edge end
points in 3D) is included so that `dirichlet_tags` constrain every DOF on
the boundary.
"""
function _structure_boundary_faces(model, Γs)
  Dc = num_cell_dims(model)
  topo = Gridap.Geometry.get_grid_topology(model)
  mfaces = Gridap.Geometry.get_glue(Γs, Val(Dc - 1)).tface_to_mface
  face_to_bfaces = Gridap.Geometry.get_faces(topo, Dc - 1, Dc - 2)
  count = Dict{Int, Int}()
  for f in mfaces, b in face_to_bfaces[f]
    count[b] = get(count, b, 0) + 1
  end
  bfaces = sort!([b for (b, c) in count if c == 1])
  faces_by_dim = Dict{Int, Vector{Int}}(Dc - 2 => bfaces)
  for d in 0:(Dc - 3)
    bface_to_dfaces = Gridap.Geometry.get_faces(topo, Dc - 2, d)
    faces_by_dim[d] = sort!(unique(reduce(vcat, (collect(bface_to_dfaces[b]) for b in bfaces); init = Int[])))
  end
  faces_by_dim
end

"""
    _relabel_faces!(labels, faces_by_dim) -> Vector{Int}

Move the given model faces to fresh entities and return the new entity ids.

For every original entity `o` touched, a new entity `n` is created and added
to every existing tag that contains `o`, so existing tags (`"surface"`,
`"boundary"`, …) keep covering exactly the same faces.
"""
function _relabel_faces!(labels, faces_by_dim)
  next = maximum(maximum(ents; init = 0) for ents in labels.d_to_dface_to_entity)
  new_entities = Int[]
  for (d, faces) in faces_by_dim
    ents = labels.d_to_dface_to_entity[d + 1]
    for old in unique(ents[faces])
      next += 1
      for f in faces
        ents[f] == old && (ents[f] = next)
      end
      for tag_entities in labels.tag_to_entities
        old in tag_entities && push!(tag_entities, next)
      end
      push!(new_entities, next)
    end
  end
  new_entities
end

"""
    _tag_structure_boundaries!(model, structures, trians) -> nothing

Add a model face tag `structure.boundary_tag` on the boundary of every
structure, plus the union tag `"structure_boundary"`.  These tags can be
passed as `dirichlet_tags` to FE spaces on the structure triangulations.
Calling it twice on the same model is a no-op.
"""
function _tag_structure_boundaries!(model, structures, trians)
  isempty(structures) && return nothing
  labels = get_face_labeling(model)
  "structure_boundary" in labels.tag_to_name && return nothing
  all_new = Int[]
  for (s, Γs) in zip(structures, trians)
    s.boundary_tag in labels.tag_to_name && error(
      "StructureDomain :$(s.domain_symbol) boundary_tag \"$(s.boundary_tag)\" " *
      "already exists in the model face labeling. Each boundary_tag must be unique.",
    )
    new_entities = _relabel_faces!(labels, _structure_boundary_faces(model, Γs))
    Gridap.Geometry.add_tag!(labels, s.boundary_tag, new_entities)
    append!(all_new, new_entities)
  end
  Gridap.Geometry.add_tag!(labels, "structure_boundary", unique(all_new))
  nothing
end


# ─────────────────────────────────────────────────────────────────────────────
# TankDomain{D} triangulations
# ─────────────────────────────────────────────────────────────────────────────

"""
    build_triangulations(domain::TankDomain{D}, model) -> TankTriangulations

Build the full set of named sub-triangulations from `model` (for `D = 2`
and `D = 3`).  The standard tank labels are applied to `model` in place.

## Surface partition

The top-surface `Boundary` (`:Γ`) is partitioned into three regions by
applying centroid masks derived from `domain.structure_domains` and
`domain.damping_zones`:

- `Γη`  — cells covered by at least one `StructureDomain`
- `Γκ`  — cells *not* covered by any `StructureDomain` (open water + damping)
- `Γfs` — cells covered by neither structure nor damping (pure open water)

Each descriptor also produces its own named triangulation stored under its
`domain_symbol` key.

## Structure skeletons and boundaries

- `:Λη` — interior skeleton of `Γη` (for C/DG beam/plate terms).  If
  `joint_domains` are present (2D only), joint facets are split off into
  `:Λ_joints`.
- `:Λ_structures[i]` — interior skeleton of structure `i` (joint facets excluded).
- `:∂Γ_structures[i]` — boundary of structure `i`: end points in 2D, edge
  curve in 3D.  `:∂Γη` is the boundary of the union `Γη`.
- Model face tags `structure.boundary_tag` and `"structure_boundary"` on
  those boundaries, for use as `dirichlet_tags`.

`domain.resonator_domains` are carried as metadata so
`get_integration_domains` can build grouped `DiracDelta` functionals under
`:δ_p`.

In 3D, `:Γlateral` holds the lateral walls (`y = min` and `y = max`).
"""
function build_triangulations(domain::TankDomain{D}, model) where {D}
  # — Label model faces ————————————————————————————————
  _label_tank_model!(model)

  # — Base triangulations ——————————————————————————————
  Ω    = Interior(model)
  Γ    = Boundary(model, tags = "surface")
  Γbot = Boundary(model, tags = "bottom")
  Γin  = Boundary(model, tags = "inlet")
  Γout = Boundary(model, tags = "outlet")

  # — Build boolean masks on surface cell coordinates ——
  xΓ = get_cell_coordinates(Γ)
  surface_partition = _partition_surface_by_zones(
    Γ,
    xΓ,
    domain.structure_domains,
    domain.damping_zones,
  )
  Γfs = surface_partition.Γfs
  Γκ = surface_partition.Γκ
  Γη = surface_partition.Γη
  structure_trians = surface_partition.structures.triangulations
  _check_nonempty_structures(domain.structure_domains, structure_trians)

  if !isempty(domain.joint_domains) && isempty(domain.structure_domains)
    error("Joint domains require at least one structure domain in TankDomain{$D}.")
  end

  # Joint skeleton triangulations from the structure skeleton
  if !isempty(domain.joint_domains)
    joint_partition = _partition_joint_skeletons(Γη, domain.joint_domains)
  elseif !isempty(domain.structure_domains)
    # No joints: use full structure skeleton.
    joint_partition = (Λη = Skeleton(Γη), Λ_joints = Any[], symbols = Symbol[])
  else
    joint_partition = (Λη = nothing, Λ_joints = Any[], symbols = Symbol[])
  end

  # Per-structure skeletons and boundaries, plus boundary tags on the model
  Λ_structures = Any[_structure_skeleton(Γs, domain.joint_domains) for Γs in structure_trians]
  ∂Γ_structures = Any[Boundary(Γs) for Γs in structure_trians]
  ∂Γη = isempty(structure_trians) ? nothing : Boundary(Γη)
  _tag_structure_boundaries!(model, domain.structure_domains, structure_trians)

  # Compose dictionary for TankTriangulations
  trian_dict = _tank_triangulation_dict(
    Ω = Ω,
    Γ = Γ,
    Γbot = Γbot,
    Γin = Γin,
    Γout = Γout,
    Γ_structures = structure_trians,
    Γ_dampings = surface_partition.dampings.triangulations,
    Γfs = Γfs,
    Γκ = Γκ,
    Γη = Γη,
    Λη = joint_partition.Λη,
    Λ_joints = joint_partition.Λ_joints,
    joint_domains = domain.joint_domains,
    resonator_domains = domain.resonator_domains,
    Λ_structures = Λ_structures,
    ∂Γ_structures = ∂Γ_structures,
    ∂Γη = ∂Γη,
  )
  if D == 3
    trian_dict[:Γlateral] = Boundary(model, tags = "lateral_walls")
  end
  _add_triangulations_by_symbol!(
    trian_dict,
    surface_partition.structures.symbols,
    structure_trians,
    "structure",
  )
  _add_triangulations_by_symbol!(
    trian_dict,
    surface_partition.dampings.symbols,
    surface_partition.dampings.triangulations,
    "damping",
  )
  _add_triangulations_by_symbol!(
    trian_dict,
    joint_partition.symbols,
    joint_partition.Λ_joints,
    "joint",
  )
  return TankTriangulations(trian_dict)
end


# ─────────────────────────────────────────────────────────────────────────────
# AbstractDomain interface for TankDomain
# ─────────────────────────────────────────────────────────────────────────────

"""
    ambient_dimension(d::TankDomain{D}) -> Int

Return the ambient dimension `D`.
"""
ambient_dimension(::TankDomain{D}) where {D} = D

"""
    manifold_dimension(d::TankDomain{D}) -> Int

Return the manifold dimension `D` for volume meshes.
"""
manifold_dimension(::TankDomain{D}) where {D} = D

"""
    boundary_tags(d::TankDomain{D}) -> Dict{String, String}

Return a dictionary mapping each [`STANDARD_TAGS`](@ref) name (plus
`"lateral_walls"` in 3D) to its Gridap face-label string on the Cartesian
mesh.

The `"structure"` key maps to `"structure"` (a virtual label resolved via
coordinate masks in [`get_boundary`](@ref) and [`build_triangulations`](@ref)).
"""
boundary_tags(::TankDomain{D}) where {D} = _tank_label_map(Val(D))

"""
    triangulation(d::TankDomain) -> Triangulation

Build the Cartesian discrete model from `d` and return the bulk-fluid
interior triangulation `Ω`.

Note: this method constructs a fresh `CartesianDiscreteModel` every call.
For performance-critical code use [`build_model`](@ref) /
[`build_triangulations`](@ref) directly.
"""
triangulation(d::TankDomain) = Interior(_labelled_tank_model(d))

"""
    get_boundary(d::TankDomain{D}, name::String) -> BoundaryTriangulation

Return the Gridap boundary triangulation for the standard region `name`.

Supported `name` values: all six [`STANDARD_TAGS`](@ref), plus
`"lateral_walls"` (and the alias `"surface"` for `"free_surface"`) in 3D.

For `"fluid"` the interior triangulation is returned.  For `"structure"` the
union of all structure-domain cells on the top surface is returned (empty
if no structure domains are defined).

Note: builds a fresh model on each call; use [`build_triangulations`](@ref)
when multiple boundaries are needed.
"""
function get_boundary(d::TankDomain{D}, name::String) where {D}
  lmap = _tank_label_map(Val(D))
  tag = _normalize_cartesian_tag(Val(D), name)
  haskey(lmap, tag) || error(
    "Unknown boundary tag \"$name\" for TankDomain{$D}. " *
    "Valid tags are: " * join(sort(collect(keys(lmap))), ", ") * ".",
  )
  model = _labelled_tank_model(d)

  if tag == "fluid"
    return Interior(model)
  end

  if tag == "structure"
    # Build the full top-surface triangulation, then sub-select structure cells
    Γ = Boundary(model, tags = lmap["free_surface"])
    if isempty(d.structure_domains)
      return Triangulation(Γ, Int[])
    end
    xΓ = get_cell_coordinates(Γ)
    structures = _partition_surface_zones(Γ, xΓ, d.structure_domains, "structure_domains")
    return Triangulation(Γ, findall(_or_bits(structures.bits, length(xΓ))))
  end

  Boundary(model, tags = lmap[tag])
end


# ─────────────────────────────────────────────────────────────────────────────
# Plate triangulation helper
# ─────────────────────────────────────────────────────────────────────────────

"""
    get_plate_triangulation(Γ, xb₀, xb₁, yb₀, yb₁)

Split a top-surface triangulation `Γ` into plate and free-surface
sub-triangulations by coordinate mask.

A cell belongs to the plate only if all its nodes satisfy:
`xb₀ ≤ x ≤ xb₁` and `yb₀ ≤ y ≤ yb₁`.

Returns `(Γb, Γf, Λb)` where `Λb = Skeleton(Γb)`.

!!! note "Legacy"
    Prefer declaring a 3D [`StructureDomain`](@ref) (`L`, `W`, `x₀`) on a
    [`TankDomain`](@ref): `build_triangulations` then produces `Γη`, `Γκ`,
    `Λη`, the structure boundary and its Dirichlet tag automatically.
"""
function get_plate_triangulation(Γ, xb₀, xb₁, yb₀, yb₁)
  function is_plate(cell_nodes)
    minimum([(xb₀ <= n[1] <= xb₁) && (yb₀ <= n[2] <= yb₁)
             for n in cell_nodes])
  end

  xΓ = get_cell_coordinates(Γ)
  mask = lazy_map(is_plate, xΓ)
  Γb_idx = findall(mask)
  Γf_idx = findall(!, mask)

  Γb = Triangulation(Γ, Γb_idx)
  Γf = Triangulation(Γ, Γf_idx)
  Λb = Skeleton(Γb)

  return Γb, Γf, Λb
end
