"""
  StructureDomain

Defines an axis-aligned structure sub-domain on a Cartesian free surface.

By default this represents a 1D segment embedded in 2D (`ambient_dim=2`,
`manifold_dim=1`), but the descriptor is dimension-aware and can be
parameterized for higher-dimensional embeddings.

Variables:
- L::Float64 = 1.0: Size of the sub-domain along each manifold axis [m]
- x₀::Vector{Float64} = [0.0, 1.0]: Lower corner/anchor point [m]
- ambient_dim::Int = 2: Ambient space dimension
- manifold_dim::Int = 1: Manifold dimension of the structure surface
"""
@with_kw struct StructureDomain
    L::Float64 = 1.0
    x₀::Vector{Float64} = [0.0, 1.0]
        ambient_dim::Int = 2
        manifold_dim::Int = 1
    domain_symbol::Symbol = :Γ_s
end

"""
  DampingZone

Defines an axis-aligned damping-zone descriptor on a Cartesian free surface.

By default this represents a 1D segment embedded in 2D (`ambient_dim=2`,
`manifold_dim=1`), but the descriptor is dimension-aware and can be
parameterized for higher-dimensional embeddings.

Variables:
- L::Float64 = 0.5: Size of the damping zone along each manifold axis [m]
- σ::Float64 = 1.0: Damping coefficient
- x₀::Vector{Float64} = [0.0, 1.0]: Lower corner/anchor point [m]
- ambient_dim::Int = 2: Ambient space dimension
- manifold_dim::Int = 1: Manifold dimension of the damping-zone surface
"""
@with_kw struct DampingZone
    L::Float64 = 0.5
    σ::Float64 = 1.0
    x₀::Vector{Float64} = [0.0, 1.0]
        ambient_dim::Int = 2
        manifold_dim::Int = 1
    domain_symbol::Symbol = :Γ_d
end

"""
    JointDomain

Declares a rotational-spring joint located at a specific point on a 1D beam
boundary embedded in the 2D tank mesh.

`JointDomain` is purely a geometry-level descriptor.  After you add it to
`TankDomain2D.joint_domains`, `build_triangulations` automatically builds the
corresponding Gridap skeleton sub-triangulation (one interior facet of
`Skeleton(Γη)` whose centroid is closest to `location`).  Then
`get_integration_domains` populates `IntegrationDomains` with the resulting
measure and outward normal, ready to be consumed by
[`EulerBernoulliBeam`](@ref) via the matching `JointRotationalSpring`.

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

tank = TankDomain2D(L=4.0, H=1.0, nx=40, ny=4,
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

struct TankDomain{D, C, SZ, DZ, JZ} <: AbstractDomain
    cartesian::C
    structure_domains::SZ
    damping_zones::DZ
    joint_domains::JZ
end

const TankDomain2D = TankDomain{2}
const TankDomain3D = TankDomain{3}

function _validate_tank_domain_inputs(
    ::Val{2},
    structure_domains,
    damping_zones,
    joint_domains,
)
    nothing
end

function _validate_tank_domain_inputs(
    ::Val{3},
    structure_domains,
    damping_zones,
    joint_domains,
)
    isempty(structure_domains) ||
        error("TankDomain3D does not yet support structure_domains.")
    isempty(damping_zones) ||
        error("TankDomain3D does not yet support damping_zones.")
    isempty(joint_domains) ||
        error("TankDomain3D does not yet support joint_domains.")
    nothing
end

function TankDomain(
    cartesian;
    structure_domains = StructureDomain[],
    damping_zones = DampingZone[],
    joint_domains = JointDomain[],
)
    D = ambient_dimension(cartesian)
    _validate_tank_domain_inputs(
        Val(D),
        structure_domains,
        damping_zones,
        joint_domains,
    )
    TankDomain{
        D,
        typeof(cartesian),
        typeof(structure_domains),
        typeof(damping_zones),
        typeof(joint_domains),
    }(
        cartesian,
        structure_domains,
        damping_zones,
        joint_domains,
    )
end

function TankDomain(;
    L = 4.0,
    H = 1.0,
    nx = 16,
    ny = 2,
    W = nothing,
    nz = nothing,
    map = x -> x,
    structure_domains = StructureDomain[],
    damping_zones = DampingZone[],
    joint_domains = JointDomain[],
)
    cartesian = CartesianDomain(
        L = L,
        H = H,
        nx = nx,
        ny = ny,
        W = W,
        nz = nz,
        map = map,
    )
    TankDomain(
        cartesian;
        structure_domains = structure_domains,
        damping_zones = damping_zones,
        joint_domains = joint_domains,
    )
end

function TankDomain2D(;
    L = 4.0,
    H = 1.0,
    nx = 16,
    ny = 2,
    map = x -> x,
    structure_domains = StructureDomain[],
    damping_zones = DampingZone[],
    joint_domains = JointDomain[],
)
    TankDomain(
        L = L,
        H = H,
        nx = nx,
        ny = ny,
        map = map,
        structure_domains = structure_domains,
        damping_zones = damping_zones,
        joint_domains = joint_domains,
    )
end

function TankDomain3D(;
    L = 4.0,
    W = 2.0,
    H = 1.0,
    nx = 8,
    ny = 4,
    nz = 2,
    map = x -> x,
    structure_domains = StructureDomain[],
    damping_zones = DampingZone[],
    joint_domains = JointDomain[],
)
    TankDomain(
        L = L,
        W = W,
        H = H,
        nx = nx,
        ny = ny,
        nz = nz,
        map = map,
        structure_domains = structure_domains,
        damping_zones = damping_zones,
        joint_domains = joint_domains,
    )
end

function _cartesian_domain(domain::TankDomain)
    getfield(domain, :cartesian)
end

function _tank_legacy_dimensions(domain::TankDomain{2})
    cartesian = _cartesian_domain(domain)
    (
        L = cartesian.maxs[1] - cartesian.mins[1],
        H = cartesian.maxs[2] - cartesian.mins[2],
        nx = cartesian.parts[1],
        ny = cartesian.parts[2],
        map = cartesian.map,
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
    )
end

function Base.getproperty(domain::TankDomain{2}, name::Symbol)
    if name in (:L, :H, :nx, :ny, :map)
        return getproperty(_tank_legacy_dimensions(domain), name)
    end
    getfield(domain, name)
end

function Base.getproperty(domain::TankDomain{3}, name::Symbol)
    if name in (:L, :W, :H, :nx, :ny, :nz, :map)
        return getproperty(_tank_legacy_dimensions(domain), name)
    end
    getfield(domain, name)
end

function Base.propertynames(::TankDomain{2}, private::Bool = false)
    public = (
        :cartesian,
        :structure_domains,
        :damping_zones,
        :joint_domains,
        :L,
        :H,
        :nx,
        :ny,
        :map,
    )
    public
end

function Base.propertynames(::TankDomain{3}, private::Bool = false)
    public = (
        :cartesian,
        :structure_domains,
        :damping_zones,
        :joint_domains,
        :L,
        :W,
        :H,
        :nx,
        :ny,
        :nz,
        :map,
    )
    public
end

"""
  build_model(domain::TankDomain{2})

Build a `CartesianDiscreteModel` from the specifications in `domain`.
The model is constructed on the rectangular domain defined by `L` and `H`,
discretized into `nx` by `ny` elements.
"""
function build_model(domain::TankDomain{2})
    build_model(_cartesian_domain(domain))
end


# ─────────────────────────────────────────────────────────────
# Surface masking
# ─────────────────────────────────────────────────────────────

"""
    _centroid(xs) -> VectorValue

Compute the centroid of a cell given its node coordinates.
"""
function _centroid(xs)
    n = length(xs)
    (1 / n) * sum(xs)
end

_always_false_mask() = xs -> false

function _centroid_axis_eq_mask(axis::Int, value; tol = 1.0e-10)
    xs -> abs(_centroid(xs)[axis] - value) < tol
end

function _centroid_axis_either_eq_mask(axis::Int, value_a, value_b; tol = 1.0e-10)
    return function (xs)
        c = _centroid(xs)
        abs(c[axis] - value_a) < tol || abs(c[axis] - value_b) < tol
    end
end

function _surface_zone_mask(x0, L, ambient_dim::Int, manifold_dim::Int)
    return function (xs)
        c = _centroid(xs)
        in_plane = all(x0[i] <= c[i] <= x0[i] + L for i in 1:manifold_dim)
        on_surface = all(c[i] ≈ x0[i] for i in (manifold_dim + 1):ambient_dim)
        in_plane && on_surface
    end
end

function _surface_zone_mask(zone)
    _check_surface_zone(zone)
    _surface_zone_mask(zone.x₀, zone.L, zone.ambient_dim, zone.manifold_dim)
end

"""
    surface_mask(zone::StructureDomain) -> Function

Return a closure `(xs) -> Bool` that tests whether the centroid of
a surface cell lies within the structure domain.

For each manifold axis `i = 1:manifold_dim`, the centroid coordinate must
belong to `[x₀[i], x₀[i]+L]`. Remaining ambient coordinates are constrained to
`x₀[i]` using approximate equality.
"""
function surface_mask(zone::StructureDomain)
    _surface_zone_mask(zone)
end

"""
    surface_mask(zone::DampingZone) -> Function

Return a closure `(xs) -> Bool` that tests whether the centroid of
a surface cell lies within the damping zone.

For each manifold axis `i = 1:manifold_dim`, the centroid coordinate must
belong to `[x₀[i], x₀[i]+L]`. Remaining ambient coordinates are constrained to
`x₀[i]` using approximate equality.
"""
function surface_mask(zone::DampingZone)
    _surface_zone_mask(zone)
end

function _check_surface_zone(zone)
    length(zone.x₀) == zone.ambient_dim ||
        error("length(x₀)=$(length(zone.x₀)) must equal ambient_dim=$(zone.ambient_dim).")
    zone.manifold_dim >= 1 || error("manifold_dim must be >= 1.")
    zone.manifold_dim < zone.ambient_dim ||
        error("manifold_dim must be < ambient_dim for a boundary sub-domain.")
    nothing
end

"""
    surface_masks(domain::TankDomain{2})

Build mask closures for every structure domain and damping zone in `domain`.

Returns `(structure_masks, damping_masks)` where each is a `Vector{Function}`,
ordered to match `domain.structure_domains` and `domain.damping_zones`.
"""
function surface_masks(domain::TankDomain{2})
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


# ─────────────────────────────────────────────────────────────
# Triangulations
# ─────────────────────────────────────────────────────────────

function _unique_domain_symbols(zones, kind)
    syms = [z.domain_symbol for z in zones]
    seen = Set{Symbol}()
    for sym in syms
        if sym in seen
            error("Duplicate domain_symbol :$sym found in $kind. Each $(kind == "structure_domains" ? "structure" : "damping zone") must have a unique domain_symbol.")
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
            error("Duplicate joint domain_symbol :$(joint.domain_symbol) found in joint_domains. Each joint must have a unique domain_symbol.")
        end
        if joint.normal_symbol in seen_joint_normal_syms
            error("Duplicate joint normal_symbol :$(joint.normal_symbol) found in joint_domains. Each joint must have a unique normal_symbol.")
        end
        push!(seen_joint_domain_syms, joint.domain_symbol)
        push!(seen_joint_normal_syms, joint.normal_symbol)

        bits = lazy_map(joint_mask(joint), xΛη)
        joint_idxs = findall(bits)
        if isempty(joint_idxs)
            error("Joint at location $(joint.location) did not match any structure skeleton cell. Check location and tolerance.")
        end

        overlap_bits = joint_bits_all .& bits
        if any(overlap_bits)
            overlap_idxs = findall(overlap_bits)
            error("Joint domain_symbol :$(joint.domain_symbol) overlaps previously assigned joint skeleton facet(s) at index/indices $(overlap_idxs). Each skeleton facet may belong to at most one joint.")
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
            error("domain_symbol :$sym for $label triangulation would overwrite an existing triangulation key in TankTriangulations. All keys must be unique.")
        end
        trian_dict[sym] = tri
    end
    trian_dict
end

function _tank_triangulation_dict(;
    Ω,
    Γ,
    Γbot,
    Γin,
    Γout,
    Γ_structures,
    Γ_dampings,
    Γfs,
    Γκ,
    Γη,
    Λη,
    Λ_joints,
    joint_domains,
)
    Dict{Symbol, Any}(
        :Ω => Ω,
        :Γ => Γ,
        :Γbot => Γbot,
        :Γin => Γin,
        :Γout => Γout,
        :Γ_structures => Γ_structures,
        :Γ_dampings => Γ_dampings,
        :Γfs => Γfs,
        :Γκ => Γκ,
        :Γη => Γη,
        :Λη => Λη,
        :Λ_joints => Λ_joints,
        :joint_domains => joint_domains,
    )
end

"""
    build_triangulations(domain::TankDomain{2}, model) -> TankTriangulations

Label the model, extract base triangulations (`Ω`, `Γ`, `Γin`, `Γot`),
then partition the top surface `Γ` using the masks derived from the
structure domains and damping zones in `domain`.

Returns a [`TankTriangulations`](@ref) holding every sub-triangulation
needed for simulation assembly.

Labelling convention (Cartesian 2D, entity ids):
- `"surface"` → entities 3, 4, 6 (top side + top corners)
- `"bottom"`  → entities 1, 2, 5 (bottom side + bottom corners)
- `"inlet"`   → entity 7 (left wall)
- `"outlet"`  → entity 8 (right wall)
- `"water"`   → entity 9 (interior)
"""
function build_triangulations(domain::TankDomain{2}, model)
    # — Label model faces ————————————————————————————————
    _label_tank_model!(model)

    # — Base triangulations ——————————————————————————————
    Ω    = Interior(model)
    Γ    = Boundary(model, tags="surface")
    Γbot = Boundary(model, tags="bottom")
    Γin  = Boundary(model, tags="inlet")
    Γout = Boundary(model, tags="outlet")

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

    if !isempty(domain.joint_domains) && isempty(domain.structure_domains)
        error("Joint domains require at least one structure domain in TankDomain{2}.")
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

    # Compose dictionary for TankTriangulations
    trian_dict = _tank_triangulation_dict(
        Ω = Ω,
        Γ = Γ,
        Γbot = Γbot,
        Γin = Γin,
        Γout = Γout,
        Γ_structures = surface_partition.structures.triangulations,
        Γ_dampings = surface_partition.dampings.triangulations,
        Γfs = Γfs,
        Γκ = Γκ,
        Γη = Γη,
        Λη = joint_partition.Λη,
        Λ_joints = joint_partition.Λ_joints,
        joint_domains = domain.joint_domains,
    )
    _add_triangulations_by_symbol!(
        trian_dict,
        surface_partition.structures.symbols,
        surface_partition.structures.triangulations,
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


# ─────────────────────────────────────────────────────────────
# IntegrationDomains bridge
# ─────────────────────────────────────────────────────────────

"""
    get_integration_domains(tri::TankTriangulations; degree::Int=4) -> IntegrationDomains

Build an `IntegrationDomains` container of integration measures and normals
from `tri`.

# Populated keys

| Key | Source | Used by |
|-----|--------|---------|
| `:dΩ` | `tri.Ω` | PotentialFlow, Resonator |
| `:dΓfs` | `tri.Γfs` | Free surface without structures or damping zones |
| `:dΓ_s` | `tri.Γη` (all structures) | Membrane, Beam, PF↔struct coupling |
| `:dΓ_in` | `tri.Γin` | Inlet boundary terms |
| `:dΓ_out` | `tri.Γout` | Outlet boundary terms |
| `:dΓ_bot` | `tri.Γbot` | Bottom boundary terms |
| `:dΓ_d_i` | `tri.Γ_dampings[i]` | Damping zone terms (`:dΓ_d_1`, `:dΓ_d_2`, …) |

# Automatic beam skeleton keys

If at least one structure domain is present, this function also adds the beam
DG skeleton entries `:dΛη`, `:n_Λ_η`, and `:h_η`.

# Example

```julia
dom = G.get_integration_domains(tank_trians; degree=4)
```
"""
function get_integration_domains(tri::TankTriangulations; degree::Union{Int, Dict{Symbol, Int}}=4)

    d = Dict{Symbol, Any}()

    # Helper to get degree for a key
    get_deg(key) = isa(degree, Dict) ? get(degree, key, 4) : degree

    # Fluid interior
    d[:dΩ]     = Measure(tri[:Ω], get_deg(:dΩ))

    # Free surface without structures or damping zones
    d[:dΓfs] = Measure(tri[:Γfs], get_deg(:dΓfs))
    d[:nΓfs] = get_normal_vector(tri[:Γfs])

    # Free surface including damping zones, but excluding structures
    d[:dΓκ]  = Measure(tri[:Γκ], get_deg(:dΓκ))

    # All-structure surface
    d[:dΓη]   = Measure(tri[:Γη], get_deg(:dΓη))

    # Per-structure measures
    for (i, Γs) in enumerate(tri[:Γ_structures])
        key = Symbol("dΓη_$i")
        d[key] = Measure(Γs, get_deg(key))
    end

    # Walls
    d[:dΓin]  = Measure(tri[:Γin], get_deg(:dΓin))
    d[:dΓout] = Measure(tri[:Γout], get_deg(:dΓout))
    d[:dΓbot] = Measure(tri[:Γbot], get_deg(:dΓbot))

    # Per-damping-zone measures (:dΓd_1, :dΓd_2, …)
    for (i, Γd) in enumerate(tri[:Γ_dampings])
        key = Symbol("dΓd_$i")
        d[key] = Measure(Γd, get_deg(key))
        d[Symbol("nΓd_$i")] = get_normal_vector(Γd)
    end

    # Per-joint skeleton measures and normals
    if haskey(tri, :joint_domains) && haskey(tri, :Λ_joints)
        for (joint, Λj) in zip(tri[:joint_domains], tri[:Λ_joints])
            d[joint.domain_symbol] = Measure(Λj, get_deg(joint.domain_symbol))
            d[joint.normal_symbol] = get_normal_vector(Λj)
        end
    end

    # Beam DG skeleton keys (:dΛη, :n_Λη, :h_η) — required by EulerBernoulliBeam
    if !isempty(tri[:Γ_structures])
        Λη = haskey(tri, :Λη) && tri[:Λη] !== nothing ? tri[:Λη] : Skeleton(tri[:Γη])
        d[:dΛη]   = Measure(Λη, get_deg(:dΛη))
        d[:n_Λ_η] = get_normal_vector(Λη)
        # mesh size: minimum cell length on the structure surface
        cell_measures = get_cell_measure(tri[:Γη])
        d[:h_η] = minimum(cell_measures)
    end

    return IntegrationDomains(d)
end


# ─────────────────────────────────────────────────────────────
# AbstractDomain interface for TankDomain{2}
# ─────────────────────────────────────────────────────────────

"""
    _tank_label_map() -> Dict{String, String}

Internal mapping from [`STANDARD_TAGS`](@ref) to the Gridap face-label
strings used by [`build_model`](@ref) for Cartesian meshes.

`"fluid"` maps to the interior label `"water"`.  `"structure"` is handled
separately via coordinate masks in [`get_boundary`](@ref).
"""
function _tank_label_map()
  Dict{String, String}(
    "fluid"        => "water",
    "free_surface" => "surface",
    "seabed"       => "bottom",
    "inlet"        => "inlet",
    "outlet"       => "outlet",
    "structure"    => "structure",  # virtual; resolved via masks
  )
end

"""
    _label_tank_model!(model) -> nothing

Apply the standard Cartesian face labels to `model` in-place.

Adds `"surface"`, `"bottom"`, `"inlet"`, `"outlet"`, and `"water"` tags
using the fixed entity-id convention for 2D Cartesian meshes.
"""
function _label_tank_model!(model)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "surface", [3, 4, 6])
  add_tag_from_tags!(labels, "bottom",  [1, 2, 5])
  add_tag_from_tags!(labels, "inlet",   [7])
  add_tag_from_tags!(labels, "outlet",  [8])
  add_tag_from_tags!(labels, "water",   [9])
  nothing
end

"""
    ambient_dimension(d::TankDomain{2}) -> Int

Return 2: all `TankDomain{2}` problems live in a 2-D ambient space (x, z).
"""
ambient_dimension(::TankDomain{2}) = 2

"""
    manifold_dimension(d::TankDomain{2}) -> Int

Return 2: the fluid volume is a 2-D manifold.
"""
manifold_dimension(::TankDomain{2}) = 2

"""
    boundary_tags(d::TankDomain{2}) -> Dict{String, String}

Return a dictionary mapping each [`STANDARD_TAGS`](@ref) name to its
corresponding Gridap face-label string for the Cartesian mesh.

The `"structure"` key maps to `"structure"` (a virtual label resolved via
coordinate masks in [`get_boundary`](@ref) and [`build_triangulations`](@ref)).
"""
function boundary_tags(d::TankDomain{2})
  _tank_label_map()
end

"""
    triangulation(d::TankDomain{2}) -> Triangulation

Build the Cartesian discrete model from `d` and return the bulk-fluid
interior triangulation `Ω`.

Note: this method constructs a fresh `CartesianDiscreteModel` every call.
For performance-critical code use [`build_model`](@ref) /
[`build_triangulations`](@ref) directly.
"""
function triangulation(d::TankDomain{2})
  model = build_model(d)
  _label_tank_model!(model)
  Interior(model)
end

"""
    get_boundary(d::TankDomain{2}, name::String) -> BoundaryTriangulation

Return the Gridap boundary triangulation for the standard region `name`.

Supported `name` values: all six [`STANDARD_TAGS`](@ref).

For `"fluid"` the interior triangulation is returned.  For `"structure"` the
union of all structure-domain cells on the top surface is returned (empty
if no structure domains are defined).

Note: builds a fresh model on each call; use [`build_triangulations`](@ref)
when multiple boundaries are needed.
"""
function get_boundary(d::TankDomain{2}, name::String)
  valid = STANDARD_TAGS
  if !(name in valid)
    error(
      "Unknown boundary tag \"$name\" for TankDomain{2}. " *
      "Valid tags are: " * join(valid, ", ") * ".",
    )
  end
  model = build_model(d)
  _label_tank_model!(model)
  lmap  = _tank_label_map()

  if name == "fluid"
    return Interior(model)
  end

  if name == "structure"
    # Build the full top-surface triangulation, then sub-select structure cells
    Γ = Boundary(model, tags=lmap["free_surface"])
    if isempty(d.structure_domains)
      return Triangulation(Γ, Int[])
    end
    xΓ = get_cell_coordinates(Γ)
    structures = _partition_surface_zones(Γ, xΓ, d.structure_domains, "structure_domains")
    return Triangulation(Γ, findall(_or_bits(structures.bits, length(xΓ))))
  end

  Boundary(model, tags=lmap[name])
end

# ─────────────────────────────────────────────────────────────
# TankDomain3D — structured 3D Cartesian rectangular tank
# ─────────────────────────────────────────────────────────────

"""
    build_model(domain::TankDomain{3})

Build a `CartesianDiscreteModel` from the specifications in `domain`.
The model spans `[0,L]×[0,W]×[0,H]` with `nx × ny × nz` elements.
"""
function build_model(domain::TankDomain{3})
    build_model(_cartesian_domain(domain))
end

"""
        build_triangulations(domain::TankDomain{3}, model) -> TankTriangulations

Label the 3D model, extract sub-triangulations for each standard boundary
face via coordinate masks, and return a `TankTriangulations` container.

Boundary faces are identified by cell centroid coordinates:
- `free_surface` : z = H
- `seabed`       : z = 0
- `inlet`        : x = 0
- `outlet`       : x = L
- `lateral_walls`: y = 0 or y = W
- `structure`    : empty (no structure in basic 3D tank)
"""
function build_triangulations(domain::TankDomain{3}, model)
    d = _cartesian_domain(domain)
    _cartesian_boundary_triangulations(Val(3), model, d.mins, d.maxs)
end

# ─────────────────────────────────────────────────────────────
# AbstractDomain interface for TankDomain{3}
# ─────────────────────────────────────────────────────────────

"""
        ambient_dimension(d::TankDomain{3}) -> Int

Return 3: all `TankDomain{3}` problems live in a 3-D ambient space (x, y, z).
"""
ambient_dimension(::TankDomain{3}) = 3

"""
        manifold_dimension(d::TankDomain{3}) -> Int

Return 3: the fluid volume is a 3-D manifold.
"""
manifold_dimension(::TankDomain{3}) = 3

"""
        boundary_tags(d::TankDomain{3}) -> Dict{String, String}

Return a dictionary mapping each [`STANDARD_TAGS`](@ref) name to a
descriptive string for the 3D Cartesian tank.  The `"structure"` key maps
to `"structure"` (always empty for `TankDomain{3}`); the extra key
`"lateral_walls"` maps to `"lateral_walls"`.
"""
function boundary_tags(d::TankDomain{3})
    boundary_tags(_cartesian_domain(d))
end

"""
        triangulation(d::TankDomain{3}) -> Triangulation

Build the Cartesian discrete model from `d` and return the bulk-fluid
interior triangulation `Ω`.

Note: this method constructs a fresh `CartesianDiscreteModel` every call.
For performance-critical code use [`build_model`](@ref) /
[`build_triangulations`](@ref) directly.
"""
function triangulation(d::TankDomain{3})
    triangulation(_cartesian_domain(d))
end

"""
        get_boundary(d::TankDomain{3}, name::String) -> Triangulation

Return the sub-triangulation for the standard region `name`.

Supported names: all six [`STANDARD_TAGS`](@ref) plus `"lateral_walls"`.

For `"fluid"` the interior triangulation is returned.
For `"structure"` an empty sub-triangulation is returned.

Note: builds a fresh model on each call; use [`build_triangulations`](@ref)
when multiple boundaries are needed.
"""
function get_boundary(d::TankDomain{3}, name::String)
    get_boundary(_cartesian_domain(d), name)
end
