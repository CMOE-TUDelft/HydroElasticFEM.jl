"""
  StructureDomain1D

Defines a rectangular domain for 1D structural problems on a conforming 2D domain.
The domain is defined by its length `L` and an optional starting point `x₀` in 2D space.

Variables:
- L::Float64 = 1.0: Length of the 1D structure domain
- x₀::Vector{Float64} = [0.0, 1.0]: Starting point of the structure domain in 2D space
"""
@with_kw struct StructureDomain1D
    L::Float64 = 1.0
    x₀::Vector{Float64} = [0.0, 1.0]
    domain_symbol::Symbol = :Γ_s
end

"""
  DampingZone1D

Defines a 1D damping zone for absorbing outgoing waves, with length `L` and a damping coefficient `σ`.
The damping is applied in the region defined by `x₀` to `x₀ + L` along the x-axis.

Variables:
- L::Float64 = 0.5: Length of the damping zone
- σ::Float64 = 1.0: Damping coefficient
- x₀::Vector{Float64} = [0.0, 1.0]: Starting point of the damping zone in 2D space
"""
@with_kw struct DampingZone1D
    L::Float64 = 0.5
    σ::Float64 = 1.0
    x₀::Vector{Float64} = [0.0, 1.0]
    domain_symbol::Symbol = :Γ_d
end

"""
    JointDomain1D

Declares a rotational-spring joint located at a specific point on a 1D beam
boundary embedded in the 2D tank mesh.

`JointDomain1D` is purely a geometry-level descriptor.  After you add it to
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
s1 = StructureDomain1D(L=1.0, x₀=[1.5, 1.0])
j1 = JointDomain1D(location=[2.0, 1.0], domain_symbol=:dΛj_1, normal_symbol=:n_Λ_j_1)

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
@with_kw struct JointDomain1D
        location::Vector{Float64}
        domain_symbol::Symbol
        normal_symbol::Symbol
        tol::Float64 = 1.0e-6
end

"""
  TankDomain2D <: AbstractDomain

Defines a rectangular domain for 2D problems, with dimensions `L` and `H`
`[0,L]×[0,H]`, discretized into `nx` by `ny` elements.

The `map` function allows for coordinate transformations if needed.
`TankDomain2D` implements the [`AbstractDomain`](@ref) interface; use
[`triangulation`](@ref), [`get_boundary`](@ref), [`boundary_tags`](@ref),
[`ambient_dimension`](@ref), and [`manifold_dimension`](@ref) for generic code.

# Fields
- `L::Float64 = 4.0`  — length of the tank domain [m]
- `H::Float64 = 1.0`  — height of the tank domain [m]
- `nx::Int = 16`      — number of elements in the x-direction
- `ny::Int = 2`       — number of elements in the y-direction
- `map::Function`     — optional coordinate mapping (default: identity)
- `structure_domains::Vector{StructureDomain1D}` — structure sub-domains
- `damping_zones::Vector{DampingZone1D}` — damping sub-zones
- `joint_domains::Vector{JointDomain1D}` — structural joint descriptors
"""
@with_kw struct TankDomain2D <: AbstractDomain
    L::Float64 = 4.0
    H::Float64 = 1.0
    nx::Int = 16
    ny::Int = 2
    map::Function = x->x
    structure_domains::Vector{StructureDomain1D} = Vector{StructureDomain1D}()
    damping_zones::Vector{DampingZone1D} = Vector{DampingZone1D}()
    joint_domains::Vector{JointDomain1D} = Vector{JointDomain1D}()
end

"""
  build_model(domain::TankDomain2D)

Build a `CartesianDiscreteModel` from the specifications in `domain`.
The model is constructed on the rectangular domain defined by `L` and `H`, 
discretized into `nx` by `ny` elements.
"""
function build_model(domain::TankDomain2D)
    x_min, x_max = 0.0, domain.L
    y_min, y_max = 0.0, domain.H
    model = CartesianDiscreteModel((x_min, x_max, y_min, y_max), (domain.nx, domain.ny),map=domain.map)
    return model
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

"""
    surface_mask(zone::StructureDomain1D) -> Function

Return a closure `(xs) -> Bool` that tests whether the centroid of
a surface cell lies within the structure domain.

The zone spans `[x₀[1], x₀[1]+L]` at `y = x₀[2]`.
"""
function surface_mask(zone::StructureDomain1D)
    x_lo = zone.x₀[1]
    x_hi = zone.x₀[1] + zone.L
    y_ref = zone.x₀[2]
    return function (xs)
        c = _centroid(xs)
        (x_lo <= c[1] <= x_hi) && (c[2] ≈ y_ref)
    end
end

"""
    surface_mask(zone::DampingZone1D) -> Function

Return a closure `(xs) -> Bool` that tests whether the centroid of
a surface cell lies within the damping zone.

The zone spans `[x₀[1], x₀[1]+L]` at `y = x₀[2]`.
"""
function surface_mask(zone::DampingZone1D)
    x_lo = zone.x₀[1]
    x_hi = zone.x₀[1] + zone.L
    y_ref = zone.x₀[2]
    return function (xs)
        c = _centroid(xs)
        (x_lo <= c[1] <= x_hi) && (c[2] ≈ y_ref)
    end
end

"""
    surface_masks(domain::TankDomain2D)

Build mask closures for every structure domain and damping zone in `domain`.

Returns `(structure_masks, damping_masks)` where each is a `Vector{Function}`,
ordered to match `domain.structure_domains` and `domain.damping_zones`.
"""
function surface_masks(domain::TankDomain2D)
    smasks = [surface_mask(s) for s in domain.structure_domains]
    dmasks = [surface_mask(d) for d in domain.damping_zones]
    return smasks, dmasks
end

"""
    joint_mask(joint::JointDomain1D) -> Function

Return a closure `(xs) -> Bool` that selects skeleton cells whose centroid is
at `joint.location` within the joint tolerance.
"""
function joint_mask(joint::JointDomain1D)
    xj = joint.location[1]
    yj = joint.location[2]
    tol = joint.tol
    return function (xs)
        c = _centroid(xs)
        (abs(c[1] - xj) <= tol) && (abs(c[2] - yj) <= tol)
    end
end


# ─────────────────────────────────────────────────────────────
# Triangulations
# ─────────────────────────────────────────────────────────────

"""
    build_triangulations(domain::TankDomain2D, model) -> TankTriangulations

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
function build_triangulations(domain::TankDomain2D, model)
    # — Label model faces ————————————————————————————————
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "surface", [3, 4, 6])
    add_tag_from_tags!(labels, "bottom",  [1, 2, 5])
    add_tag_from_tags!(labels, "inlet",   [7])
    add_tag_from_tags!(labels, "outlet",  [8])
    add_tag_from_tags!(labels, "water",   [9])

    # — Base triangulations ——————————————————————————————
    Ω    = Interior(model)
    Γ    = Boundary(model, tags="surface")
    Γbot = Boundary(model, tags="bottom")
    Γin  = Boundary(model, tags="inlet")
    Γout = Boundary(model, tags="outlet")

    # — Build boolean masks on surface cell coordinates ——
    xΓ = get_cell_coordinates(Γ)
    smasks, dmasks = surface_masks(domain)

    s_bits = [lazy_map(m, xΓ) for m in smasks]  # per-structure
    d_bits = [lazy_map(m, xΓ) for m in dmasks]  # per-damping

    # Build structure triangulations and map to their domain_symbol
    Γ_structures = Vector{Any}(undef, length(s_bits))
    structure_syms = Symbol[]
    structure_trians = Dict{Symbol, Any}()
    seen_structure_syms = Set{Symbol}()
    for (i, b) in enumerate(s_bits)
        tri = Triangulation(Γ, findall(b))
        Γ_structures[i] = tri
        sym = domain.structure_domains[i].domain_symbol
        if sym in seen_structure_syms
            error("Duplicate domain_symbol :$sym found in structure_domains. Each structure must have a unique domain_symbol.")
        end
        push!(seen_structure_syms, sym)
        structure_trians[sym] = tri
        push!(structure_syms, sym)
    end

    # Damping triangulations and map to their domain_symbol
    Γ_dampings = Vector{Any}(undef, length(d_bits))
    damping_syms = Symbol[]
    damping_trians = Dict{Symbol, Any}()
    seen_damping_syms = Set{Symbol}()
    for (i, b) in enumerate(d_bits)
        tri = Triangulation(Γ, findall(b))
        Γ_dampings[i] = tri
        sym = domain.damping_zones[i].domain_symbol
        if sym in seen_damping_syms
            error("Duplicate domain_symbol :$sym found in damping_zones. Each damping zone must have a unique domain_symbol.")
        end
        push!(seen_damping_syms, sym)
        damping_trians[sym] = tri
        push!(damping_syms, sym)
    end

    # — Composed masks ——————————————————————————————————
    n = length(xΓ)
    any_structure = _or_bits(s_bits, n)
    any_damping   = _or_bits(d_bits, n)
    any_zone      = any_structure .| any_damping

    # Free surface = surface cells in no zone at all
    Γfs = Triangulation(Γ, findall(!, any_zone))
    # κ surface = everything except structures (free surface + damping)
    Γκ  = Triangulation(Γ, findall(!, any_structure))
    # η surface = all structures
    Γη  = Triangulation(Γ, findall(any_structure))

    if !isempty(domain.joint_domains) && isempty(domain.structure_domains)
        error("Joint domains require at least one structure domain in TankDomain2D.")
    end

    # Joint skeleton triangulations from the structure skeleton
    Λ_joints = Any[]
    Λη_no_joints = nothing
    joint_domain_syms = Symbol[]
    if !isempty(domain.joint_domains)
        Λη = Skeleton(Γη)
        xΛη = get_cell_coordinates(Λη)
        joint_bits_all = falses(length(xΛη))
        non_joint_bits = trues(length(xΛη))
        seen_joint_domain_syms = Set{Symbol}()
        seen_joint_normal_syms = Set{Symbol}()
        for joint in domain.joint_domains
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
            Λj = Triangulation(Λη, joint_idxs)
            push!(Λ_joints, Λj)
            push!(joint_domain_syms, joint.domain_symbol)
        end

        # Beam DG skeleton excludes joint facets via the complementary mask
        # to avoid double counting with explicit rotational-spring terms.
        Λη_no_joints = Triangulation(Λη, findall(non_joint_bits))

        # Partition sanity check: non-joint + joint facets must cover Λη.
        if (num_cells(Λη_no_joints) + count(joint_bits_all)) != length(xΛη)
            error("Joint/non-joint skeleton partition is inconsistent for Γη.")
        end
    elseif !isempty(domain.structure_domains)
        # No joints: use full structure skeleton.
        Λη_no_joints = Skeleton(Γη)
    end

    # Compose dictionary for TankTriangulations
    trian_dict = Dict(
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
        :Λη => Λη_no_joints,
        :Λ_joints => Λ_joints,
        :joint_domains => domain.joint_domains,
    )
    # Add each structure triangulation under its domain_symbol, error if duplicate
    for sym in structure_syms
        if haskey(trian_dict, sym)
            error("domain_symbol :$sym for structure triangulation would overwrite an existing triangulation key in TankTriangulations. All keys must be unique.")
        end
        trian_dict[sym] = structure_trians[sym]
    end
    # Add each damping triangulation under its domain_symbol, error if duplicate
    for sym in damping_syms
        if haskey(trian_dict, sym)
            error("domain_symbol :$sym for damping triangulation would overwrite an existing triangulation key in TankTriangulations. All keys must be unique.")
        end
        trian_dict[sym] = damping_trians[sym]
    end
    # Add each joint skeleton triangulation under its domain_symbol, error if duplicate
    for (i, sym) in enumerate(joint_domain_syms)
        if haskey(trian_dict, sym)
            error("domain_symbol :$sym for joint triangulation would overwrite an existing triangulation key in TankTriangulations. All keys must be unique.")
        end
        trian_dict[sym] = Λ_joints[i]
    end
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
| `:dΓ_s` | `tri.Γη` (all structures) | Membrane2D, Beam, PF↔struct coupling |
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
# AbstractDomain interface for TankDomain2D
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
    ambient_dimension(d::TankDomain2D) -> Int

Return 2: all `TankDomain2D` problems live in a 2-D ambient space (x, z).
"""
ambient_dimension(::TankDomain2D) = 2

"""
    manifold_dimension(d::TankDomain2D) -> Int

Return 2: the fluid volume is a 2-D manifold.
"""
manifold_dimension(::TankDomain2D) = 2

"""
    boundary_tags(d::TankDomain2D) -> Dict{String, String}

Return a dictionary mapping each [`STANDARD_TAGS`](@ref) name to its
corresponding Gridap face-label string for the Cartesian mesh.

The `"structure"` key maps to `"structure"` (a virtual label resolved via
coordinate masks in [`get_boundary`](@ref) and [`build_triangulations`](@ref)).
"""
function boundary_tags(d::TankDomain2D)
  _tank_label_map()
end

"""
    triangulation(d::TankDomain2D) -> Triangulation

Build the Cartesian discrete model from `d` and return the bulk-fluid
interior triangulation `Ω`.

Note: this method constructs a fresh `CartesianDiscreteModel` every call.
For performance-critical code use [`build_model`](@ref) /
[`build_triangulations`](@ref) directly.
"""
function triangulation(d::TankDomain2D)
  model = build_model(d)
  _label_tank_model!(model)
  Interior(model)
end

"""
    get_boundary(d::TankDomain2D, name::String) -> BoundaryTriangulation

Return the Gridap boundary triangulation for the standard region `name`.

Supported `name` values: all six [`STANDARD_TAGS`](@ref).

For `"fluid"` the interior triangulation is returned.  For `"structure"` the
union of all structure-domain cells on the top surface is returned (empty
if no structure domains are defined).

Note: builds a fresh model on each call; use [`build_triangulations`](@ref)
when multiple boundaries are needed.
"""
function get_boundary(d::TankDomain2D, name::String)
  valid = STANDARD_TAGS
  if !(name in valid)
    error(
      "Unknown boundary tag \"$name\" for TankDomain2D. " *
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
    smasks, _ = surface_masks(d)
    n = length(xΓ)
    any_structure = _or_bits([lazy_map(m, xΓ) for m in smasks], n)
    return Triangulation(Γ, findall(any_structure))
  end

  Boundary(model, tags=lmap[name])
end
