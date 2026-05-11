# ─────────────────────────────────────────────────────────────────────────────
# Triangulations — TankTriangulations container and low-level assembly helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    TankTriangulations

Dict-based container for all triangulations produced by
[`build_triangulations`](@ref).

Access entries with `trians[:key]` or `haskey(trians, :key)`.

## Standard keys

| Key              | Type             | Description                                          |
|------------------|------------------|------------------------------------------------------|
| `:Ω`             | `Triangulation`  | Interior (full fluid volume)                         |
| `:Γ`             | `Triangulation`  | Raw top-surface `Boundary` (all surface cells)       |
| `:Γfs`           | `Triangulation`  | Free surface (no structure, no damping zone)         |
| `:Γκ`            | `Triangulation`  | Free surface ∪ damping zones (excludes structure)    |
| `:Γη`            | `Triangulation`  | All-structure surface (union of all structure zones) |
| `:Γbot`          | `Triangulation`  | Bottom (seabed) boundary                             |
| `:Γin`           | `Triangulation`  | Inlet (left wall / wave-generation) boundary         |
| `:Γout`          | `Triangulation`  | Outlet (right wall / absorbing) boundary             |
| `:Γ_structures`  | `Vector`         | One entry per `StructureDomain`, same order          |
| `:Γ_dampings`    | `Vector`         | One entry per `DampingZone` / damping group          |
| `:Λη`            | `Triangulation`  | Beam skeleton on `Γη` (excluding joint facets)       |
| `:Λ_joints`      | `Vector`         | One skeleton sub-triangulation per `JointDomain`     |
| `:joint_domains` | `Vector`         | `JointDomain` objects in the same order as `:Λ_joints` |

Named sub-triangulations are also stored directly under their
`domain_symbol` keys (e.g. `:Γ_s_a` for a structure, `:Γ_d_1` for a damping
zone) so they can be accessed by physics entities without index arithmetic.
"""
struct TankTriangulations
  data::Dict{Symbol, Any}
end

TankTriangulations(; kwargs...) =
  TankTriangulations(Dict{Symbol, Any}(k => v for (k, v) in pairs(kwargs)))

Base.getindex(t::TankTriangulations, k::Symbol) = t.data[k]
Base.haskey(t::TankTriangulations, k::Symbol) = haskey(t.data, k)
Base.get(t::TankTriangulations, k::Symbol, default) = get(t.data, k, default)
Base.setindex!(t::TankTriangulations, val, k::Symbol) = (t.data[k] = val)
Base.keys(t::TankTriangulations) = keys(t.data)

# ─────────────────────────────────────────────────────────────────────────────
# Internal dict constructor used by both CartesianDomain and TankDomain paths
# ─────────────────────────────────────────────────────────────────────────────

"""
    _tank_triangulation_dict(; Ω, Γ, Γbot, Γin, Γout,
                               Γ_structures, Γ_dampings,
                               Γfs, Γκ, Γη, Λη, Λ_joints, joint_domains)
                             -> Dict{Symbol,Any}

Assemble the canonical `Dict{Symbol,Any}` used to construct a
[`TankTriangulations`](@ref).  All keyword arguments are required.

This helper is called by both the Cartesian and TankDomain assembly paths so
the key set is always consistent.
"""
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
    :Ω            => Ω,
    :Γ            => Γ,
    :Γbot         => Γbot,
    :Γin          => Γin,
    :Γout         => Γout,
    :Γ_structures => Γ_structures,
    :Γ_dampings   => Γ_dampings,
    :Γfs          => Γfs,
    :Γκ           => Γκ,
    :Γη           => Γη,
    :Λη           => Λη,
    :Λ_joints     => Λ_joints,
    :joint_domains => joint_domains,
  )
end
