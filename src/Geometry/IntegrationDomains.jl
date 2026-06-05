# ─────────────────────────────────────────────────────────────────────────────
# IntegrationDomains — measure container and assembly bridge
# ─────────────────────────────────────────────────────────────────────────────

"""
    IntegrationDomains

Dict-based container for Gridap measures, outward normals, `DiracDelta`
functionals, and any other domain-level data consumed by weak forms.

Each physics entity reads only the keys it needs via `dom[:key]`.
No schema is enforced; new keys can be added without changing this type.

## Standard key conventions

### Fluid interior
| Key    | Type      | Description                                   |
|--------|-----------|-----------------------------------------------|
| `:dΩ`  | `Measure` | Fluid interior quadrature measure             |

### Boundary measures and normals
| Key       | Type             | Description                                |
|-----------|------------------|--------------------------------------------|
| `:dΓfs`   | `Measure`        | Free surface (no structure / damping zone) |
| `:nΓfs`   | `GenericCellField` | Outward normal on `:Γfs`                 |
| `:dΓκ`    | `Measure`        | Free surface ∪ damping zones              |
| `:dΓη`    | `Measure`        | All-structure surface                      |
| `:dΓη_i`  | `Measure`        | Per-structure measure (`:dΓη_1`, `:dΓη_2`, …)|
| `:dΓin`   | `Measure`        | Inlet boundary                             |
| `:dΓout`  | `Measure`        | Outlet boundary                            |
| `:dΓbot`  | `Measure`        | Bottom (seabed) boundary                   |

### Damping zones
| Key        | Type             | Description                              |
|------------|------------------|------------------------------------------|
| `:dΓd_i`   | `Measure`        | `i`-th damping-zone measure (1-indexed)  |
| `:nΓd_i`   | `GenericCellField` | Outward normal on `i`-th damping zone  |

### Beam DG skeleton
| Key       | Type             | Description                                |
|-----------|------------------|--------------------------------------------|
| `:dΛη`    | `Measure`        | Interior-facet (skeleton) measure on `Γη`  |
| `:n_Λ_η`  | `GenericCellField` | Skeleton outward normal                  |
| `:h_η`    | `Float64`        | Minimum cell length on `Γη` (mesh size)    |

### Joint skeleton (one per `JointDomain`)
Stored under `joint.domain_symbol` / `joint.normal_symbol` as declared in
[`JointDomain`](@ref).

### Resonators
| Key    | Type     | Description                                |
|--------|----------|--------------------------------------------|
| `:δ_p` | `Vector` | `DiracDelta` functionals (one per resonator)|

## Usage

```julia
model  = G.build_model(tank)
trians = G.build_triangulations(tank, model)
dom    = G.get_integration_domains(trians; degree=4)

# physics assembly
a(u, v) = ∫(u * v) * dom[:dΩ]
```
"""
struct IntegrationDomains
  data::Dict{Symbol, Any}
end

# Keyword-arg constructor: build from named measure pairs (e.g., IntegrationDomains(dΩ=dΩ, dΓ=dΓ)).
IntegrationDomains(; kwargs...) =
  IntegrationDomains(Dict{Symbol, Any}(k => v for (k, v) in pairs(kwargs)))

# Dict-like interface delegating to the internal `data` dictionary.
Base.getindex(d::IntegrationDomains, k::Symbol) = d.data[k]
Base.haskey(d::IntegrationDomains, k::Symbol) = haskey(d.data, k)
Base.get(d::IntegrationDomains, k::Symbol, default) = get(d.data, k, default)
Base.setindex!(d::IntegrationDomains, val, k::Symbol) = (d.data[k] = val)
Base.keys(d::IntegrationDomains) = keys(d.data)

# Convert a triangulation key into the matching measure key used by weak forms.
# Examples: :Ω -> :dΩ, :Γ_s -> :dΓ_s, :Λη -> :dΛη.
_measure_key(domain_symbol::Symbol) = Symbol("d", domain_symbol)

# Convert a measure key back to its triangulation key.  This lets callers pass
# either :Γ_s or :dΓ_s in a degree dictionary and get the same quadrature order.
_domain_key(measure_symbol::Symbol) = Symbol(string(measure_symbol)[2:end])

# Recognize the measure-key families created by this module.  Avoid treating an
# arbitrary user domain beginning with lowercase `d` as a measure key.
function _is_measure_key(key::Symbol)
  s = string(key)
  startswith(s, "dΩ") || startswith(s, "dΓ") || startswith(s, "dΛ")
end

_degree_for(degree::Int, ::Symbol) = degree

# Resolve quadrature degree for both supported dictionary styles:
#   Dict(:Γ_s => 6)  and  Dict(:dΓ_s => 6).
# Missing entries keep the historical fallback degree of 4.
function _degree_for(degree::Dict{Symbol, Int}, key::Symbol)
  domain_key = _is_measure_key(key) ? _domain_key(key) : key
  measure_key = _is_measure_key(key) ? key : _measure_key(key)
  return get(degree, key, get(degree, measure_key, get(degree, domain_key, 4)))
end

# Extract triangulation-domain keys from the degree dictionary.  These are the
# domains requested by physics entities through `space_domain_symbol`.
function _degree_domain_keys(degree::Dict{Symbol, Int})
  keys = Symbol[]
  for key in Base.keys(degree)
    push!(keys, _is_measure_key(key) ? _domain_key(key) : key)
  end
  unique(keys)
end

_degree_domain_keys(::Int) = Symbol[]

# Only scalar triangulations can be wrapped in a Measure here.  Container and
# metadata entries such as :Γ_structures, :Γ_dampings, and :joint_domains are
# handled explicitly above or intentionally skipped.
_can_define_measure(trian) = !(trian === nothing || trian isa AbstractVector || trian isa Tuple)

# ─────────────────────────────────────────────────────────────────────────────
# Assembly bridge
# ─────────────────────────────────────────────────────────────────────────────

"""
    get_integration_domains(
        tri::TankTriangulations;
        degree::Union{Int, Dict{Symbol,Int}} = 4,
    ) -> IntegrationDomains

Build an [`IntegrationDomains`](@ref) from a [`TankTriangulations`](@ref).

Pass `degree` as a single `Int` to use the same quadrature order everywhere,
or as a `Dict{Symbol,Int}` to override it per key (missing keys fall back
to `4`).

## Populated keys

| Key(s)              | Source triangulation           | Consumer            |
|---------------------|-------------------------------|---------------------|
| `:dΩ`               | `:Ω`                          | PotentialFlow, Resonator |
| `:dΓfs`, `:nΓfs`    | `:Γfs`                        | FreeSurface, damping-free BCs |
| `:dΓκ`              | `:Γκ`                         | FreeSurface (open water + damping) |
| `:dΓη`, `:dΓη_i`    | `:Γη`, `:Γ_structures[i]`     | Structure physics   |
| `:dΓin`, `:dΓout`, `:dΓbot` | `:Γin`, `:Γout`, `:Γbot` | Wall/radiation BCs  |
| `:dΓd_i`, `:nΓd_i`  | `:Γ_dampings[i]`               | `DampingZoneBC`     |
| `:dΛη`, `:n_Λ_η`, `:h_η` | skeleton of `:Γη`         | `EulerBernoulliBeam` DG |
| `joint.domain_symbol`, `joint.normal_symbol` | `:Λ_joints` | `JointRotationalSpring` |

## Example

```julia
model  = G.build_model(tank)
trians = G.build_triangulations(tank, model)
dom    = G.get_integration_domains(trians)
dom    = G.get_integration_domains(trians; degree=Dict(:dΩ => 6, :dΓη => 8))
```
"""
function get_integration_domains(
  tri::TankTriangulations;
  degree::Union{Int, Dict{Symbol, Int}} = 4,
)
  d = Dict{Symbol, Any}()

  get_deg(key) = _degree_for(degree, key)

  # Fluid interior
  d[:dΩ] = Measure(tri[:Ω], get_deg(:dΩ))

  # Free surface (no structure, no damping)
  d[:dΓfs] = Measure(tri[:Γfs], get_deg(:dΓfs))
  d[:nΓfs] = get_normal_vector(tri[:Γfs])

  # Free surface ∪ damping zones (excludes structure)
  d[:dΓκ] = Measure(tri[:Γκ], get_deg(:dΓκ))

  # All-structure surface
  d[:dΓη] = Measure(tri[:Γη], get_deg(:dΓη))

  # Per-structure measures (:dΓη_1, :dΓη_2, …)
  for (i, Γs) in enumerate(tri[:Γ_structures])
    key = Symbol("dΓη_$i")
    d[key] = Measure(Γs, get_deg(key))
  end

  # Walls
  d[:dΓin]  = Measure(tri[:Γin],  get_deg(:dΓin))
  d[:dΓout] = Measure(tri[:Γout], get_deg(:dΓout))
  d[:dΓbot] = Measure(tri[:Γbot], get_deg(:dΓbot))

  # Per-damping-zone measures (:dΓd_1, :dΓd_2, …)
  for (i, Γd) in enumerate(tri[:Γ_dampings])
    d[Symbol("dΓd_$i")] = Measure(Γd, get_deg(Symbol("dΓd_$i")))
    d[Symbol("nΓd_$i")] = get_normal_vector(Γd)
  end

  # Per-joint skeleton measures and normals (stored under domain_symbol /
  # normal_symbol declared in each JointDomain)
  if haskey(tri, :joint_domains) && haskey(tri, :Λ_joints)
    for (joint, Λj) in zip(tri[:joint_domains], tri[:Λ_joints])
      d[joint.domain_symbol] = Measure(Λj, get_deg(joint.domain_symbol))
      d[joint.normal_symbol] = get_normal_vector(Λj)
    end
  end

  # Beam DG skeleton (:dΛη, :n_Λ_η, :h_η) — required by EulerBernoulliBeam
  if !isempty(tri[:Γ_structures])
    Λη = (haskey(tri, :Λη) && tri[:Λη] !== nothing) ? tri[:Λη] : Skeleton(tri[:Γη])
    d[:dΛη]   = Measure(Λη, get_deg(:dΛη))
    d[:n_Λ_η] = get_normal_vector(Λη)
    d[:h_η]   = minimum(get_cell_measure(tri[:Γη]))
  end

  # User-defined space domains declared by physics entities arrive through
  # `degree` as triangulation keys (e.g. `:Γ_s`).  Mirror each such scalar
  # triangulation into the matching measure key (`:dΓ_s`) unless a standard
  # measure above already populated it.
  for domain_symbol in _degree_domain_keys(degree)
    haskey(tri, domain_symbol) || continue
    measure_symbol = _measure_key(domain_symbol)
    haskey(d, measure_symbol) && continue
    trian = tri[domain_symbol]
    _can_define_measure(trian) || continue
    d[measure_symbol] = Measure(trian, get_deg(domain_symbol))
  end

  IntegrationDomains(d)
end
