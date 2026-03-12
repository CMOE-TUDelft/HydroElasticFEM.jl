# Geometry Module — Development Plan

## 1. Motivation

Every script (legacy and current) repeats the same boilerplate:
domain rectangle → `CartesianDiscreteModel` → label entities → mask surface
into sub-triangulations → create measures → build `WeakFormDomains`.
This plan captures a module that encapsulates those steps so that a user
(or `SimManager`) only declares *what* is in the tank and gets back a
fully-populated `WeakFormDomains` ready for `Physics`.

---

## 2. Scope — What the module owns

| Responsibility | Currently in… |
|---|---|
| Rectangular domain definition (x₀, LΩ, H₀) | every script |
| Mesh grading map (geometric progression in y) | every script (`f_y`) |
| `CartesianDiscreteModel` construction | every script |
| Face-labelling (`surface`, `bottom`, `inlet`, `outlet`, `water`) | every script |
| Surface masking (structure zone, damping zones, free surface) | every script |
| Sub-triangulation creation (`Γm`, `Γfs`, `Γd1`, `Γd2`, `Γκ`, `Γη`) | every script |
| Structural-boundary labels (`mem_bnd`, `beam_bnd`) | every script |
| Skeleton triangulation + normals for DG beams | beam scripts |
| Integration measure creation (`dΩ`, `dΓm`, …) | every script |
| DiracDelta construction from resonator positions | LRHS/resonator scripts |
| Packing everything into `WeakFormDomains` | every script |
| Optional VTK mesh dump | every script |

**Not in scope (stays in user scripts / other modules):**
wave input functions, physics parameters, FE-space construction
(`FESpaceAssembly`), weak-form assembly, time-stepping, solvers,
post-processing.

---

## 3. Proposed module structure

```
src/Geometry/
├── Geometry.jl            # module wrapper, exports
├── TankDomain2D.jl        # parameter struct + model builder
├── SurfaceLayout.jl       # masking / triangulation partitioning
├── MeshGrading.jl         # f_y geometric-progression map + future extensions
├── DampingZone.jl         # inlet / outlet zone descriptor
└── GeometryToPhysics.jl   # build WeakFormDomains from geometry output
```

Integration into the package:

```julia
# src/HydroElasticFEM.jl
include("Geometry/Geometry.jl")
using .Geometry
```

`Geometry` depends only on `Gridap` and `Physics.Domains`
(for `WeakFormDomains`). It does **not** depend on `Entities` or
`FESpaces` — it sits beside `Physics`, not inside it.

---

## 4. Key types

### 4.1 `TankDomain2D`  — top-level config

```julia
@with_kw struct TankDomain2D
    H0::Float64                        # Still-water depth
    Ld_in::Float64                     # Inlet damping zone length
    Ld_out::Float64                    # Outlet damping zone length
    L_free::Float64                    # Free-water length before structure
    structures::Vector{StructureZone}  # Ordered left→right
    L_after::Float64 = L_free          # Free-water length after last structure
    nx::Int                            # Horizontal element count
    ny::Int                            # Vertical element count
    mesh_ry::Float64 = 1.2            # Vertical grading ratio (1.0 = uniform)
    order::Int = 2                     # Polynomial order
end
```

Derived quantities computed on construction:
- `x₀ = -Ld_in`
- `LΩ` = sum of all zone widths
- `domain = (x₀, x₀+LΩ, -H0, 0.0)`
- element size `h = LΩ / nx`

### 4.2 `StructureZone` — one structural member on the surface

```julia
struct StructureZone
    label::Symbol            # e.g. :membrane, :beam
    L::Float64               # length along x
    L_gap::Float64           # gap before next structure (0.0 if last)
end
```

Multiple structures are placed sequentially starting after
`L_free` from the inlet damping edge.

### 4.3 `DampingZone` — spatial damping profile

```julia
@with_kw struct DampingZone
    μ₀::Float64 = 2.5                # Peak damping amplitude
    profile::Symbol = :sinusoidal     # :sinusoidal | :polynomial | :custom
end
```

The module provides ready-made damping-coefficient functions
`μ₁(x)` and `μ₂(x) = μ₁(x)·k` (wave-number scaled) for each zone.

### 4.4 `TankGeometry` — output of the build step

```julia
struct TankGeometry
    model::CartesianDiscreteModel
    labels::FaceLabeling

    # Base triangulations
    Ω::Triangulation       # Interior
    Γ::Triangulation       # Full surface
    Γin::Triangulation     # Inlet boundary
    Γot::Triangulation     # Outlet boundary (if needed)

    # Surface partition
    Γfs::Triangulation     # Free surface (no structure, no damping)
    Γd_in::Triangulation   # Inlet damping zone on surface
    Γd_out::Triangulation  # Outlet damping zone on surface
    structures::OrderedDict{Symbol, StructureTriangulations}

    # Composed (κ = surface − all structures, η = union of structures)
    Γκ::Triangulation
    Γη::Triangulation

    # Measures (degree = 2·order)
    measures::Dict{Symbol, Any}  # :dΩ, :dΓ_fs, etc.

    # Normals for skeleton / boundary terms
    normals::Dict{Symbol, Any}   # :n_Λ_s, :n_Λ_sb, etc.
end
```

Where each structure gets:

```julia
struct StructureTriangulations
    Γs::Triangulation      # Surface of this structure
    Λs::Triangulation      # Interior skeleton (DG beams)
    Λsb::Triangulation     # Boundary edges (fixed-BC Neumann)
    dΓs::Measure
    dΛs::Measure
    dΛsb::Measure
    n_Λs                   # Normal on skeleton
    n_Λsb                  # Normal on boundary
    h::Float64             # Element size on structure surface
end
```

---

## 5. Public API

### 5.1 Build geometry

```julia
tank  = TankDomain2D(H0=10, Ld_in=150, Ld_out=150,
          L_free=80,
          structures=[StructureZone(:membrane, 20.0, 0.0)],
          nx=1000, ny=20)
geom  = build_geometry(tank)           # → TankGeometry
```

### 5.2 Produce `WeakFormDomains` for Physics

```julia
dom = to_weakform_domains(geom)        # → WeakFormDomains
# Populates :dΩ, :dΓ_fs, :dΓ_s, :dΛ_s, :n_Λ_s, :h_s,
#   :dΛ_sb, :n_Λ_sb, :dΓ_in, :dΓ_ot
```

For resonator points:

```julia
dom = to_weakform_domains(geom, resonators)
# Adds :δ_p => [DiracDelta(…), …]
```

### 5.3 Damping coefficient functions

```julia
μ_in, μ_out = damping_functions(tank; k=k)
# μ_in(x::VectorValue) -> Float64   (inlet zone)
# μ_out(x::VectorValue) -> Float64  (outlet zone)
```

### 5.4 VTK mesh dump

```julia
write_geometry_vtk(geom, "output/mesh")
# Writes model, Ω, Γ, Γm, Γfs, Γd1, Γd2, Λmb, etc.
```

### 5.5 Accessor helpers

```julia
get_model(geom)                        # CartesianDiscreteModel
get_triangulation(geom, :Ω)           # Interior
get_triangulation(geom, :Γ_fs)        # Free surface
get_triangulation(geom, :membrane)    # Structure by label
get_measure(geom, :dΩ)                # Measure(Ω, degree)
element_size(geom)                     # h = LΩ/nx
```

---

## 6. Implementation phases

### Phase 1 — Core geometry builder (MVP)

Files: `Geometry.jl`, `TankDomain2D.jl`, `MeshGrading.jl`, `SurfaceLayout.jl`

Tasks:
1. Define `TankDomain2D` and `StructureZone` structs.
2. Implement `build_model(tank)` — constructs `CartesianDiscreteModel`
   with `f_y` grading map.
3. Implement `label_model!(model)` — applies the standard five labels
   (surface, bottom, inlet, outlet, water).
4. Implement surface masking: given a list of `StructureZone` and the
   computed x-positions, create `is_<zone>(xs)` mask closures and
   partition `Γ` into sub-triangulations.
5. Structural-boundary labelling: find structure edge vertices, create
   `"<label>_bnd"` tags.
6. Skeleton triangulation for DG beams (skip if no beam entity).
7. Return `TankGeometry` with all triangulations, measures, and normals.

Tests (Phase 1):
- Round-trip: `TankDomain2D → build_geometry → TankGeometry` with a
  single membrane produces the correct number of cells in each
  sub-triangulation.
- Label consistency: every surface cell belongs to exactly one of
  `Γfs ∪ Γd_in ∪ Γd_out ∪ ⋃ Γ_structures`.
- Measure creation: `sum(∫(1)dΩ) ≈ LΩ * H0`.
- Grading: first and last vertical element sizes differ by `mesh_ry^ny`.
- Skeleton normals exist and have unit length.

### Phase 2 — Physics bridge + damping

Files: `GeometryToPhysics.jl`, `DampingZone.jl`

Tasks:
1. `to_weakform_domains(geom)` — map `TankGeometry` fields to the
   `WeakFormDomains` key conventions expected by entities
   (`:dΩ`, `:dΓ_fs`, `:dΓ_s`, `:dΛ_s`, `:n_Λ_s`, `:h_s`,
   `:dΛ_sb`, `:n_Λ_sb`, `:dΓ_in`, `:dΓ_ot`).
2. `to_weakform_domains(geom, resonators)` — also add `:δ_p` from
   `ResonatorSingle.XZ` positions.
3. `DampingZone` struct + `damping_functions(tank; k)` returning
   μ₁ᵢₙ, μ₁ₒᵤₜ (and wave-number-scaled variants).

Tests (Phase 2):
- `to_weakform_domains` populates all mandatory keys.
- DiracDelta vector length matches resonator count.
- Damping coefficient μ equals 0 at inner edge and μ₀ at domain
  boundary.
- Assembled `WeakFormDomains` plugs into existing `WeakFormAssembly`
  integration tests without change.

### Phase 3 — Multi-structure & advanced layouts

Tasks:
1. Support multiple `StructureZone`s (e.g., two membranes, or a membrane
   followed by a beam) with gaps between them.
2. Each structure gets its own labelled sub-triangulation and
   skeleton/boundary, with its own key in `WeakFormDomains`.
3. Variable bathymetry: optional `bathy::Function` field on
   `TankDomain2D` that modifies the grading map in x as well as y
   (cf. `mem_freq_bathy_damp.jl`).
4. `write_geometry_vtk` utility.

Tests (Phase 3):
- Two structures: triangulations are disjoint and cover
  the correct x-ranges.
- Bathymetry: bottom boundary changes depth smoothly.
- VTK files are written without error.

---

## 7. Design decisions & rationale

| Decision | Rationale |
|---|---|
| Module sits **outside** Physics | Geometry is a *pre-processing* concern; Physics should stay physics-only. Prevents circular deps. |
| Struct-based config (`TankDomain2D`) instead of keyword args | Enables serialisation, testing, and re-use of identical geometry across frequency/time-domain runs. |
| `TankGeometry` stores triangulations **and** measures | Measures are cheap and always needed; avoids a second construction step. |
| Mapping to `WeakFormDomains` is a separate function | Keeps geometry independent of the Physics key conventions; easy to update if conventions change. |
| `StructureZone` is physics-agnostic (just a label + length) | The actual entity type (Membrane2D vs EulerBernoulliBeam) is determined elsewhere; geometry only needs extent. |
| `f_y` grading extracted into `MeshGrading.jl` | Allows future grading strategies (e.g., tanh, cosine, adaptive) without touching the rest. |

---

## 8. Migration path (existing scripts)

After Phase 2, a typical frequency-domain script reduces from ~100 lines
of geometry boilerplate to:

```julia
using HydroElasticFEM
using HydroElasticFEM.Geometry

tank = TankDomain2D(H0=10, Ld_in=150, Ld_out=150,
          L_free=80,
          structures=[StructureZone(:membrane, 20.0, 0.0)],
          nx=1000, ny=20)
geom = build_geometry(tank)
dom  = to_weakform_domains(geom)

# Physics + assembly unchanged
pf   = PotentialFlow()
fs   = FreeSurface()
mem  = Membrane2D(L=20.0, mᵨ=0.9, Tᵨ=98.1)
X, Y, fmap = build_fe_spaces(pf => geom.Ω, fs => geom.Γκ, mem => geom.structures[:membrane].Γs)
# … assemble and solve as before
```

All existing `WeakFormDomains` consumers remain unchanged.
