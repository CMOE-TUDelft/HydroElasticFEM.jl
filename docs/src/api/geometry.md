# Geometry

The `Geometry` sub-module provides domain specification types, Cartesian mesh
construction, and the integration-domain containers used by weak forms.

```@docs
HydroElasticFEM.Geometry
```

## Integration domains

```@docs
HydroElasticFEM.Geometry.STANDARD_TAGS
HydroElasticFEM.Geometry.AbstractDomain
HydroElasticFEM.Geometry.IntegrationDomains
HydroElasticFEM.Geometry.TankTriangulations
HydroElasticFEM.Geometry.ambient_dimension
HydroElasticFEM.Geometry.manifold_dimension
HydroElasticFEM.Geometry.boundary_tags
HydroElasticFEM.Geometry.triangulation
HydroElasticFEM.Geometry.get_boundary
```

## Cartesian geometry

```@docs
HydroElasticFEM.Geometry.CartesianDomain
HydroElasticFEM.Geometry.TankDomain
HydroElasticFEM.Geometry.StructureDomain
HydroElasticFEM.Geometry.DampingZone
HydroElasticFEM.Geometry.JointDomain
HydroElasticFEM.Geometry.build_model
HydroElasticFEM.Geometry.build_triangulations
HydroElasticFEM.Geometry.get_integration_domains
HydroElasticFEM.Geometry.surface_mask
HydroElasticFEM.Geometry.surface_masks
HydroElasticFEM.Geometry.joint_mask
```

## 3D Cartesian helpers

```@docs
HydroElasticFEM.Geometry.f_z
HydroElasticFEM.Geometry.map_fn
HydroElasticFEM.Geometry.get_plate_triangulation
```

## Gmsh geometry

```@docs
HydroElasticFEM.Geometry.GmshDomain
HydroElasticFEM.Geometry.validate_gmsh_tags
```

## Setting up structural joints

Joints are declared at the geometry level using `JointDomain` inside a
2D `TankDomain`, then automatically converted into
skeleton triangulations and integration-domain keys by
`build_triangulations` and `get_integration_domains`.

On the physics side, each joint maps one-to-one to a `JointRotationalSpring`
attached to the `EulerBernoulliBeam` via the same symbol keys.

```julia
using HydroElasticFEM

# ── 1. Geometry: beam at y = 1, x ∈ [0, 4] with a joint at x = 2 ──
beam_dom = StructureDomain(L=4.0, x₀=[0.0, 1.0])
joint    = JointDomain(
    location      = [2.0, 1.0],   # joint position in 2D (must lie on the beam)
    domain_symbol = :dΛj_1,       # key for the skeleton Measure
    normal_symbol = :n_Λ_j_1,     # key for the skeleton normal
)

tank = TankDomain(
    L=8.0, H=1.0, nx=80, ny=4,
    structure_domains = [beam_dom],
    joint_domains     = [joint],     # ← declare joint here
)

model  = build_model(tank)
trians = build_triangulations(tank, model)  # joint skeleton facet extracted
dom    = get_integration_domains(trians)    # :dΛj_1 and :n_Λ_j_1 now in dom

# ── 2. Physics: beam with a rotational spring at the same symbols ──
kᵣ   = 1.0e4   # rotational stiffness [N·m² / ρw]
beam = EulerBernoulliBeam(
    L    = 4.0,
    mᵨ   = 0.5,
    EIᵨ  = 200.0,
    joints = [JointRotationalSpring(:dΛj_1, :n_Λ_j_1, kᵣ)],
)

# ── 3. Simulate ──
config  = FreqDomainConfig(ω = 1.2)
problem = build_problem(tank, [beam], config)
result  = simulate(problem)
```

!!! tip
    The `domain_symbol`/`normal_symbol` pair in `JointDomain` **must match**
    those in the corresponding `JointRotationalSpring`.  Any mismatch will
    cause a `KeyError` at assembly time.

!!! note
    Multiple joints are supported — add one `JointDomain` per joint location
    to `joint_domains` and one `JointRotationalSpring` per joint to
    `EulerBernoulliBeam.joints`.  Each pair must use unique symbol names.
