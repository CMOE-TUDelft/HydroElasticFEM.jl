# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `TensionedEulerBernoulliBeam`, a structural model combining membrane
  pre-tension stiffness with Euler-Bernoulli bending stiffness (governing
  equation `m·ηₜₜ + EI·Δ²η - ∇·(T∇η) = p`). Reduces exactly to `Membrane`
  when `EIᵨ = 0` and to `EulerBernoulliBeam` when `Tᵨ = 0`, verified at the
  assembled-matrix level. Supports the same `JointRotationalSpring`
  connections as `EulerBernoulliBeam`.
- `AbstractHydroelasticStructure`, a shared base type for `Membrane`,
  `EulerBernoulliBeam`, and `TensionedEulerBernoulliBeam` that factors out
  their common `mass`, hydrostatic `stiffness`, `rhs`, and
  stiffness-proportional Rayleigh `damping` weak forms. New structures
  following this pattern now only implement `mass_density`,
  `damping_parameter`, and `stiffness_operator` (see the "Shortcut for
  standard hydroelastic structures" note in the "Adding a New Structure"
  guide).
- `examples/FloatingTensionedBeamExample.jl`, comparing `Membrane`,
  `EulerBernoulliBeam`, and `TensionedEulerBernoulliBeam` under identical
  wave conditions.

### Changed
- `Membrane` and `EulerBernoulliBeam` now subtype `AbstractHydroelasticStructure`
  instead of `Structure` directly. Their assembled weak forms are unchanged;
  `mass`, `damping`, `stiffness`, and `rhs` are now provided by the shared
  base type instead of being implemented per-structure. `EulerBernoulliBeam`'s
  C/DG bending operator and joint-penalty assembly are now shared helper
  functions (`_eb_bending_stiffness_operator`, `_joint_stiffness_form`),
  reused by `TensionedEulerBernoulliBeam`.
- `RadiationBC` now supports both frequency-domain and time-domain assembly contexts, fixing issue [#44](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/issues/44). Since [PR#46](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/46).

### Fixed

## [0.1.2 - 2026-06-17]

### Added

### Changed

### Fixed
- Enabled normal vectors of user-defined domains. Since [PR#40](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/40)

## [0.1.1] - 2026-06-05

### Added
- `GmshDomain` type for hydroelastic simulations on unstructured Gmsh meshes. Boundaries are identified by physical-group names; damping zones are detected automatically from `"damping_in"` / `"damping_out"` groups. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `CartesianDomain{D}` and `TankDomain{D}` are now dimension-generic, replacing the old `TankDomain2D` / `TankDomain3D` specialised types. 3-D sloshing and hydroelastic frequency-domain simulations are now supported. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- Liu (2004) floating-beam benchmark on an unstructured Gmsh mesh (`examples/LiuBenchmarkGmsh.jl`). Since [PR#16](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/16).
- `KirchhoffLovePlate` struct and SIPG weak forms (`mass`, `stiffness`, `rhs`) in `src/Physics/Structures/KirchhoffLovePlate.jl`, with validation test against the Timoshenko simply-supported square plate reference. Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).
- `build_kl_tensor` / `build_KL_tensor` helpers for the KL constitutive fourth-order tensor. Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).
- 3-D frequency-domain sloshing example (`examples/YagoBenchmark3DFreq.jl`). Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).
- `TimoshenkoBeam` struct with two-field (`w`, `θ`) formulation and mixed-order interpolation (order 2 / order 1) for shear-locking-free behaviour. Since [PR#18](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/18).
- Bibliography and academic references added to the documentation. Since [PR#23](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/23).
- `is_periodic` flag support in `CartesianDomain` and `TankDomain` to enable periodic boundary conditions along the horizontal direction. Since [PR#27](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/27).

### Changed
- `Membrane2D` renamed to `Membrane` throughout `src/`, `test/`, and `docs/`. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `src/Geometry/` rewritten around an `AbstractDomain` interface that unifies `TankDomain` and `GmshDomain` under a single `build_model` → `build_triangulations` → `get_integration_domains` pipeline. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `src/Physics/` reorganised: fluid entities (`PotentialFlow`, `FreeSurface`) moved to `src/Physics/Fluid/`; structural entities (`EulerBernoulliBeam`, `Membrane`, `Resonator`, `KirchhoffLovePlate`) moved to `src/Physics/Structures/`. The `Plate/` sub-subfolder is removed; `KirchhoffLovePlate` sits directly in `Structures/`. Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).
- `FESpaceAssembly.build_fe_spaces` extended with `variable_symbols` / `field_fe_configs` protocol to support multi-field entities (e.g. `TimoshenkoBeam`). Single-field entities are unaffected (backward-compatible default implementations provided in `Physics.jl`). Since [PR#18](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/18).
- Improved documentation: expanded content, better structure, and corrected repository URLs. Since [PR#19](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/19), [PR#21](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/21), [PR#22](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/22).
- CI workflow updated: documentation now built and deployed via CI pipeline; `julia-actions/setup-julia` bumped from v2 to v3. Since [PR#20](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/20), [PR#24](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/24).

### Fixed
- Hardcoded `src/Physics/PotentialFlow.jl` and `src/Physics/FreeSurface.jl` paths in `test/Simulation/SimulationTests.jl` updated to reflect the new `Fluid/` subfolder location.

## [0.1.0] - 2026-05-03

### Changed
- Previous changes were not tracked in the CHANGELOG.md file. Please, see commit history.
