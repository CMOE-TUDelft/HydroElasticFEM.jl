# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `GmshDomain` type for hydroelastic simulations on unstructured Gmsh meshes. Boundaries are identified by physical-group names; damping zones are detected automatically from `"damping_in"` / `"damping_out"` groups. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `CartesianDomain{D}` and `TankDomain{D}` are now dimension-generic, replacing the old `TankDomain2D` / `TankDomain3D` specialised types. 3-D sloshing and hydroelastic frequency-domain simulations are now supported. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- Liu (2004) floating-beam benchmark on an unstructured Gmsh mesh (`examples/LiuBenchmarkGmsh.jl`). Since [PR#16](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/16).
- `KirchhoffLovePlate` struct and SIPG weak forms (`mass`, `stiffness`, `rhs`) in `src/Physics/Structures/KirchhoffLovePlate.jl`, with validation test against the Timoshenko simply-supported square plate reference. Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).
- `build_kl_tensor` / `build_KL_tensor` helpers for the KL constitutive fourth-order tensor. Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).
- 3-D frequency-domain sloshing example (`examples/YagoBenchmark3DFreq.jl`). Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).

### Changed
- `src/Physics/` reorganised: fluid entities (`PotentialFlow`, `FreeSurface`) moved to `src/Physics/Fluid/`; structural entities (`EulerBernoulliBeam`, `Membrane`, `Resonator`, `KirchhoffLovePlate`) moved to `src/Physics/Structures/`. The `Plate/` sub-subfolder is removed; `KirchhoffLovePlate` sits directly in `Structures/`. Since [PR#17](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/17).
- `src/Geometry/` rewritten around an `AbstractDomain` interface that unifies `TankDomain` and `GmshDomain` under a single `build_model` → `build_triangulations` → `get_integration_domains` pipeline. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `Membrane2D` renamed to `Membrane` throughout `src/`, `test/`, and `docs/`. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).

### Fixed
- Hardcoded `src/Physics/PotentialFlow.jl` and `src/Physics/FreeSurface.jl` paths in `test/Simulation/SimulationTests.jl` updated to reflect the new `Fluid/` subfolder location.

## [0.1.0] - 2026-05-03

### Changed
- Previous changes were not tracked in the CHANGELOG.md file. Please, see commit history.
