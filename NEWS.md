# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `GmshDomain` type for hydroelastic simulations on unstructured Gmsh meshes. Boundaries are identified by physical-group names; damping zones are detected automatically from `"damping_in"` / `"damping_out"` groups. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `CartesianDomain{D}` and `TankDomain{D}` are now dimension-generic, replacing the old `TankDomain2D` / `TankDomain3D` specialised types. 3-D sloshing and hydroelastic frequency-domain simulations are now supported. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `KirchhoffLovePlate` weak form in `src/Physics/Structures/Plate/`. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- Liu (2004) floating-beam benchmark on an unstructured Gmsh mesh (`examples/LiuBenchmarkGmsh.jl`). Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- 3-D frequency-domain sloshing example (`examples/YagoBenchmark3DFreq.jl`). Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).

### Changed
- `src/Geometry/` rewritten around an `AbstractDomain` interface that unifies `TankDomain` and `GmshDomain` under a single `build_model` → `build_triangulations` → `get_integration_domains` pipeline. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).
- `Membrane2D` renamed to `Membrane` throughout `src/`, `test/`, and `docs/`. Since [feature/gmsh_geometry](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feature/gmsh_geometry).

### Fixed
-

## [0.1.0] - 2026-05-03

### Changed
- Previous changes were not tracked in the CHANGELOG.md file. Please, see commit history.
