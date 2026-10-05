# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `TensionedEulerBernoulliBeam`, a structural model combining membrane pre-tension stiffness with Euler-Bernoulli bending stiffness (governing equation `m·ηₜₜ + EI·Δ²η - ∇·(T∇η) = p`). Reduces exactly to `Membrane` when `EIᵨ = 0` and to `EulerBernoulliBeam` when `Tᵨ = 0`, verified at the assembled-matrix level. Supports the same `JointRotationalSpring` connections as `EulerBernoulliBeam`.
- `AbstractHydroelasticStructure`, a shared base type for `Membrane`, `EulerBernoulliBeam`, and `TensionedEulerBernoulliBeam` that factors out their common `mass`, hydrostatic `stiffness`, `rhs`, and stiffness-proportional Rayleigh `damping` weak forms. New structures following this pattern now only implement `mass_density`, `damping_parameter`, and `stiffness_operator` (see the "Shortcut for standard hydroelastic structures" note in the "Adding a New Structure" guide).
- `examples/FloatingTensionedBeamExample.jl`, comparing `Membrane`, `EulerBernoulliBeam`, and `TensionedEulerBernoulliBeam` under identical wave conditions.
- `get_integration_domains` stores `:dΓlateral` / `:nΓlateral` for the lateral walls of 3D tanks. Since [PR#59](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/59).
- `TimeConfig.scheme` selects the time integrator: `:generalized_alpha` (default, `ρ∞`) or `:newmark` (`γ`, `β`). Since [feat/time-scheme-newmark](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/tree/feat/time-scheme-newmark).
- Optional `mask` field in `PrescribedInletPotentialBC` to weight the forcing pointwise (e.g. apply it only on the generation part of a wall). Since Since [PR#60](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/60).
- `PrescribedInletPotentialBC(quantity = :velocity)`: prescribes the Neumann flux `∂ₙϕ = v_in ⋅ n` from an incident velocity vector, e.g. on 3D inlet and lateral walls. Since [PR#60](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/60).
- `JointLineDomain`, line joints on 3D plates: structure-skeleton facets on the given segments are split off into a joint skeleton with its own measure and normal, and removed from the C/DG skeletons. Misaligned or overlapping segments are rejected. Since [PR#61](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/61).
- `hinge_grid` returns the interior connection lines of an `nfx × nfy` floater array as a `JointLineDomain`. Since [PR#61](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/61).
- `AbstractJointDomain`, the common supertype of `JointDomain` and `JointLineDomain`. [PR#61](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/61).
- `KirchhoffLovePlate.joints`: rotational springs (`JointRotationalSpring`, `kᵣ` per unit joint length divided by ρ) on plate line joints declared with a `JointLineDomain`. A free hinge only needs the `JointLineDomain`: its facets are excluded from the C/DG slope-continuity terms. Since [PR#61](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/61).
- `get_integration_domains` also stores each structure's C/DG skeleton data under its `domain_symbol` (`:dΛ_<sym>`, `:n_Λ_<sym>`, `:h_<sym>`), and `skeleton_keys(sym)` returns these keys. Since [PR#62](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/62).
- `time_integration_parameters` and `stabilization_αₕ` helpers in `ParameterHandler`; `stabilization_αₕ` returns `γ/(βΔt)/g · (1-βₕ)/βₕ` for the configured scheme. Since [PR#63](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/63).
- Spike test for coupling two separate plates across a shared edge (`test/Geometry/PlateInterfaceSkeletonSpikeTests.jl`): an interface `SkeletonTriangulation` with the plus side on one plate and the minus side on the other works for traces, gradients, normals and multi-field assembly (go decision for plate connections). Since [PR#68](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/68).
- `StructureConnection` and `TankDomain(structure_connections = …)`: the interface skeleton between two structures that have separate fields, with the plus side always in structure `a` and `n⁺` pointing from `a` to `b`; its measure and normal are stored under the connection's `domain_symbol` / `normal_symbol`. Since [PR#69](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/69).
- `get_integration_domains` also stores each structure's C/DG skeleton data under its `domain_symbol` (`:dΛ_<sym>`, `:n_Λ_<sym>`, `:h_<sym>`), and `skeleton_keys(sym)` returns these keys. Since [PR#71](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/71).
- `PlateConnection`, a field-less coupling entity between two `KirchhoffLovePlate`s with separate fields across a `StructureConnection`: shear as a spring or rigid (penalty), and rotation as rigid (C/DG slope continuity), free (hinge) or a rotational spring. Since [PR#71](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/71).
- `IncidentSea`, a linear Airy sea whose elevation, potential and velocity fields can drive time-domain boundary conditions. Forcing from one sea is precomputed per frequency and reused at each time step.

### Changed
- CI runs the test suite at `-O1` (`Pkg.test(julia_args = ["-O1"])`). Almost all test time was LLVM optimisation of Gridap's generated code at the default `-O2`: the full suite drops from ~38 min to ~5 min locally with identical results. The duplicated Khabakhpasheva run in `TensionedEulerBernoulliBeamWeakFormTests.jl` was removed (the same case is checked in `KhabakhpashevaBeamJointTests.jl`).
- `Membrane` and `EulerBernoulliBeam` now subtype `AbstractHydroelasticStructure` instead of `Structure` directly. Their assembled weak forms are unchanged; `mass`, `damping`, `stiffness`, and `rhs` are now provided by the shared base type instead of being implemented per-structure. `EulerBernoulliBeam`'s C/DG bending operator and joint-penalty assembly are now shared helper functions (`_eb_bending_stiffness_operator`, `_joint_stiffness_form`), reused by `TensionedEulerBernoulliBeam`.
- `RadiationBC` now supports both frequency-domain and time-domain assembly contexts, fixing issue [#44](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/issues/44). Since [PR#46](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/46).
- Potential-flow boundary conditions fall back to the normal of the measure's own triangulation when no `:nΓ…` key is stored (e.g. `:dΓin`, `:dΓout`). Since [PR#60](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/60).
- `TankDomain{3}` accepts `JointLineDomain` joints; `JointDomain` remains 2D only. Since [PR#61](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/61).
- `KirchhoffLovePlate` reads its skeleton measure, normal and element size from its own `space_domain_symbol` (global `:dΛη`/`:n_Λ_η`/`:h_η` for `:Γη`, unchanged), so several plates with separate fields each get only their own interior facets. Since [PR#62](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/62).
- `TimeConfig.αₕ = nothing` is now computed automatically for time-domain problems with damping zones (previously an error); `αₕ = :auto` opts in for other problems. Problems without damping zones are unchanged. Since [PR#63](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/63).
- `TimeConfig.u0 = nothing` starts the simulation from rest (all initial states zero). Since [PR#63](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/63).
- Space/time forcing functions are resolved once per assembly instead of wrapping every quadrature-point evaluation in `try`/`catch`, which reduces the per-step cost of time-dependent forcing. Behaviour is unchanged. Since [PR#64](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/64).
- Zero body forcing (the default when no `rhs_fn` is given) no longer assembles `∫ v·0` volume and structure integrals at every step, and a problem without any forcing now gets an empty right-hand side instead of an error. Since [PR#65](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/65).
- `ResonatorArray` self and resonator↔structure forms no longer add a dummy `∫(ξ⋅q·0)` integral over the host surface; they are now pure point (Dirac) terms. Since [PR#66](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/66).
- Time-domain linear residuals rely on Gridap's constant-form optimization; HydroElasticFEM only adds a residual adapter when it can replace an `IncidentSea` forcing with precomputed vectors.
- `KirchhoffLovePlate` reads its skeleton measure, normal and element size from its own `space_domain_symbol` (global `:dΛη`/`:n_Λ_η`/`:h_η` for `:Γη`, unchanged), so several plates with separate fields each get only their own interior facets. Since [PR#71](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/71).
- `build_fe_spaces` skips physics entities without fields of their own (e.g. `PlateConnection`). Since [PR#71](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/71).
- Resonator ↔ structure coupling is restricted to the resonators whose `host_domain_symbol` equals the structure's `space_domain_symbol` (`attached_resonators`). With several structures that have separate fields, each resonator couples only to its own structure, and one `ResonatorArray` may span several structures (mixed hosts are now allowed). Since [PR#73](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/73).

### Fixed
- Integration domains not covered by a physics entity (C/DG skeletons such as `:Λη`, damping zones, `:Γlateral`, joints) now use the highest entity quadrature degree instead of 2, so higher-order plate skeletons are no longer under-integrated. Derived measure keys (`:dΓd_i`, `:dΛη_i`, joint measures) use the same value via a new `:default` entry of the degree dictionary instead of a hard-coded 4. Since [PR#59](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/59).
- `build_fe_spaces` builds the resonator `ConstantFESpace`s on the volume `:Ω`. With the default host `:Γη`, assembling resonator terms failed with a Gridap type-instability error. [PR#66](https://github.com/CMOE-TUDelft/HydroElasticFEM.jl/pull/66).

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
