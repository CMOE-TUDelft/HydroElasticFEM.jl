# HydroElasticFEM.jl

**HydroElasticFEM.jl** is a Julia package for finite-element simulation of
hydro-elastic wave–structure interaction problems.  It couples potential-flow
fluid models with thin structural models (2-D membranes, Euler–Bernoulli beams,
resonators) using the Gridap.jl FE framework.

## Package features

- Frequency-domain and time-domain simulations.
- Modular physics entities: `PotentialFlow`, `FreeSurface`, `Membrane2D`,
  `EulerBernoulliBeam`, `ResonatorSingle`.
- Automated multi-field FE-space construction and operator assembly.
- Symbol-based field access via `FieldMap`.
- Cartesian mesh generation with structure and damping-zone support.
- Coupling terms detected automatically from physics entity traits.

## Quick navigation

| Section | Description |
|---------|-------------|
| [Getting Started](@ref) | Installation and first simulation |
| [Examples](@ref) | Worked example scripts |
| [Theory](@ref) | Governing equations and discretisation |
| [API Reference](@ref "API Overview") | Full docstring reference |

## Contents

```@contents
Pages = [
    "guide/getting_started.md",
    "guide/examples.md",
    "guide/theory.md",
    "api/index.md",
    "api/physics.md",
    "api/geometry.md",
    "api/simulation.md",
    "api/assembly_contexts.md",
    "api/parameter_handler.md",
]
Depth = 2
```
