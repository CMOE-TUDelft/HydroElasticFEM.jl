# HydroElasticFEM.jl

```@raw html
<div style="text-align: center; margin: 0 0 1.5rem;">
    <img src="assets/logo.svg" alt="HydroElasticFEM.jl logo" style="max-width: 360px; width: 100%; height: auto;">
</div>
```

**HydroElasticFEM.jl** is a Julia package for finite element (FE) simulation of
hydro-elastic wave–structure interaction problems.  It couples potential flow
fluid models with structural models (membranes, Euler–Bernoulli beams,
resonators) using the [Gridap.jl](https://www.github.com/gridap/Gridap.jl) FE framework.

## Package features

- Frequency-domain and time-domain simulations.
- Modular and exxpandable physics entities: `PotentialFlow`, `FreeSurface`, `Membrane2D`,
  `EulerBernoulliBeam`, `ResonatorSingle`.
- Automated multi-field FE space construction and operator assembly.
- Automatic differentiation for nonlinear problems.
- Cartesian mesh generation with structure and damping-zone support (unstructured mesh under development).
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
