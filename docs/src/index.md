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
- Modular and expandable physics entities: `PotentialFlow`, `FreeSurface`, `Membrane`,
  `EulerBernoulliBeam`, `ResonatorSingle`.
- Automated multi-field FE space construction and operator assembly.
- Automatic differentiation for nonlinear problems.
- Cartesian mesh generation with structure and damping-zone support (unstructured mesh under development).
- Coupling terms detected automatically from physics entity traits.

## Quick navigation

| Section | Description |
|---------|-------------|
| [Getting Started](@ref) | Installation and first simulation |
| [First Simulation](@ref "Run Your First Hydroelastic Simulation") | Step-by-step walkthrough: empty tank and floating membrane |
| [Examples](@ref) | Worked example scripts |
| [Theory](@ref) | Governing equations and discretisation |
| [Architecture Guide](@ref "Architecture Guide") | Three-layer design and extension points |
| [Adding a New Structural Entity](@ref "How to Add a New Structural Entity") | Developer recipe for new physics entities |
| [Debugging](@ref "Debugging Guide: When Simulations Go Wrong") | Error messages, physics checks, convergence tips |
| [API Reference](@ref "API Overview") | Full docstring reference |
| [References](references.md) | Primary literature and software dependencies |

## Contents

```@contents
Pages = [
    "references.md",
    "guide/index.md",
    "guide/getting_started.md",
    "guide/first_simulation.md",
    "guide/examples.md",
    "guide/theory.md",
    "guide/architecture.md",
    "guide/adding_structure.md",
    "guide/debugging.md",
    "api/index.md",
    "api/physics.md",
    "api/geometry.md",
    "api/simulation.md",
    "api/assembly_contexts.md",
    "api/parameter_handler.md",
]
Depth = 2
```
