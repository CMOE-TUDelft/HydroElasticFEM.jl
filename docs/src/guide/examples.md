# Examples

Examples are located in the `examples/` directory of the repository.

## Empty tank (frequency domain)

`examples/EmptyTankExample.jl` — sets up a rectangular fluid domain with a
free surface and runs a single-frequency potential-flow simulation.

## Empty tank with time-domain damping

`examples/EmptyTankTimeDomainDampingExample.jl` — extends the frequency-domain
tank to a time-domain problem using Generalized-α time integration and a
damping-zone radiation condition.

## Floating membrane

`examples/FloatingMembraneExample.jl` — couples a `PotentialFlow` fluid with a
`Membrane2D` structural entity and a `FreeSurface` free-surface condition.
Demonstrates how coupling terms are detected automatically.

## Euler-Bernoulli beam with joints

`scripts/EulerBernoulliBeam_example.jl` — standalone beam problem using the
DG weak form.  For an example of a **beam with rotational-spring joints**,
see the [Setting up structural joints](@ref) section of the Geometry API
reference, which walks through the end-to-end setup:

1. Declare a `JointDomain1D` (location + symbol keys) inside `TankDomain2D`.
2. Call `build_triangulations` and `get_integration_domains` — the joint
   skeleton facet and its measure/normal are built automatically.
3. Pass a `JointRotationalSpring` (matching symbols + `kᵣ`) to
   `EulerBernoulliBeam.joints`.
4. Build and solve with `build_problem` / `simulate` as usual.

## Running the examples

```julia
using HydroElasticFEM
include("examples/EmptyTankExample.jl")
```

Each example is self-contained and can be run from the Julia REPL or via
`julia --project examples/<ExampleFile>.jl`.
