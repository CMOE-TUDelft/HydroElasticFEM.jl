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

## Running the examples

```julia
using HydroElasticFEM
include("examples/EmptyTankExample.jl")
```

Each example is self-contained and can be run from the Julia REPL or via
`julia --project examples/<ExampleFile>.jl`.
