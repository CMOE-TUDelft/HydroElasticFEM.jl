# Getting Started

## Installation

HydroElasticFEM.jl is not yet registered in the General registry.
Install it directly from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/CMOE/HydroElasticFEM.jl")
```

The package requires Julia ≥ 1.10 and depends on
[Gridap.jl](https://github.com/gridap/Gridap.jl) and
[WaveSpec.jl](https://github.com/CMOE/WaveSpec.jl).

## Minimal frequency-domain example

```julia
using HydroElasticFEM

# 1. Define geometry
domain = TankDomain2D(
    Lx = 10.0, Lz = 1.0,
    nx = 40,   nz = 10,
)

# 2. Define physics entities
fluid  = PotentialFlow(ρw = 1025.0, g = 9.81)
fsurf  = FreeSurface(ρw = 1025.0, g = 9.81)

# 3. Configure simulation
config = FreqDomainConfig(ω = 1.0)

# 4. Build and solve
problem = build_problem(domain, [fluid, fsurf], config)
result  = simulate(problem)
```

## Time-domain example

```julia
tconfig = TimeConfig(Δt = 0.01, tf = 10.0,
                     u0 = (x -> 0.0, x -> 0.0))
tdconfig = TimeDomainConfig()
problem  = build_problem(domain, [fluid, fsurf], tdconfig; tconfig)
result   = simulate(problem, tconfig)
```

## Package layout

| Module | Purpose |
|--------|---------|
| `HydroElasticFEM.Physics` | Physics entity types and weak forms |
| `HydroElasticFEM.Geometry` | Mesh and integration domain helpers |
| `HydroElasticFEM.Simulation` | FE space assembly and time/frequency solvers |
| `HydroElasticFEM.AssemblyContexts` | Context objects for assembly |
| `HydroElasticFEM.ParameterHandler` | Configuration structs |
