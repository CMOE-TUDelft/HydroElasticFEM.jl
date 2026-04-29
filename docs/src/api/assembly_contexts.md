# Assembly Contexts

Assembly contexts carry analysis-dependent metadata (domains, frequency, time,
stabilisation parameter) and are passed to all weak-form methods.

## Types

```@docs
HydroElasticFEM.AssemblyContexts.AbstractAssemblyContext
HydroElasticFEM.AssemblyContexts.FrequencyAssemblyContext
HydroElasticFEM.AssemblyContexts.TimeAssemblyContext
```

## Accessors

```@docs
HydroElasticFEM.AssemblyContexts.domains
HydroElasticFEM.AssemblyContexts.frequency
HydroElasticFEM.AssemblyContexts.current_time
HydroElasticFEM.AssemblyContexts.stabilization_parameter
HydroElasticFEM.AssemblyContexts.has_stabilization
HydroElasticFEM.AssemblyContexts.with_time
```
