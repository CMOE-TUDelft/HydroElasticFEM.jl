# PostProcessing

The `PostProcessing` sub-module computes standard postprocessing functions. 
For example, determination of reflection, transmission, and absorption coefficients 
(and related energy quantities) from a known, single-frequency free-surface 
elevation field. It has no dependency on `Physics`/`Simulation` internals — 
it operates on whatever field and probe points you give it, whether that came 
from `result.solution`, reloaded VTK output, or another solver entirely.

See the [Reflection, Transmission & Absorption](@ref) guide for a full
worked example on a `TankDomain`.

```@docs
HydroElasticFEM.PostProcessing
```

## Wave decomposition (physics-agnostic building blocks)

```@docs
HydroElasticFEM.PostProcessing.WaveProbeFit
HydroElasticFEM.PostProcessing.fit_wave_components
HydroElasticFEM.PostProcessing.suggest_probe_offsets
HydroElasticFEM.PostProcessing.probe_positions
HydroElasticFEM.PostProcessing.sample_probe_line
```

## Reflection, transmission & absorption

```@docs
HydroElasticFEM.PostProcessing.ReflectionTransmissionResult
HydroElasticFEM.PostProcessing.reflection_transmission_coefficients
HydroElasticFEM.PostProcessing.absorption_coefficient
```

## Energy flux and dissipation-based absorption

```@docs
HydroElasticFEM.PostProcessing.group_velocity
HydroElasticFEM.PostProcessing.energy_flux
HydroElasticFEM.PostProcessing.absorption_from_dissipated_power
HydroElasticFEM.PostProcessing.resonator_dissipated_power
```
