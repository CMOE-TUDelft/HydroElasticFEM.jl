# Physics

The `Physics` sub-module defines all physical-entity types together with their
weak-form methods (`mass`, `damping`, `stiffness`, `rhs`) and the composed
frequency/time-domain forms (`weakform`, `residual`, `jacobian`, …).

## Abstract base

```@docs
HydroElasticFEM.Physics.PhysicsParameters
HydroElasticFEM.Physics.print_parameters
HydroElasticFEM.Physics.variable_symbol
```

## Weak-form interface

```@docs
HydroElasticFEM.Physics.mass
HydroElasticFEM.Physics.damping
HydroElasticFEM.Physics.stiffness
HydroElasticFEM.Physics.rhs
HydroElasticFEM.Physics.weakform
HydroElasticFEM.Physics.residual
HydroElasticFEM.Physics.jacobian
HydroElasticFEM.Physics.jacobian_t
HydroElasticFEM.Physics.jacobian_tt
```

## Form-presence traits

```@docs
HydroElasticFEM.Physics.has_mass_form
HydroElasticFEM.Physics.has_damping_form
HydroElasticFEM.Physics.has_stiffness_form
HydroElasticFEM.Physics.has_rhs_form
HydroElasticFEM.Physics.active_forms
```

## Fluid: potential flow

```@docs
HydroElasticFEM.Physics.PotentialFlow
HydroElasticFEM.Physics.AbstractPotentialFlowBC
HydroElasticFEM.Physics.RadiationBC
HydroElasticFEM.Physics.PrescribedInletPotentialBC
HydroElasticFEM.Physics.DampingZoneBC
```

## Free surface

```@docs
HydroElasticFEM.Physics.FreeSurface
```

## Structures

```@docs
HydroElasticFEM.Physics.Membrane2D
HydroElasticFEM.Physics.EulerBernoulliBeam
```

## Resonators

```@docs
HydroElasticFEM.Physics.ResonatorSingle
HydroElasticFEM.Physics.resonator_array
```
