# Physics

The `Physics` sub-module defines all physical-entity types together with their
weak-form methods (`mass`, `damping`, `stiffness`, `rhs`) and the composed
frequency/time-domain forms (`weakform`, `residual`, `jacobian`, …).

```@docs
HydroElasticFEM.Physics
```

## Abstract base

```@docs
HydroElasticFEM.Physics.PhysicsParameters
HydroElasticFEM.Physics.print_parameters
HydroElasticFEM.Physics.variable_symbol
HydroElasticFEM.Physics.variable_symbols
HydroElasticFEM.Physics.field_fe_configs
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

## Hydroelastic structure abstraction

`Membrane`, `EulerBernoulliBeam`, and `TensionedEulerBernoulliBeam` share one
implementation of `mass`, `damping`, `stiffness`, and `rhs` via
`AbstractHydroelasticStructure`; each only implements `mass_density`,
`damping_parameter`, and `stiffness_operator` (and, optionally,
`extra_stiffness_form` for contributions not subject to Rayleigh damping,
such as `EulerBernoulliBeam`/`TensionedEulerBernoulliBeam` joints). See [How
to Add a New Structural Entity](@ref) for a worked example.

```@docs
HydroElasticFEM.Physics.AbstractHydroelasticStructure
HydroElasticFEM.Physics.mass_density
HydroElasticFEM.Physics.damping_parameter
HydroElasticFEM.Physics.gravitational_acceleration
HydroElasticFEM.Physics.stiffness_operator
HydroElasticFEM.Physics.extra_stiffness_form
```

## Structures

```@docs
HydroElasticFEM.Physics.Membrane
HydroElasticFEM.Physics.EulerBernoulliBeam
HydroElasticFEM.Physics.JointRotationalSpring
HydroElasticFEM.Physics.TensionedEulerBernoulliBeam
HydroElasticFEM.Physics.KirchhoffLovePlate
HydroElasticFEM.Physics.TimoshenkoBeam
HydroElasticFEM.Physics.build_kl_tensor
HydroElasticFEM.Physics.build_KL_tensor
HydroElasticFEM.Physics.check_major_symmetry
HydroElasticFEM.Physics.equivalent_beam_rigidity
```

## Resonators

```@docs
HydroElasticFEM.Physics.ResonatorSingle
HydroElasticFEM.Physics.resonator_array
```
