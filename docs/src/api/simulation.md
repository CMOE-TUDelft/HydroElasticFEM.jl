# Simulation

The `Simulation` sub-module orchestrates FE-space construction, operator
assembly, and time/frequency-domain solving.

## Entry points

```@docs
HydroElasticFEM.Simulation.build_problem
HydroElasticFEM.Simulation.simulate
HydroElasticFEM.Simulation.SimResult
```

## FE-space assembly

```@docs
HydroElasticFEM.Simulation.FESpaceAssembly.build_fe_spaces
HydroElasticFEM.Simulation.FESpaceAssembly.build_test_fe_space
HydroElasticFEM.Simulation.FESpaceAssembly.build_trial_fe_space
```

## FE-operator construction

```@docs
HydroElasticFEM.Simulation.FEOperators.FieldMap
HydroElasticFEM.Simulation.FEOperators.detect_couplings
HydroElasticFEM.Simulation.FEOperators.build_fe_operator
HydroElasticFEM.Simulation.FEOperators.build_frequency_fe_operator
HydroElasticFEM.Simulation.FEOperators.build_time_fe_operator
```

## Assembly context builders

```@docs
HydroElasticFEM.Simulation.build_frequency_context
HydroElasticFEM.Simulation.build_time_context
HydroElasticFEM.Simulation.get_assembly_context
```

## Low-level form assemblers

```@docs
HydroElasticFEM.Simulation.FEOperators.assemble_weakform
HydroElasticFEM.Simulation.FEOperators.assemble_mass
HydroElasticFEM.Simulation.FEOperators.assemble_damping
HydroElasticFEM.Simulation.FEOperators.assemble_stiffness
HydroElasticFEM.Simulation.FEOperators.assemble_rhs
HydroElasticFEM.Simulation.FEOperators.assemble_residual
HydroElasticFEM.Simulation.FEOperators.assemble_jacobian
HydroElasticFEM.Simulation.FEOperators.assemble_jacobian_t
HydroElasticFEM.Simulation.FEOperators.assemble_jacobian_tt
```
