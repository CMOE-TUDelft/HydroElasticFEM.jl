# References

## Primary formulation references

### [C23] Colomes et al. (2023)
O. Colomes, F. Verdugo, and I. Akkerman,
"A monolithic finite element formulation for the hydroelastic analysis
of very large floating structures,"
*International Journal for Numerical Methods in Engineering*,
124(3):714-751, 2023.
DOI: [10.1002/nme.7140](https://doi.org/10.1002/nme.7140)

This is the primary reference for HydroElasticFEM.jl. It provides:
- The monolithic potential-flow formulation (Section 2)
- The Euler-Bernoulli beam C/DG weak form (Section 3.1)
- The Kirchhoff-Love plate C/DG weak form (Section 3.2)
- Fluid-structure coupling terms (Sections 3.1-3.2)
- Numerical examples and validation studies

### [A24] Agarwal et al. (2024)
S. Agarwal, O. Colomes, and A. V. Metrikine,
"Dynamic analysis of viscoelastic floating membranes using monolithic
finite element method,"
*Journal of Fluids and Structures*, 129:104167, 2024.
DOI: [10.1016/j.jfluidstructs.2024.104167](https://doi.org/10.1016/j.jfluidstructs.2024.104167)

This reference supports membrane-oriented formulation details and the
monolithic treatment of viscoelastic floating membranes used by the
`Membrane` implementation.

## Post-processing references

### Goda & Suzuki (1976)
Y. Goda and Y. Suzuki,
"Estimation of incident and reflected waves in random wave experiments,"
*Proceedings of the 15th International Conference on Coastal Engineering*,
828-845, 1976.

Basis for the three-gauge reflection/transmission separation method
generalized by [`HydroElasticFEM.PostProcessing.fit_wave_components`](@ref)
to an arbitrary number of probes.

### Zelt & Skjelbreia (1992)
J. A. Zelt and J. E. Skjelbreia,
"Estimating incident and reflected wave fields using an arbitrary number of
wave gauges,"
*Proceedings of the 23rd International Conference on Coastal Engineering*,
777-789, 1992.

Least-squares generalization of the Goda & Suzuki method to `N > 3` gauges,
used by the `PostProcessing` module.

## Software dependencies

- **Gridap.jl**: F. Verdugo and S. Badia, "The software design of Gridap:
  a finite element package based on the Julia JIT compiler,"
  *Computer Physics Communications*, 276:108341, 2022.
  DOI: [10.1016/j.cpc.2022.108341](https://doi.org/10.1016/j.cpc.2022.108341)
- **GridapGmsh.jl**: <https://github.com/gridap/GridapGmsh.jl>
- **Gmsh**: C. Geuzaine and J.-F. Remacle, "Gmsh: A 3-D finite element
  mesh generator with built-in pre- and post-processing facilities,"
  *International Journal for Numerical Methods in Engineering*, 79(11):1309-1331,
  2009. DOI: [10.1002/nme.2579](https://doi.org/10.1002/nme.2579)
