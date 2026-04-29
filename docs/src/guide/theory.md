# Theory

## Governing equations

HydroElasticFEM solves the linearised water-wave/structure interaction problem
in either the frequency domain (harmonic steady state) or the time domain.

### Fluid: potential flow

The fluid is governed by Laplace's equation for the velocity potential $\phi$:

```math
\nabla^2 \phi = 0 \quad \text{in } \Omega,
```

subject to the linearised free-surface condition on $\Gamma_\kappa$:

```math
\frac{\partial \phi}{\partial t} + g\,\kappa = 0, \qquad
\frac{\partial \kappa}{\partial t} = \frac{\partial \phi}{\partial z},
```

where $\kappa$ is the free-surface elevation auxiliary variable.

### Structures

**Membrane (1D in a 2D domain):**

```math
m_\rho \ddot{\eta} - T_\rho \nabla^2 \eta + \text{coupling} = 0
\quad \text{on } \Gamma_\eta.
```

**Euler–Bernoulli beam:**

```math
m_\rho \ddot{\eta} + EI_\rho \nabla^4 \eta + \text{coupling} = 0
\quad \text{on } \Gamma_\eta.
```

### Resonators

Point-mass resonators are modelled as locally resonant mass–spring–damper
systems coupled to the free surface through Dirac-delta functionals.

## Frequency-domain discretisation

After Fourier transformation (``\partial_t \to -i\omega``), the system
reduces to a complex-valued linear FE problem.  The bilinear form is:

```math
a(u,v) = -\omega^2 m(u,v) - i\omega c(u,v) + k(u,v).
```

## Time-domain discretisation

Time integration uses the **Generalised-``\alpha``** method (second-order,
unconditionally stable for ``\rho_\infty \in [0,1]``) as implemented by
`Gridap.ODEs.GeneralizedAlpha2`.

## Stabilisation

A stabilised free-surface formulation is optionally available, controlled by
``\beta_h \in (0,1]`` in `FreeSurface` and the derived ``\alpha_h`` parameter
stored in the assembly context.

## References

- Colomes, Agarwal et al. — *HydroElasticFEM: a Gridap-based FE solver for
  wave–structure interaction* (in preparation).
- Donéa & Huerta, *Finite Element Methods for Flow Problems*, Wiley, 2003.
