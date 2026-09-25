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

**Tensioned Euler–Bernoulli beam:**

`TensionedEulerBernoulliBeam` combines the two operators above — membrane
pre-tension and Euler–Bernoulli bending — in a single structure:

```math
m_\rho \ddot{\eta} + EI_\rho \nabla^4 \eta - \nabla\cdot(T_\rho \nabla \eta)
+ \text{coupling} = 0 \quad \text{on } \Gamma_\eta.
```

Setting `EIᵨ = 0` recovers the `Membrane` operator exactly; setting `Tᵨ = 0`
recovers the `EulerBernoulliBeam` operator exactly (both are verified at the
assembled-matrix level, to machine precision, in
`test/examples/TensionedEulerBernoulliBeamWeakFormTests.jl`). For a mode of
wavenumber $k$, the two contributions scale as $EI_\rho k^2$ (bending) and
$T_\rho$ (tension); the crossover wavenumber $k^\star = \sqrt{T_\rho / EI_\rho}$
separates a tension-dominated regime ($EI_\rho k^2 \ll T_\rho$, membrane-like)
from a bending-dominated regime ($EI_\rho k^2 \gg T_\rho$, beam-like).

### Resonators

Point-mass resonators are modelled as locally resonant mass–spring–damper
systems coupled to the free surface through Dirac-delta functionals.

## Frequency-domain discretisation

After Fourier transformation (``\partial_t \to -i\omega``), the system
reduces to a complex-valued linear FE problem.  The bilinear form is:

```math
a(u,v) = -\omega^2 m(u,v) - i\omega c(u,v) + k(u,v).
```

For `Membrane`, `EulerBernoulliBeam`, and `TensionedEulerBernoulliBeam`, the
mass form $m$, the hydrostatic part of the stiffness form $k$, and the
right-hand-side form share one common implementation
(`AbstractHydroelasticStructure`, see [How to Add a New Structural
Entity](@ref)), and the damping form $c$ is obtained automatically as the
structure's elastic stiffness operator scaled by its Rayleigh damping
coefficient $\tau$. Each structure only supplies its own elastic stiffness
operator — pre-tension for `Membrane`, bending for `EulerBernoulliBeam`, or
their sum for `TensionedEulerBernoulliBeam`.

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
