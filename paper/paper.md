---
title: 'HydroElasticFEM.jl: A Julia package for monolithic finite element
        analysis of hydroelastic wave-structure interaction'
tags:
  - Julia
  - finite element method
  - hydroelasticity
  - wave-structure interaction
  - very large floating structures
  - potential flow
authors:
  - name: Oriol Colomés
    orcid: 0000-XXXX-XXXX-XXXX   # TODO: fill in from https://orcid.org
    affiliation: 1
  - name: Shagun Agarwal
    orcid: 0000-XXXX-XXXX-XXXX   # TODO: fill in from https://orcid.org
    affiliation: 1
affiliations:
  - name: Coastal and Marine Offshore Engineering, Delft University of Technology,
          the Netherlands
    index: 1
date: 2026-05-31
bibliography: paper.bib
---

# Summary

Very large floating structures (VLFS) — offshore wind platforms, floating
breakwaters, flexible wave energy converters, and proposed floating airports —
experience wave-induced deformation that can determine both their structural
integrity and their hydrodynamic performance.  Predicting this *hydroelastic*
response requires a coupled model: the fluid pressure drives structural motion,
and structural motion reshapes the fluid boundary conditions, so the two
systems must be solved together.  Traditional approaches treat the fluid and
structure as separate sub-problems, exchanging forces and displacements in a
staggered loop, which introduces splitting errors and limits achievable
accuracy for strongly coupled configurations.

HydroElasticFEM.jl is an open-source Julia package [@Bezanson2017] for finite
element (FEM) simulation of hydroelastic wave–structure interaction.  It
assembles the linearised fluid velocity potential, the free-surface elevation,
and all structural displacement fields into a single multi-field linear system
and solves it in one step — a *monolithic* strategy that avoids splitting
errors by design.  The package is built on Gridap.jl [@Verdugo2022], a
high-level Julia FEM framework, and supports both structured Cartesian
meshes and three-dimensional unstructured meshes imported from Gmsh
[@Geuzaine2009].

The package provides frequency-domain and time-domain simulation modes. Frequency-domain mode solves a single complex-valued linear system for monochromatic wave excitation and is the standard tool for VLFS resonance studies and transfer-function calculations. Time-domain mode integrates the equations with the Generalised-$\alpha$ method, giving access to transient effects and irregular sea states. Both modes share the same physics entities and mesh infrastructure; switching between them requires changing only a configuration object.

# Statement of Need

Hydroelastic analysis of VLFS demands an open, extensible, and
computationally accessible simulation tool.  Existing commercial packages such
as WAMIT and ANSYS Aqua provide robust hydrodynamic solvers but are closed
source, treat structures as rigid bodies or rely on loosely coupled structural
sub-solvers, and do not expose the governing equations as user-extensible code.
High-fidelity CFD-based FSI frameworks such as OpenFOAM carry computational
costs orders of magnitude higher than linear potential-flow models for the
wave-frequency regime where VLFS analysis is most relevant.  Open-source
boundary-element codes address the fluid side of the problem but provide
limited structural flexibility and typically lack 3D or unstructured-mesh
support.

HydroElasticFEM.jl targets ocean engineers, offshore structure designers, and
academic researchers who need: (i) a physically transparent model in which
the governing equations appear as readable Julia code; (ii) full 2D and 3D
support on both structured and unstructured meshes; and (iii) a plugin
architecture that allows new structural models to be added without modifying
the assembly core.  By building on Gridap.jl, the package inherits automatic
differentiation, a broad finite element space library, and an active
open-source community.

# State of the Field

Commercial tools for hydroelastic analysis — WAMIT, ANSYS Aqua, and HydroDyn
(part of OpenFAST) — are widely used but carry licensing restrictions and, in
most cases, implement staggered fluid–structure coupling.  OpenFOAM-based FSI
approaches are fully open but are impractical for parametric frequency-domain
studies.

The direct predecessor of HydroElasticFEM.jl is MonolithicFEMVLFS.jl
[@Colomes2022], which demonstrated monolithic FEM for a 2D Euler–Bernoulli
beam on a Cartesian mesh and validated against the Liu floating-beam benchmark
[@Liu1992] and the submerged plate results of [@Yago2021].  HydroElasticFEM.jl
extends that work in four ways: (i) 3D support via a dimension-generic
`TankDomain{D}` and a `GmshDomain` interface for external Gmsh meshes
[@Geuzaine2009]; (ii) Kirchhoff–Love plate and Timoshenko beam physics under
the same plugin interface as the original Euler–Bernoulli beam; (iii)
point-mass resonators for locally-resonant metamaterial applications; and (iv)
a validated rotational-spring joint formulation reproducing the beam-joint
benchmark of [@Khabakhpasheva2002].

# Software Design

HydroElasticFEM.jl is organized into three strictly layered modules: Geometry,
Physics, and Simulation.

The **Geometry layer** converts a user-specified domain into named Gridap
triangulations and quadrature measures.  `TankDomain{D}` generates Cartesian
meshes for 2D and 3D numerical wave tanks, including embedded structural
sub-domains and sponge-layer damping zones.  `GmshDomain` wraps an external
`.msh` file; boundaries are identified by physical-group names, so the same
physics code runs on any mesh that defines the required groups.

The **Physics layer** defines the equations.  Each physics entity is a
parameter struct that declares its field symbols, finite element configuration,
and which bilinear forms — mass, damping, stiffness, and right-hand side — it
contributes.  Cross-entity coupling (e.g., the kinematic condition at the
wetted surface that couples fluid pressure to structural velocity) is declared
by specializing trait functions such as `has_damping_form(::PotentialFlow,
::EulerBernoulliBeam)` in a dedicated coupling file.  The Simulation layer
detects all active coupling pairs automatically by evaluating these traits over
all entity combinations, so new structural models require no changes to the
assembly core.

The Euler–Bernoulli beam and Kirchhoff–Love plate formally require $C^1$
inter-element continuity, which standard Lagrange elements cannot provide.
HydroElasticFEM.jl uses the Symmetric Interior Penalty Galerkin (SIPG)
discontinuous Galerkin method [@Colomes2022] to enforce slope continuity weakly
across interior facet skeletons, achieving spectral accuracy without
specialised $C^1$ elements.

The **Simulation layer** orchestrates four steps: building the discrete mesh
model, extracting named sub-triangulations, constructing multi-field trial and
test spaces for all active physics entities, and assembling the global FE
operator.  Users interact with two functions: `build_problem`, which accepts a
domain, a list of physics entities, and a configuration; and `simulate`, which
returns a `SimResult` from which individual solution fields are recovered by
name.

# Research Impact Statement

HydroElasticFEM.jl enables research on VLFS deflection under irregular seas,
parametric design of floating breakwaters and wave energy converters with
compliant elastic components, and investigation of resonator-enhanced wave
attenuation in locally resonant metamaterial coatings.  The package is the
principal simulation tool for hydroelastic FSI research at the Coastal and
Marine Offshore Engineering (CMOE) group at TU Delft.  Its predecessor
MonolithicFEMVLFS.jl was used in the study of [@Colomes2022]; the extended
3D and Gmsh capabilities of HydroElasticFEM.jl support ongoing projects
on three-dimensional floating structures and wave energy devices within the
CMOE group.

# AI Usage Disclosure

AI assistance was used during the development of HydroElasticFEM.jl and the
preparation of this paper.  Specifically, Claude (Anthropic) and GitHub
Copilot (Microsoft/OpenAI) were used to assist with code generation,
documentation drafting, and code review.  All AI-generated content was
reviewed, verified, and modified by the authors.  The authors take full
responsibility for the accuracy and correctness of all submitted materials.

# Acknowledgements

The authors thank the Gridap.jl development team for their framework and
for responsive support.  This work is part of the research activities of the
Coastal and Marine Offshore Engineering (CMOE) group at Delft University of
Technology.

# References
