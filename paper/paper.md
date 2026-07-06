---
title: 'HydroElasticFEM.jl: A Julia package for the finite element analysis of hydroelastic wave-structure interaction problems'
tags:
  - Julia
  - finite element method
  - hydroelasticity
  - wave-structure interaction
  - very large floating structures
  - potential flow
  - fully differentiable
authors:
  - name: Oriol Colomés
    orcid: 0000-0002-5552-9695   
    affiliation: 1
  - name: Shagun Agarwal
    orcid: 0000-0003-1922-4242   
    affiliation: 1
affiliations:
  - name: Civil Engineering and Geosciences Faculty, Delft University of Technology, the Netherlands
    index: 1
date: 2026-06-19
bibliography: paper.bib
---

# Summary

Very large floating structures (VLFS) are gaining attention in many applications, for example in offshore wind platforms, floating breakwaters, flexible wave energy converters, and futuristic floating cities. The interaction of VLFS with free surface waves induces structural deformation and structural stresses that are critical for their performance. Predicting the *hydroelastic* response of VLFS in waves requires a coupled model: the fluid pressure drives structural motion, and structural motion affects the wave propagation through fluid boundary conditions. Therefore, the two systems must be solved together.  Traditional approaches treat the fluid and structure as separate sub-problems, exchanging forces and displacements in a
staggered loop, which might introduce splitting errors and limit achievable accuracy for strongly coupled configurations. In additioned partitioned schemes suffer to converge when the added mass effect is large, which is particularly relevant for floating flexible thin structures.

HydroElasticFEM.jl is an open-source Julia package [@bezanson2017julia] for finite element (FE) simulation of hydroelastic wave–structure interaction problems.  It solves a coupled problem for the linearised fluid velocity potential, the free-surface elevation, and, possibly, structural deflections of multi-body flexible structures. HydroelasticFEM.jl uses a monolithic approach to solve the coupled problem, as originally proposed in [@colomes2023monolithic]. The package is built on Gridap.jl [@verdugo2022software], a high-level Julia FE framework, and supports both structured Cartesian meshes and three-dimensional unstructured meshes imported from Gmsh [@geuzaine2009gmsh].

The package provides frequency-domain and time-domain simulation modes. Frequency-domain mode solves a single complex-valued linear system for monochromatic wave excitation and is the standard tool for VLFS eigenvalue analysis and transfer-function calculations. Time-domain mode integrates the equations using the Generalised-$\alpha$ method, giving access to transient effects and irregular sea states. Both modes share the same physics entities and mesh infrastructure. The architecture of HydroElasticFEM.jl enables the user to switch between frequency-domain and time-domain modes by changing only a configuration object.

# Statement of Need

Hydroelastic analysis of VLFS demands an open, extensible, scalable, and computationally accessible simulation tool. Existing commercial packages provide robust hydrodynamic solvers, but mostly using a boundary element method (BEM) approach, and often are closed source. By nature of the solution methodology, these frameworks are limitted in terms of scalability and level of fidelity. Furthermore, available hydrodynamic analysis software treat structures as rigid bodies or rely on loosely coupled structural sub-solvers, and do not expose the governing equations as user-extensible code. On the other hand, high-fidelity tools that solve the Navier-Stokes equations, such as OpenFOAM or DualSPHysics, carry computational costs orders of magnitude higher than linear potential-flow models for the wave-frequency regime where VLFS analysis is most relevant.  

HydroElasticFEM.jl targets ocean engineers, offshore structure designers, and academic researchers who need: (i) a physically transparent model in which the governing equations appear as readable Julia code; (ii) full 2D and 3D support on both structured and unstructured meshes; and (iii) a plugin architecture that allows new structural models to be added without modifying the assembly core. By building on Gridap.jl, the package inherits automatic differentiation, a broad finite element space library, and an active open-source community.

# State of the Field

Commercial tools for potential flow-based hydrodynamic analysis, e.g. WAMIT, ANSYS Aqua, ORCAWAVE, and HydroDyn (part of OpenFAST), are widely used but carry licensing restrictions and, in most cases, implement staggered fluid–structure coupling.  Fully viscous models such as OpenFOAM or DualSPHysics for wave-structure interaction problems are fully open but are impractical for parametric frequency-domain studies.

In the past years the authors of this package developed a formulation for the hydroelastic analysis of flexible floating structures based on a monolithic FE approach [@colomes2023monolithic;@agarwal2024dynamic]. The numerical tests were implemented in the open-souce library MonolithicFEMVLFS.jl [Colomes_MonolithicFEMVLFS_2022], which demonstrated the new monolithic FE formulation for a wide variety if benchmarks.  HydroElasticFEM.jl extends that work in four ways: (i) 2D/3D support via a dimension-generic `TankDomain{D}` and a `GmshDomain` interface for external Gmsh meshes [@geuzaine2009gmsh]; (ii) Membrane, Euler-Bernoulli, Kirchhoff–Love plate and Timoshenko beam physics under the same plugin interface; (iii) point-mass resonators for locally-resonant metamaterial applications; and (iv) arbitrary number of structures with arbitrary number of rotational-spring joint formulation. In addition, HydroElasticFEM.jl can take advantage of differentiable programming intrinsic capabilities to compute gradients using automatic differentiation, as exploited in adjoint-based optimization problems such as in [el2026adjoint].

# Software Design

HydroElasticFEM.jl is organized into three layered modules: Geometry, Physics, and Simulation.

Figure 1 summarizes the software architecture and highlights the four main
design contributions of HydroElasticFEM.jl.

![Figure 1: Software architecture of HydroElasticFEM.jl, showing the Geometry,
Physics, and Simulation layers.](figure1_software_design.svg)

The **Geometry layer** converts a user-specified domain into named Gridap triangulations and quadrature measures.  `TankDomain{D}` generates Cartesian meshes for 2D and 3D numerical wave tanks, including embedded structural sub-domains and sponge-layer damping zones.  `GmshDomain` wraps an external `.msh` file where boundaries are identified by physical-group names, so the same
physics code runs on any mesh that defines the required groups.

The **Physics layer** defines the weak forms of the system.  Each physics entity is a parameter struct that declares its field symbols, finite element configuration, and which bilinear forms it contributes, e.g. mass, damping, stiffness, and right-hand.  Cross-physics coupling (e.g., the kinematic condition at the wetted surface that couples fluid pressure to structural velocity) is declared by specializing trait functions such as `has_damping_form(::PotentialFlow,::EulerBernoulliBeam)` in a dedicated coupling file.  Later, the Simulation layer detects all active coupling pairs automatically by evaluating these traits over
all entity combinations, so new structural models require no changes to the assembly core.

It is important to highlight that the Euler–Bernoulli beam and Kirchhoff–Love plate formulations formally require $C^1$
inter-element continuity, which standard Lagrange elements cannot provide. HydroElasticFEM.jl uses the Symmetric Interior Penalty Galerkin (SIPG) discontinuous Galerkin method [@colomes2023monolithic] to enforce slope continuity weakly across interior facet skeletons, achieving optimal accuracy without specialised $C^1$ elements.

The **Simulation layer** orchestrates four steps: building the discrete mesh model, extracting named sub-triangulations, constructing multi-field trial and test spaces for all active physics entities, and assembling the global FE
operator.  Users interact with two functions: `build_problem`, which accepts a domain, a list of physics entities, and a configuration; and `simulate`, which returns a `SimResult` from which individual solution fields are recovered by
name.

# Research Impact Statement

HydroElasticFEM.jl enables research on VLFS response under irregular seas, parametric design of floating breakwaters and wave energy converters with compliant elastic components, and investigation of resonator-enhanced wave attenuation in locally resonant metamaterial super-structures.  The package is the principal simulation tool for hydroelastic FSI research at the Computational Multiphysics (CMOE) group at TU Delft. HydroElasticFEM.jl supports ongoing projects on three-dimensional floating structures and wave energy devices within the CMOE group and is the basis of several research collaborations on the topic with international partners.

# AI Usage Disclosure

AI assistance was used during the development of HydroElasticFEM.jl and the preparation of this paper.  Specifically, Claude (Anthropic) and GitHub Copilot (Microsoft/OpenAI) were used to assist with code generation, documentation drafting, and code review.  All AI-generated content was reviewed, verified, and modified by the authors.  The authors take full responsibility for the accuracy and correctness of all submitted materials.

# Acknowledgements

The authors thank the Gridap.jl development team for their framework and for responsive support.  This work is part of the research activities of the Computational Multiphysics in Offshore Engineering (CMOE) group at Delft University of
Technology. This software is part of the project DigiOcean4Solar with file number 21225 of the research programme NWO Talent Programme Vidi AES 2023, which is financed by the Dutch Research Council (NWO) under the grant ID https://doi.org/10.61686/OPCTU16570. O. Colomés gratefully acknowledges this NWO funding.

# References
