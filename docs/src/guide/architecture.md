# Architecture Guide

This guide is written for contributors with a finite-element background who are new to
HydroElasticFEM.jl and [Gridap.jl](https://github.com/gridap/Gridap.jl).
It explains how the package is structured, why it is structured that way, and how to
extend it with new physics entities.

## 1. Big Picture

HydroElasticFEM.jl solves hydroelastic wave–body interaction problems using a monolithic
finite-element approach on a single fluid–structure mesh.
The code is organized into three layers that communicate strictly top-to-bottom:

```
┌──────────────────────────────────────────────────────────────┐
│  Geometry                                                    │
│  TankDomain / GmshDomain                                     │
│    → build_model      → DiscreteModel                        │
│    → build_triangulations → TankTriangulations               │
│    → get_integration_domains → IntegrationDomains            │
└─────────────────────────┬────────────────────────────────────┘
                          │  triangulations + measures
┌─────────────────────────▼────────────────────────────────────┐
│  Physics                                                     │
│  PotentialFlow, FreeSurface, EulerBernoulliBeam, …           │
│  Each entity declares: variable_symbol, FESpaceConfig,       │
│    space_domain_symbol, and weak form methods                 │
└─────────────────────────┬────────────────────────────────────┘
                          │  FE spaces + assembled operators
┌─────────────────────────▼────────────────────────────────────┐
│  Simulation                                                  │
│  build_problem → build_fe_spaces → build_fe_operator         │
│  simulate → SimResult                                        │
└──────────────────────────────────────────────────────────────┘
```

The **Geometry** layer converts a user-specified tank configuration or external mesh file
into a set of Gridap triangulations and quadrature measures.
The **Physics** layer defines what equations live on which sub-domain and what FE
discretization they use, entirely independently of the mesh format.
The **Simulation** layer stitches everything together: it assembles a single multi-field FE
operator from all entity contributions and calls Gridap's linear or transient solver.

This separation means that switching from a 2D Cartesian tank to a 3D Gmsh mesh, or
adding a new structural model, requires changes in exactly one layer without touching the
others.

## 2. Layer 1: Geometry

Every domain type implements the `AbstractDomain` interface:

| Function | Returns | Meaning |
|---|---|---|
| `triangulation(d)` | `Triangulation` | Bulk mesh |
| `boundary_tags(d)` | `Dict{Symbol,Any}` | Named boundary tags |
| `ambient_dimension(d)` | `Int` | Embedding space dimension |
| `manifold_dimension(d)` | `Int` | Intrinsic mesh dimension |
| `get_boundary(d, tag)` | `Boundary` | Named boundary triangulation |

### TankDomain{D}

`TankDomain{D}` is the built-in Cartesian domain for 2D (`D = 2`) and 3D (`D = 3`)
numerical wave tanks.
In 2D it supports embedded structural sub-domains (`StructureDomain`), sponge layers
(`DampingZone`), and rotational-spring joint locations (`JointDomain`).
`build_model(domain)` generates a `DiscreteModel` from the Cartesian description, and
`build_triangulations(domain, model)` extracts and names all sub-triangulations into a
`TankTriangulations` struct.

### GmshDomain

`GmshDomain` wraps an external Gmsh `.msh` file.
It requires the following physical group names to be defined in the file:
`"fluid"`, `"seabed"`, `"inlet"`, `"outlet"`, `"free_surface"`, and `"structure"`.
The function `validate_gmsh_tags` checks for these at load time and reports any missing
groups before the simulation starts.
Optional damping-zone physical groups may be added with arbitrary names and referenced
through the `boundary_conditions` field of `PotentialFlow`.

### TankTriangulations and IntegrationDomains

`TankTriangulations` is a Dict-based container holding all named sub-triangulations.
The standard keys are:

| Key | Sub-domain |
|---|---|
| `:Ω` | Fluid bulk |
| `:Γκ` | Free-surface (elevation field) |
| `:Γη` | Structure wetted surface |
| `:Γbot` | Seabed |
| `:Γin`, `:Γout` | Inlet / outlet lateral boundaries |
| `:Λη` | DG skeleton facets on the structural surface |
| `:Γ_dampings` | Damping-zone boundary segments |

`IntegrationDomains` wraps a `TankTriangulations` and converts each triangulation into a
Gridap `Measure` object (`:dΩ`, `:dΓκ`, `:dΓη`, `:dΛη`, etc.) at user-specified
quadrature orders.
Physics entities receive an `IntegrationDomains` during assembly and index it by symbol,
for example `dom[:dΓη]`.

## 3. Layer 2: Physics

### The PhysicsParameters Interface

Every physics entity is a concrete subtype of `PhysicsParameters`.
The required interface consists of form-presence traits and two identification functions:

| Function | Signature | Default |
|---|---|---|
| `variable_symbol` | `(s) → Symbol` | must implement |
| `variable_symbols` | `(s) → NTuple{N,Symbol}` | wraps `variable_symbol(s)` |
| `field_fe_configs` | `(s) → NTuple{N,FESpaceConfig}` | returns `(s.fe,)` |
| `has_mass_form` | `(s) → Bool` | `true` |
| `has_damping_form` | `(s) → Bool` | `true` |
| `has_stiffness_form` | `(s) → Bool` | `true` |
| `has_rhs_form` | `(s) → Bool` | `true` |

Every entity also carries a `space_domain_symbol::Symbol` field that tells `build_fe_spaces`
which key in `TankTriangulations` to use when building the FE space for that entity.

### Built-in Entities

| Entity | `variable_symbols` | `space_domain_symbol` |
|---|---|---|
| `PotentialFlow` | `(:ϕ,)` | `:Ω` |
| `FreeSurface` | `(:κ,)` | `:Γκ` |
| `EulerBernoulliBeam` | `(:η_b,)` | `:Γη` |
| `TimoshenkoBeam` | `(:w, :θ)` | `:Γη` |
| `KirchhoffLovePlate` | `(:η,)` | `:Γη` |
| `Membrane` | `(:η_m,)` | `:Γη` |
| `ResonatorSingle` | `(:q,)` | `:Ω` |

### Coupling Detection

Cross-entity coupling is detected automatically by `detect_couplings`.
The function iterates over all ordered pairs `(a, b)` with `a ≠ b` and checks whether
any of `has_mass_form(a, b)`, `has_damping_form(a, b)`, `has_stiffness_form(a, b)`, or
`has_rhs_form(a, b)` returns `true`.
Cross-entity weak form methods (e.g. the fluid–structure kinematic condition on `Γη`) are
defined in `src/Physics/CouplingTerms.jl` by specializing those four trait functions for
the relevant type pair.
All coupling traits default to `false` for unknown pairs, so adding a new entity that
does not couple to any existing one requires no changes to `CouplingTerms.jl`.

### C/DG Formulation for High-Order Beam and Plate Entities

The Euler–Bernoulli beam and Kirchhoff–Love plate require $C^1$ inter-element continuity,
which is unavailable in standard Gridap element libraries.
HydroElasticFEM.jl uses the Symmetric Interior Penalty Galerkin (SIPG) method to enforce
slope continuity weakly across the DG skeleton $\Lambda_\eta$:

```math
\frac{\gamma \, EI_\rho}{h} \int_{\Lambda_\eta}
    \llbracket \nabla w \rrbracket \cdot \llbracket \nabla v \rrbracket \, \mathrm{d}\Lambda.
```

The dimensionless penalty coefficient $\gamma$ is stored in `FESpaceConfig.γ` and should
be set to $p(p-1)$ for polynomial order $p$.
The skeleton measure `:dΛη` is provided by `IntegrationDomains` and accessed inside the
`stiffness` method of each C/DG entity.

## 4. Layer 3: Simulation

### The Four-Step Pipeline

`build_problem(domain, physics, config)` orchestrates the following steps:

```
build_problem(domain, physics, config)
 │
 ├─ 1.  build_model(domain)               → DiscreteModel
 │
 ├─ 2.  build_triangulations(domain, model)  → TankTriangulations
 │
 ├─ 3.  get_integration_domains(trians)   → IntegrationDomains
 │
 └─ 4.  build_fe_spaces(physics, trians, config)
 │         for each entity e:
 │           trian = trians[e.space_domain_symbol]
 │           Vᵢ    = TestFESpace(trian, ReferenceFE(…); …)
 │           Uᵢ    = TrialFESpace(Vᵢ, …)
 │         X    = MultiFieldFESpace([U₁, U₂, …])
 │         Y    = MultiFieldFESpace([V₁, V₂, …])
 │         fmap = Dict(:ϕ => 1, :κ => 2, :η_b => 3, …)
 │
 └─ 5.  build_fe_operator(…)              → FEOperator / TransientFEOperator
           detect_couplings(entities) → coupling pairs
           assemble bilinear + linear forms using fmap
           return Gridap operator
```

`simulate(problem)` wraps Gridap's `solve()` and returns a `SimResult`.
Fields are extracted in the order set by `fmap`:

```julia
problem = build_problem(tank, [potential, free_surface, beam], config)
result  = simulate(problem)
ϕh, κh, ηh = result.solution
```

### FieldMap

Inside every weak form method, trial and test functions arrive as a `FieldMap` — a struct
wrapping a positional Gridap tuple and a `Dict{Symbol, Int}` index map.
Symbol-based indexing (`x[:ϕ]`, `y[:η_b]`) decouples each entity's weak form from its
global position in the multi-field system, so inserting a new entity never requires
rewriting existing forms.
`assemble_weakform` loops over entities and coupling pairs, calling
`P.weakform(entity, ctx, x, y)` for each and summing contributions.

## 5. Frequency-Domain Monolithic System

After the substitution $\partial_t \to -i\omega$, the monolithic bilinear form for
each entity $s$ is composed from its linear forms:

```math
a_s(u, v) = -\omega^2 \, m_s(u, v) + (-i\omega) \, c_s(u, v) + k_s(u, v).
```

The full system sums contributions from all entities and all detected coupling pairs:

```math
\sum_{s} a_s(u_h, v) + \sum_{(a,b)} a_{ab}(u_h, v) = \sum_{s} l_s(v)
    \qquad \forall\, v \in Y.
```

The `FrequencyAssemblyContext` carries $\omega$ and the `IntegrationDomains`.
Entities access both through `AC.frequency(ctx)` and `AC.domains(ctx)` respectively.
The same `mass`, `damping`, and `stiffness` methods are reused for the time-domain path;
only the assembly context type changes, and Gridap's `GeneralizedAlpha2` integrator calls
the residual and Jacobian forms directly.

## 6. Recipe: Adding a New Structural Entity

The following checklist adds a single-field structural entity called `MyBeam`.

1. Create `src/Physics/Structures/MyBeam.jl`.

2. Define a `@with_kw` parameter struct that subtypes `Structure <: PhysicsParameters`:
   ```julia
   @with_kw struct MyBeam <: Structure
       EIᵨ::Float64
       mᵨ::Float64
       fe::FESpaceConfig
       space_domain_symbol::Symbol = :Γη
       symbol::Symbol              = :η_mb
   end
   ```

3. Implement `variable_symbol`:
   ```julia
   variable_symbol(s::MyBeam) = s.symbol
   ```

4. Implement the relevant weak form methods, accessing fields via `FieldMap`:
   ```julia
   function mass(b::MyBeam, dom::IntegrationDomains, x_tt, y)
       η_tt = x_tt[b.symbol];  v = y[b.symbol]
       ∫(b.mᵨ * v * η_tt)dom[:dΓη]
   end

   function stiffness(b::MyBeam, dom::IntegrationDomains, x, y)
       η = x[b.symbol];  v = y[b.symbol]
       γ = b.fe.γ;  h = b.fe.order   # penalty and mesh size from FESpaceConfig
       ∫(b.EIᵨ * Δ(η) ⊙ Δ(v))dom[:dΓη] +
           ∫(γ * b.EIᵨ * (∇∇(v) ⊙ mean(∇∇(η))))dom[:dΛη] # SIPG penalty
   end
   ```

5. Override `has_damping_form(::MyBeam) = false` (and any other absent form traits)
   so that `detect_couplings` and the assembler skip those contributions.

6. Add `include("Structures/MyBeam.jl")` in `src/Physics/Physics.jl`, inside the
   structural entity block.

7. If `MyBeam` couples to `PotentialFlow` (kinematic condition on `Γη`), add
   specializations to `src/Physics/CouplingTerms.jl`:
   ```julia
   has_damping_form(::PotentialFlow, ::MyBeam) = true

   function damping(pf::PotentialFlow, b::MyBeam, dom::IntegrationDomains, x_t, y)
       ϕₜ = x_t[variable_symbol(pf)];  ηₜ = x_t[b.symbol]
       w  = y[variable_symbol(pf)];    v  = y[b.symbol]
       ∫(v * ϕₜ - w * ηₜ)dom[:dΓη]
   end
   ```

8. Export `MyBeam` from `src/HydroElasticFEM.jl`.

9. Write a docstring for `MyBeam` that lists every field with its SI unit and a
   minimal example (required by `checkdocs = :exports` in `docs/make.jl`).

10. Add `HydroElasticFEM.Physics.MyBeam` to the `@docs` block in
    `docs/src/api/physics.md`.

11. Add tests to `test/Physics/` that verify: (a) the stiffness matrix is symmetric,
    (b) the weak form reproduces a known analytical solution or converges at the
    expected polynomial rate under mesh refinement.
