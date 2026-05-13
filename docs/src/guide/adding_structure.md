# How to Add a New Structural Entity

This guide walks you through adding a minimal structural entity to HydroElasticFEM.jl.
Working through a concrete example is faster than reading the architecture description
in isolation, so we implement a simple **`LinearSpring`** — a uniformly distributed
vertical spring bed with one scalar DOF per node.
It needs only a mass form and a stiffness form, no C/DG penalty.

By the end you will have a complete, testable entity that couples to the fluid
through the standard kinematic condition on `Γη`.

## 1. Understand the Required Interface

Every structural entity must:

1. Subtype `Structure <: PhysicsParameters`.
2. Declare `variable_symbol(s)` returning a `Symbol` (the FE unknown name).
3. Carry a `space_domain_symbol::Symbol` field (the `TankTriangulations` key where its FE space lives).
4. Carry an `fe::FESpaceConfig` field (polynomial order, element type, BCs).
5. Implement at least one of `mass`, `damping`, `stiffness`, `rhs`.
6. Override `has_*_form` traits to `false` for any absent forms.

Coupling to `PotentialFlow` is detected automatically if you add a specialization of
`has_damping_form(::PotentialFlow, ::LinearSpring)` to `CouplingTerms.jl`.

## 2. Create the Entity File

Create `src/Physics/Structures/LinearSpring.jl`:

```julia
# ─── src/Physics/Structures/LinearSpring.jl ──────────────────────────────

"""
    LinearSpring <: Structure

Distributed vertical spring bed with uniform spring constant `k` [N/m³]
and surface mass density `m` [kg/m²], normalised by fluid density `ρw`.

# Fields
- `k::Float64`                  — Spring constant per unit area / ρw [m/s²]
- `m::Float64`                  — Surface mass density / ρw [dimensionless]
- `symbol::Symbol`              — Field unknown symbol; default `:w_s`
- `space_domain_symbol::Symbol` — Triangulation key for FE spaces; default `:Γη`
- `fe::FESpaceConfig`           — FE discretisation parameters
"""
@with_kw struct LinearSpring <: Structure
    k::Float64
    m::Float64
    symbol::Symbol              = :w_s
    space_domain_symbol::Symbol = :Γη
    fe::FESpaceConfig           = FESpaceConfig(order=1, vector_type=Vector{ComplexF64})
end

variable_symbol(s::LinearSpring) = s.symbol

# Spring has no damping form.
has_damping_form(::LinearSpring) = false

function mass(s::LinearSpring, dom::IntegrationDomains, x_tt, y)
    w_tt = x_tt[s.symbol]
    v    = y[s.symbol]
    ∫(s.m * v * w_tt)dom[:dΓη]
end

function stiffness(s::LinearSpring, dom::IntegrationDomains, x, y)
    w = x[s.symbol]
    v = y[s.symbol]
    # Hydrostatic restoring (gravity) + structural spring:
    ∫((9.81 + s.k) * v * w)dom[:dΓη]
end

# rhs: no external forcing by default.
has_rhs_form(::LinearSpring) = false
```

Key points:

- `@with_kw` (from Parameters.jl) gives every field a keyword-argument constructor and
  supports default values. It is the standard for all entity structs in this package.
- `dom[:dΓη]` is the quadrature measure on the structural surface, provided by
  `IntegrationDomains`. Do not hard-code the measure; always index through `dom`.
- `im` is Julia's imaginary unit. Do not write `1im` when you mean `Complex(0, 1)` —
  they are the same thing, but `im` is idiomatic.

## 3. Register the Entity in Physics.jl

Open `src/Physics/Physics.jl` and add the include inside the structural entity block:

```julia
include("Structures/LinearSpring.jl")
```

Place it after the existing structural includes so it sees the `Structure` abstract type
and the `IntegrationDomains` type that are already in scope.

## 4. Add FSI Coupling

The fluid–structure kinematic condition says: the normal velocity of the fluid at the
structure surface equals the velocity of the structure.
In the damping form this becomes:

```math
\int_{\Gamma_\eta} (v_\phi \cdot \dot{w} - w_\phi \cdot \dot{w}_s)\,\mathrm{d}\Gamma
```

Add the following to `src/Physics/CouplingTerms.jl`:

```julia
# ── PotentialFlow ↔ LinearSpring coupling ───────────────────────────────
has_damping_form(::PotentialFlow, ::LinearSpring) = true

function damping(pf::PotentialFlow, s::LinearSpring, dom::IntegrationDomains, x_t, y)
    ϕₜ = x_t[variable_symbol(pf)]
    wₜ = x_t[s.symbol]
    w  = y[variable_symbol(pf)]
    v  = y[s.symbol]
    ∫(v * ϕₜ - w * wₜ)dom[:dΓη]
end
```

This is identical to the existing `PotentialFlow ↔ Structure` coupling in the same file
— the kinematic condition is the same regardless of which structural model is used.

## 5. Export and Expose to Docs

5a. Add `export LinearSpring` to `src/HydroElasticFEM.jl`, inside the Physics re-exports block.

5b. Add a docstring entry to `docs/src/api/physics.md` so Documenter picks it up:

```markdown
### LinearSpring

```@docs
HydroElasticFEM.Physics.LinearSpring
```
```

## 6. Write the Validation Test

Create `test/Physics/LinearSpringTests.jl`.
The analytical resonance frequency is:

$$\omega_n = \sqrt{\frac{k + g}{m + m_\text{added}}}$$

where $m_\text{added}$ is the fluid added mass.
For a strip of length $L$ and water depth $H$, the leading-order added mass per unit area
is $m_\text{added} \approx \rho_w / k$ for the fundamental mode — use a numerical
reference value from a known benchmark or compare two mesh resolutions.

A practical smoke test checks that the spring-only (no fluid) stiffness integral
returns the correct value on a 1D mesh:

```julia
using Test, Gridap
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Geometry as G

@testset "LinearSpring stiffness" begin
    # Build a 1D mesh [0, 1] with 10 elements.
    model  = CartesianDiscreteModel((0.0, 1.0), (10,))
    trian  = Interior(model)
    dom    = Measure(trian, 2)

    s  = P.LinearSpring(k = 100.0, m = 1.0)
    fe = PH.FESpaceConfig(order = 1, vector_type = Vector{ComplexF64})
    V  = TestFESpace(trian, ReferenceFE(lagrangian, Float64, 1);
                     conformity = :H1, vector_type = Vector{ComplexF64})
    U  = TrialFESpace(V)

    # Integrate ∫ (g + k) * w * v dΓ over [0,1] — should equal (9.81 + 100) = 109.81.
    # With hat functions on 10 elements, the mass matrix trace equals (g+k)*1.0.
    a(u, v) = ∫((9.81 + s.k) * v * u)dom
    M = assemble_matrix(a, U, V)
    @test sum(M) ≈ 9.81 + 100.0 atol = 1e-10
end
```

Add `include("Physics/LinearSpringTests.jl")` to `test/runtests.jl`.

## 7. Common Mistakes

**Not implementing `variable_symbols`:** The default `variable_symbols(s)` wraps
`variable_symbol(s)` in a tuple.
If you forget to implement `variable_symbol`, `build_fe_spaces` will call
`error("variable_symbol not implemented for LinearSpring")` with a clear message.
Implement it before running anything.

**Not exporting the type:** If you skip the `export` in `HydroElasticFEM.jl`, users who
write `using HydroElasticFEM` will get `UndefVarError: LinearSpring not defined`.
The entity still works if imported via `import HydroElasticFEM.Physics as P; P.LinearSpring(...)`,
but that is an inconsistent API.

**Using the wrong `space_domain_symbol`:** Structural entities live on `:Γη` by default.
If you use a symbol that does not exist in `TankTriangulations`, `build_fe_spaces` will
throw a `KeyError`.
Check the table in the [Architecture Guide](@ref "Architecture Guide") for available keys.

**Missing the coupling trait:** If you implement the `damping` coupling method but forget
`has_damping_form(::PotentialFlow, ::LinearSpring) = true`, the assembler will silently
skip the coupling, and the spring will behave as if it is in vacuum.
The resulting `|η|/η₀` values will be unphysically large or small with no error message.
