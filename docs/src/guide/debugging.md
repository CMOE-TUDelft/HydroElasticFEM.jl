# Debugging Guide: When Simulations Go Wrong

This guide covers the most common failure modes in HydroElasticFEM.jl simulations —
from Julia error messages you will see in the REPL, to physically wrong but silently
passing results, to convergence problems you only detect by comparing with a reference.

## A. Reading Error Messages

### `MethodError: no matching method for variable_symbol`

```
MethodError: no matching method for `variable_symbol(::MyNewEntity)`
```

You defined a new entity but did not implement `variable_symbol`.
The fix is always the same:

```julia
variable_symbol(s::MyNewEntity) = s.symbol   # or whatever field name you chose
```

If you are defining a **multi-field entity** (like `TimoshenkoBeam`), also override
`variable_symbols`:

```julia
variable_symbols(s::MyNewEntity) = (s.symbol_a, s.symbol_b)
```

### `KeyError` during `build_fe_spaces`

```
KeyError: key :Γstrut not found
```

The `space_domain_symbol` field of your entity points to a key that does not exist in
`TankTriangulations`.
Either you misspelled the key, or the triangulation for that region was not registered.
Check the table in the [Architecture Guide](@ref "Architecture Guide") for standard keys,
and verify that the `StructureDomain` or `DampingZone` you declared in `TankDomain` uses
the matching `domain_symbol`.

### `DimensionMismatch` during assembly

```
DimensionMismatch: dimensions must match: a has dims (240,), b has dims (480,)
```

This usually means two FE spaces that are being added together have incompatible sizes.
Common causes:

1. One entity uses `order = 1` and another uses `order = 2` on the same triangulation,
   with coupling between them — the assembled cross-block has the wrong shape.
2. The `space_domain_symbol` points to the bulk mesh (`:Ω`) when it should point to a
   surface mesh (`:Γη`), so the DOF count is wrong.

Print the DOF counts of each space to diagnose:

```julia
problem = build_problem(tank, physics, config)
Y = get_test_fe_space(problem)
for (i, V) in enumerate(Y)
    println("Field $i: $(num_free_dofs(V)) free DOFs")
end
```

### `LinearAlgebra.SingularException`

```
SingularException(42)
```

The assembled stiffness matrix is singular.
This happens when:

1. **A Dirichlet BC is missing** — e.g., the potential field `ϕ` has no Dirichlet or
   Robin condition on any boundary, leaving a pure Neumann problem with a one-dimensional
   null space (constant mode).
   Fix: ensure `RadiationBC` is attached to both inlet and outlet, or add a point
   constraint.
2. **An entity contributes only to the stiffness form at `ω = 0`** — the
   frequency-domain system uses `-ω² m + (-iω) c + k`.
   At `ω = 0` the mass and damping terms vanish, leaving only `k`, which may be
   rank-deficient for a pure spring without a rigid-body constraint.
3. **Wrong `space_domain_symbol`** — an entity's FE space is built on an empty
   triangulation, producing a zero-DOF space that makes the block row/column all zeros.

### `BoundsError` near resonator assembly

```
BoundsError: attempt to access 0-element Vector{...} at index [1]
```

You passed an empty `Vector{ResonatorSingle}` to the physics array.
The resonator assembler expects at least one element.
Either remove the resonator from the `physics` vector, or initialize it with at least
one `ResonatorSingle`.

## B. Physics Sanity Checks

Run these checks after every new simulation to catch silent errors before publishing.

### 1. Flux conservation

In a lossless empty tank, the horizontal flux through any vertical cross-section must
equal the inlet flux.
Compute the inlet flux and the mid-tank flux and compare:

```julia
Γin   = get_triangulations(problem)[:Γin]
dΓin  = Measure(Γin, 4)
nin   = get_normal_vector(Γin)
flux_in  = sum(∫(∇(ϕh) ⋅ nin)dΓin)

# Create a cross-section at x = L/2 and repeat.
# For a simple check, verify |flux_in| is approximately ω * η0 * H0 * cosh(kH0)/sinh(kH0).
```

### 2. Symmetry check

For a tank that is symmetric about `x = L/2` with a symmetric incident wave, the
free-surface amplitude `|κh(x)|` should equal `|κh(L-x)|` to within numerical precision.
If symmetry is broken, check whether the mesh is symmetric and whether the damping zones
have the same parameters on both sides.

### 3. Low-frequency rigid-body limit

At very low `ω → 0`, a free-floating structure must move with the wave — its
displacement amplitude should approach `η₀`.
Run the simulation at `ω = 0.1` (much lower than the structural resonance) and verify
that `maximum(abs.(ηh(probes))) / η0 ≈ 1`.

### 4. High-frequency cut-off

At very high `ω`, short waves cannot penetrate significantly under a large floating
structure.
The amplitude directly below the structure's midpoint should decay towards zero.
Plotting `|ηh|` versus frequency is the quickest way to see this transition and to
locate the natural frequencies.

## C. Convergence Debugging

### Is the error from mesh resolution or formulation?

Double the number of cells in each direction and rerun.
If the error in a chosen quantity (e.g., `|ηh(x_probe)|`) drops by a factor of
$2^{p+1}$ (where $p$ is the polynomial order), the solution is converging at the
expected rate and the formulation is correct.
If the error does not decrease, or decreases at a rate slower than $O(h)$, the
formulation has a bug.

A practical script:

```julia
errors = Float64[]
for nx in [20, 40, 80]
    tank   = TankDomain(L = 60.0, H = 10.0, nx = nx, ny = nx÷6)
    problem = build_problem(tank, physics, config)
    result  = simulate(problem)
    _, κh  = result.solution
    push!(errors, abs(κh(Point(30.0, 0.0)) - κ_reference))
end
println(errors[1] / errors[2], "  (expect ≈ 2^(p+1) = 4 for p=1)")
println(errors[2] / errors[3])
```

### Polynomial order check

Increase `order` by 1 in `FESpaceConfig` and rerun on the same mesh.
If the error drops by more than a factor of 2, the solution has not yet converged and you
need either more cells or a higher order.
If the error does not change, the solution is dominated by something other than
polynomial approximation error — check boundary conditions and coupling terms.

### C/DG oscillations: penalty too small

If a beam or plate solution oscillates wildly between elements — easily visible in the
VTK output as a checker-board pattern in `eta_re` — the SIPG penalty `γ` is too small.
The recommended value is $\gamma = p(p-1)$ where $p$ is the polynomial order.
Increase `γ` in `FESpaceConfig`:

```julia
fe = FESpaceConfig(order = 4, vector_type = Vector{ComplexF64}, γ = 12.0)
# 4 * (4-1) = 12
```

If increasing `γ` makes the solution worse (introduces locking), the mesh is too coarse
for the chosen order — refine first, then tune `γ`.

## D. Using VTK Output for Visual Debugging

### Enable output

After calling `simulate`, write VTK files with:

```julia
trians = get_triangulations(problem)

writevtk(trians[:Ω],  "fluid",
    cellfields = ["phi_re" => real(ϕh), "phi_im" => imag(ϕh)])
writevtk(trians[:Γκ], "free_surface",
    cellfields = ["kap_re" => real(κh), "kap_im" => imag(κh)])
writevtk(trians[:Γη], "structure",
    cellfields = ["eta_re" => real(ηh), "eta_im" => imag(ηh)])
```

Each call produces a `.vtu` file in the current directory.
For time-domain simulations, Gridap can write a `.pvd` collection automatically.

### What to look for in Paraview

1. **`phi_re` in `Ω`:** Should show a smooth wave pattern propagating in the x-direction.
   Discontinuities or sharp local extrema indicate a mesh resolution problem or a wrong
   boundary condition.

2. **`kap_re` on the free surface:** Should match the incident wave amplitude far from
   the structure (`|κ| ≈ η₀`) and show the expected standing-wave or transmitted-wave
   pattern near the structure.
   A flat `kap_re ≈ 0` everywhere means the free-surface coupling is not active —
   check that `FreeSurface` is in the physics array and that coupling is detected.

3. **`eta_re` on the structure:** Look for the expected mode shape.
   For a simply supported beam at resonance, you expect a half-sine shape.
   A step function or zero field usually means a Dirichlet BC is incorrectly applied, or
   the coupling to the fluid is missing.

4. **Velocity field from `∇ϕ`:** Add the gradient as a vector field to spot reversed
   flow or incorrect inlet/outlet conditions:
   ```julia
   writevtk(trians[:Ω], "velocity",
       cellfields = ["v_re" => real(∇(ϕh)), "v_im" => imag(∇(ϕh))])
   ```
   The `v_re` vector in the bulk should point in the direction of wave propagation and
   be zero at the seabed (no-penetration condition satisfied).
