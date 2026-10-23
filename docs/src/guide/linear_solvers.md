# Linear Solvers and Large Problems

Every HydroElasticFEM simulation ends in a sparse linear solve. In frequency
domain there is one complex system per frequency. In time domain there is one
real system per time step. This page explains which solver is used, when its
factorization is reused, and how to plug in a different solver for large 3D
problems.

## The default: sparse LU

Both `FreqDomainConfig` and `TimeDomainConfig` have a `solver` field. If it
is left as `nothing`, `simulate` uses Gridap's `LUSolver()`, a direct sparse LU
factorization (UMFPACK through Julia's `SparseArrays`).

```julia
config  = FreqDomainConfig(ω = 1.2)                  # LUSolver()
config  = TimeDomainConfig(tf = 100.0)               # LUSolver()
config  = TimeDomainConfig(tf = 100.0, solver = LUSolver())  # explicit
```

## Time domain: the matrix is factorized once

The time-domain operator built by `build_problem` is a Gridap
`TransientLinearFEOperator` with constant stiffness, damping and mass forms.
For such an operator, Gridap's `GeneralizedAlpha2` and `Newmark` solvers build
the system matrix

```math
J = w_0 K + w_1 C + w_2 M
```

once. With a fixed time step the weights ``w_k`` never change, so Gridap
reuses the factorization (`numerical_setup`) for every step. Each step then
costs:

1. assembling the right-hand side (forcing terms at the current time),
2. forming the residual of the current state, and
3. one forward/backward substitution with the stored factors.

The factorization is paid once, at the first step. This is why the first step
of a time-domain run is much slower than the others, and why a time step that
changes during the run would be expensive.

## When LU becomes the bottleneck

In 2D, and for moderate 3D meshes, sparse LU is the fastest option.
For large 3D meshes (roughly beyond a few ×10⁵ unknowns), fill-in makes the LU
factors grow faster than the matrix. Memory then becomes the limit before time
does. Two options, neither of which needs a change in HydroElasticFEM:

### A different direct solver

Any Gridap `LinearSolver` can be passed through `solver`. Gridap's companion
packages provide parallel and out-of-core direct solvers, for example:

- **MUMPS** or **PETSc** solvers through
  [GridapPETSc.jl](https://github.com/gridap/GridapPETSc.jl)
  (`PETScLinearSolver`, configured via PETSc options such as
  `-pc_type lu -pc_factor_mat_solver_type mumps`);
- **Intel MKL Pardiso** through
  [GridapPardiso.jl](https://github.com/gridap/GridapPardiso.jl)
  (`PardisoSolver`).

```julia
using GridapPETSc
options = "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps"
GridapPETSc.with(args = split(options)) do
    config = TimeDomainConfig(tf = 100.0, solver = PETScLinearSolver())
    prob   = build_problem(tank, physics, config; tconfig = tconfig)
    result = simulate(prob, tconfig)
    # iterate result.solution inside the `with` block
end
```

A direct solver keeps the one-factorization-per-run behaviour described above.

### An iterative solver

An iterative Krylov solver (e.g. GMRES through PETSc) avoids the fill-in. It
needs a good preconditioner, though: the monolithic fluid–free-surface–structure
system is indefinite and poorly scaled. The solve also has to be repeated at
every step instead of reusing a factorization. Prefer a direct solver unless
memory rules it out.

## Checklist for large 3D runs

- Keep `Δt` constant, so the factorization is reused.
- Run a short smoke test first: the first step shows the factorization cost
  and memory, and the next steps show the per-step cost.
- Use the lowest polynomial order that resolves the waves. Q2 in 3D has about
  8× the unknowns of Q1 on the same mesh, and much denser factors.
- Grade the mesh in depth (`map` in `TankDomain`) rather than refining it
  uniformly: most of the resolution is needed near the free surface.
- If LU runs out of memory, switch `solver` to MUMPS or Pardiso before
  changing the formulation.
