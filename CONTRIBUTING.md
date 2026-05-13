# Contributing to HydroElasticFEM.jl

## 1. Development Environment Setup

```bash
git clone https://github.com/CMOE/HydroElasticFEM.jl
cd HydroElasticFEM.jl
julia --project -e "using Pkg; Pkg.instantiate()"
julia --project -e "using Pkg; Pkg.test()"
```

All 865 tests should pass (2 broken on macOS ARM64 — the Liu Gmsh benchmark —
are expected and are not regressions).

## 2. Code Style

We follow the [Gridap.jl style guide](https://github.com/gridap/Gridap.jl/blob/master/CONTRIBUTING.md).

- 2 spaces for indentation level
- No trailing whitespace
- CamelCase for type names
- Pluralized CamelCase for files that implement a type (`EulerBernoulliBeam.jl`)
- `CamelCaseTests` for test file names (`EulerBernoulliBeamTests.jl`)
- lowercase with underscores for methods (`build_fe_spaces`)
- Use whitespace for readability
- 80-character line length limit
- Use `method!` for mutating methods
- Wrap multiline expressions in parentheses to avoid ambiguity
- `import Module as M` in source files — never bare `using` inside `src/`
- `@with_kw` (Parameters.jl) for all parameter structs with default fields
- Docstrings: use `[unit]; default value` (not `[unit] (default value)`) to avoid
  Documenter.jl parsing the parenthetical as a Markdown link

See [the Julia CONTRIBUTING.md](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md)
for further guidance.

## 3. Adding a New Physics Entity — Checklist

Before opening a PR for a new structural entity, complete every step:

1. Create `src/Physics/Structures/MyEntity.jl` with a `@with_kw` struct subtyping
   `Structure <: PhysicsParameters`.
2. Implement `variable_symbol(s::MyEntity) = s.symbol`.
   Override `variable_symbols` and `field_fe_configs` only for multi-field entities.
3. Implement the active weak form methods (`mass`, `damping`, `stiffness`, `rhs`).
4. Override `has_damping_form(::MyEntity) = false` (and other absent form traits).
5. Add `include("Structures/MyEntity.jl")` inside `src/Physics/Physics.jl`.
6. If the entity couples to `PotentialFlow`, add `has_damping_form(::PotentialFlow, ::MyEntity) = true`
   plus the `damping(pf, s, ...)` method to `src/Physics/CouplingTerms.jl`.
7. Add `export MyEntity` to `src/HydroElasticFEM.jl`.
8. Write a docstring with all fields in SI units using `; default value` style.
9. Add `HydroElasticFEM.Physics.MyEntity` to `docs/src/api/physics.md`.
10. Write tests in `test/Physics/MyEntityTests.jl` and include them in `test/runtests.jl`.
11. Add `include("Physics/Structures/MyEntity.jl")` to `src/Physics/Structures/AGENT.md`
    entity table.

## 4. Test Conventions

| Test type | Description | Tolerance |
|---|---|---|
| Unit test | One entity in isolation, no fluid coupling | exact or 1e-10 |
| Integration test | Full FSI pipeline end-to-end | ≤ 1% for well-resolved meshes |
| Benchmark (coarse) | Published result, fast CI run | ≤ 5% |
| Benchmark (fine) | Published result, production run | ≤ 1% |

**CI budget**: each test file must complete in under 60 seconds.
Use `nx=4, ny=4` (or equivalent) for structural unit tests.
Expensive convergence studies belong in `benchmark/` or `scripts/`, not `test/`.

Every new entity needs at least one quantitative test with a stated tolerance and
a reference (analytical formula or published paper result).

## 5. Branch Naming

| Prefix | Use for |
|---|---|
| `feature/X` | New functionality |
| `fix/X` | Bug fixes |
| `docs/X` | Documentation-only changes |
| `test/X` | New or improved tests |
| `refactor/X` | Internal refactoring without behaviour change |

## 6. Commit Message Convention

Use a short prefix followed by a colon and an imperative sentence:

```
feat: add KirchhoffLovePlate entity with SIPG weak form
fix: correct FESpace type in build_fe_spaces accumulator
docs: architecture guide for new developers
test: add Timoshenko thin-limit convergence test
refactor: extract _assemble_form helper in FEOperators
```

Keep the subject line under 72 characters.
Add a body paragraph if the change is non-obvious.

## 7. Pull Request Checklist

Before requesting review, verify:

- [ ] `julia --project -e "using Pkg; Pkg.test()"` passes (865 passed, 2 broken)
- [ ] `julia --project=docs docs/make.jl` builds with zero errors and zero warnings
- [ ] New entity has a docstring listing all fields with SI units
- [ ] New entity appears in `docs/src/api/physics.md` under `@docs`
- [ ] `src/HydroElasticFEM.jl` exports the new type
- [ ] At least one quantitative test with tolerance is in `test/`

## 8. Legacy Code Policy

- Do not add new code under `archive/`.
- Legacy scripts are removed once equivalent behavior exists in `src/` + `examples/`.
- New features must be implemented in `src/` and covered by tests under `test/`.

## 9. Scripts Policy

The `scripts/` directory is for thin reproducible local runners that call current `src/`
modules. Do not duplicate physics or assembly logic inside scripts.
When a script workflow becomes a supported user workflow, promote it to `examples/` and
add test coverage.