We follow the following styling guide from [Gridap.jl CONTRIBUTING.md](https://github.com/gridap/Gridap.jl/blob/master/CONTRIBUTING.md)

Gridap Style Guides
===

* 2 spaces for indentation level
* No trailing white spaces
* CamelCase for typenames
* Pluralized CamelCase for files that implement a type
* CamelCasesTests for CamelCase type test file
* Use lowercase for methods, with underscores only when necessary
* Use whitespace for readability
* 80 characterl line length limit
* Use method! for muting methods
* Wrap multiline expressions in parentheses to avoid errors

See [the Julia CONTRIBUTING.md](https://github.com/JuliaLang/julia/blob/master/CONTRIBUTING.md) for further information.

---

Legacy Code Policy
===

We are migrating all maintained workflows to the structured API under `src/` and
user-facing examples under `examples/`.

* Do not add new code under `archive/`.
* Legacy scripts are removed once equivalent behavior exists in `src/` +
	`examples/`.
* New features must be implemented in `src/` and covered by tests under `test/`.

Scripts Policy
===

The `scripts/` directory is retained only for reproducible local runners that
call current `src/` modules.

* Keep scripts thin: parameter setup + invocation of maintained `src/` code.
* Do not duplicate physics/assembly logic inside scripts.
* When a script workflow becomes a supported user workflow, promote it to
	`examples/` and add test coverage.