"""
    SimResult

Container for simulation outputs returned by `simulate`.

# Fields
- `X` — trial multi-field FE space (Gridap `MultiFieldFESpace`)
- `Y` — test multi-field FE space
- `fmap::Dict{Symbol,Int}` — map from field symbol (e.g. `:phi`) to positional
  index in the multi-field solution
- `op` — assembled FE operator (`AffineFEOperator` for frequency domain,
  `TransientLinearFEOperator` for time domain)
- `solution` — solved FE function (frequency domain) or ODE solution iterator
  (time domain)

# Example
```julia
result = S.simulate(problem)
# Frequency domain: unpack fields by position
phi_h, kappa_h = result.solution
# Time domain: iterate over time steps
for (t, u_h) in result.solution
    phi_h, kappa_h = u_h
end
```
"""
struct SimResult
    X
    Y
    fmap::Dict{Symbol, Int}
    op
    solution
end

"""
    body_motion(result::SimResult, body::AbstractRigidBody) -> NamedTuple

Complex amplitudes (frequency domain) of the active modes of `body`, keyed by
mode name, e.g. `body_motion(res, body).heave`.
"""
function body_motion(result::SimResult, body::P.AbstractRigidBody)
    ξh = result.solution[result.fmap[P.variable_symbol(body)]]
    NamedTuple{body.modes}(Tuple(get_free_dof_values(ξh)))
end
