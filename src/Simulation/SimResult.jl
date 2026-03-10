"""
    SimResult

Container for simulation outputs.

# Fields
- `X` — trial multi-field FE space
- `Y` — test multi-field FE space
- `fmap::Dict{Symbol,Int}` — field-symbol to positional-index map
- `op` — FE operator (`AffineFEOperator` or `TransientLinearFEOperator`)
- `solution` — solved FE function (frequency) or ODE solution (time)
"""
struct SimResult
    X
    Y
    fmap::Dict{Symbol, Int}
    op
    solution
end
