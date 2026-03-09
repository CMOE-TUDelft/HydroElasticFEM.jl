# ─────────────────────────────────────────────────────────────
# FESpaceConfig: numerical FE discretisation parameters
# stored inside each physics entity.
# ─────────────────────────────────────────────────────────────

"""
    FESpaceConfig

Numerical FE discretisation parameters stored inside each physics entity.
Controls how `build_fe_spaces` constructs the entity's FE space.

# Fields
- `reffe_type`          — ReferenceFE family, e.g. `lagrangian` (default `lagrangian`)
- `space_type::DataType` — field type for the ReferenceFE (default `Float64`)
- `order::Int`          — polynomial order of the reference FE (default 1)
- `conformity::Symbol`  — FE conformity, e.g. `:H1`, `:L2` (default `:H1`)
- `vector_type::DataType` — Gridap vector type (default `Vector{ComplexF64}`)
- `γ::Float64`          — Nitsche penalty parameter (default `10.0 * order^2`);
                           only used by entities with DG / skeleton terms
- `dirichlet_tags`      — `nothing`, `String`, or `Vector{String}` (default `nothing`)
- `dirichlet_value`     — Dirichlet BC value: `nothing`, a function, or a constant
                           (default `nothing`)
"""
@with_kw struct FESpaceConfig
    reffe_type             = lagrangian
    space_type::DataType  = Float64
    order::Int             = 1
    conformity::Symbol     = :H1
    vector_type::DataType  = Vector{ComplexF64}
    γ::Float64             = 10.0 * order^2
    dirichlet_tags         = nothing
    dirichlet_value        = nothing
end
