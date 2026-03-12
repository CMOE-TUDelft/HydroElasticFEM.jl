"""
    module ParameterHandler

Configuration structs for HydroElasticFEM.

Provides lightweight, `@with_kw`-constructed parameter containers
consumed by other modules:

- **`FESpaceConfig`** — numerical FE discretisation parameters (order, conformity,
  vector type, Dirichlet tags/values) stored inside each physics entity.
- **`SimConfig`** — simulation run settings (domain type, frequency, solver).
- **`TimeConfig`** — time-domain integration parameters (time step, final time,
  initial conditions, spectral radius).

Loaded early in the module dependency chain so that both `PhysicsCore`
(entities) and `Simulation` can depend on these types.
"""
module ParameterHandler

using Parameters
using Gridap

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

"""
    TimeConfig

Time-domain integration parameters.

# Fields
- `Δt::Float64` — time step
- `t₀::Float64` — start time (default 0.0)
- `tf::Float64` — final time
- `ρ∞::Float64` — spectral radius for Generalized-α (default 1.0)
- `u0` — initial condition(s); tuple/vector of interpolatable objects per field
- `u0t` — initial velocity (optional)
- `u0tt` — initial acceleration (optional)
"""
@with_kw struct TimeConfig
    Δt::Float64
    t₀::Float64 = 0.0
    tf::Float64
    ρ∞::Float64 = 1.0
    u0 = nothing
    u0t = nothing
    u0tt = nothing

    @assert Δt > 0 "Δt must be positive"
    @assert tf > t₀ "tf must be greater than t₀"
end

"""
    SimConfig

Configuration for a simulation run.

# Fields
- `domain::Symbol` — `:frequency` or `:time`
- `ω` — angular frequency (required when `domain == :frequency`)
- `degree::Int` — quadrature degree (default 4)
- `solver` — optional solver override (e.g. `LUSolver()`)
"""
@with_kw struct SimConfig
    domain::Symbol
    ω::Union{Float64, Nothing} = nothing
    degree::Int = 4
    solver = nothing

    @assert domain in (:frequency, :time) "domain must be :frequency or :time"
    @assert domain == :time || ω !== nothing "ω is required for frequency-domain"
end

export FESpaceConfig
export TimeConfig
export SimConfig

end # module ParameterHandler
