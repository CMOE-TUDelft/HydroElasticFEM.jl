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

Loaded early in the module dependency chain so that both `Physics`
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
    space_type::DataType   = Float64
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
- `scheme::Symbol` — time integrator: `:generalized_alpha` (default, uses `ρ∞`)
  or `:newmark` (uses `γ`, `β`)
- `ρ∞::Float64` — spectral radius for Generalized-α (default 1.0)
- `γ::Float64`, `β::Float64` — Newmark parameters (default 0.5, 0.25; only
  used when `scheme = :newmark`)
- `αₕ` — stabilized free-surface parameter.  `nothing` (default) computes it
  automatically when the problem has damping zones (where it is required) and
  leaves the free surface unstabilised otherwise; `:auto` always computes it;
  a number is used as given.  The automatic value is
  `γ / (β Δt) / g · (1 - βₕ) / βₕ`, see [`time_integration_parameters`](@ref).
- `u0` — initial condition(s); tuple/vector of interpolatable objects per field.
  `nothing` starts from rest (all fields zero).
- `u0t` — initial velocity (optional; defaults to `u0`, or zero when `u0 = nothing`)
- `u0tt` — initial acceleration (optional; defaults to `u0`, or zero when `u0 = nothing`)
"""
@with_kw struct TimeConfig
    Δt::Float64
    t₀::Float64 = 0.0
    tf::Float64
    scheme::Symbol = :generalized_alpha
    ρ∞::Float64 = 1.0
    γ::Float64 = 0.5
    β::Float64 = 0.25
    αₕ::Union{Nothing, Float64, Symbol} = nothing
    u0 = nothing
    u0t = nothing
    u0tt = nothing

    @assert Δt > 0 "Δt must be positive"
    @assert tf > t₀ "tf must be greater than t₀"
    @assert scheme in (:generalized_alpha, :newmark) "scheme must be :generalized_alpha or :newmark"
    @assert αₕ isa Union{Nothing, Float64} || αₕ === :auto "αₕ must be nothing, a number, or :auto"
end

"""
    time_integration_parameters(tc::TimeConfig) -> (; γ, β, αf, αm)

Newmark-family parameters of the configured scheme, in the convention of
Gridap's `GeneralizedAlpha2` (`Newmark` is the case `αf = αm = 0`).  For
`:generalized_alpha`, `αf = ρ∞/(1+ρ∞)`, `αm = (2ρ∞-1)/(1+ρ∞)`,
`γ = 1/2 - αm + αf` and `β = (1 - αm + αf)²/4`.

For every member of the family, `∂u̇/∂u = γ/(β Δt)` at the new step, which
gives the free-surface stabilisation `αₕ = γ/(β Δt)/g · (1 - βₕ)/βₕ`.
"""
function time_integration_parameters(tc::TimeConfig)
    if tc.scheme === :newmark
        return (γ = tc.γ, β = tc.β, αf = 0.0, αm = 0.0)
    end
    ρ∞ = tc.ρ∞
    αf = ρ∞ / (1 + ρ∞)
    αm = (2 * ρ∞ - 1) / (1 + ρ∞)
    return (γ = 1 / 2 - αm + αf, β = (1 - αm + αf)^2 / 4, αf = αf, αm = αm)
end

"""
    stabilization_αₕ(tc::TimeConfig, g, βₕ) -> Float64

Free-surface stabilisation parameter `αₕ = γ/(β Δt)/g · (1 - βₕ)/βₕ` for the
configured time integrator.
"""
function stabilization_αₕ(tc::TimeConfig, g::Real, βₕ::Real)
    p = time_integration_parameters(tc)
    return p.γ / (p.β * tc.Δt) / g * (1 - βₕ) / βₕ
end

abstract type SimulationConfig end

"""
    FreqDomainConfig

Configuration for a frequency-domain simulation run.

# Fields
- `ω` — angular frequency (required when `domain == :frequency`)
- `solver` — optional solver override (e.g. `LUSolver()`)
"""
@with_kw struct FreqDomainConfig <: SimulationConfig
    ω::Union{Float64, Nothing} = nothing
    solver = nothing
end

"""
    TimeDomainConfig

Configuration for a time-domain simulation run.

# Fields
- `t₀` — start time (default 0.0)
- `tf` — final time (default 1.0)
- `solver` — optional solver override (e.g. `LUSolver()`)
"""
@with_kw struct TimeDomainConfig <: SimulationConfig
    t₀::Float64 = 0.0
    tf::Float64 = 1.0
    solver = nothing
end

export FESpaceConfig
export TimeConfig
export time_integration_parameters, stabilization_αₕ
export FreqDomainConfig
export TimeDomainConfig
export SimulationConfig

end # module ParameterHandler
