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
