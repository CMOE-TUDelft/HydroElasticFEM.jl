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
