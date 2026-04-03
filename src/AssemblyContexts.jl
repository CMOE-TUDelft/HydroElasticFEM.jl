module AssemblyContexts

import ..Geometry as G

"""
    AbstractAssemblyContext

Abstract supertype for context objects passed to element and system assembly routines.
Contexts carry analysis-dependent metadata such as active domains, time/frequency
coordinates, and optional stabilization parameters.
"""
abstract type AbstractAssemblyContext end

"""
    FrequencyAssemblyContext(domains, ω, αₕ)

Assembly context for frequency-domain problems.

- `domains`: geometric/physical domains involved in the integration.
- `ω`: angular frequency.
- `αₕ`: optional stabilization parameter (`nothing` when disabled).
"""
struct FrequencyAssemblyContext{D,O,A} <: AbstractAssemblyContext
    domains::D
    ω::O
    αₕ::A
end

"""
    TimeAssemblyContext(domains, t, αₕ)

Assembly context for time-domain problems.

- `domains`: geometric/physical domains involved in the integration.
- `t`: current physical time.
- `αₕ`: optional stabilization parameter (`nothing` when disabled).
"""
struct TimeAssemblyContext{D,T,A} <: AbstractAssemblyContext
    domains::D
    t::T
    αₕ::A
end

"""
    domains(ctx)

Return the domain collection used in integration and stored in an 
assembly context.
"""
domains(ctx::AbstractAssemblyContext) = ctx.domains

"""
    frequency(ctx)

Return the angular frequency `ω` from a frequency-domain context.
"""
frequency(ctx::FrequencyAssemblyContext) = ctx.ω

"""
    current_time(ctx)

Return the current physical time `t` from a time-domain context.
"""
current_time(ctx::TimeAssemblyContext) = ctx.t

"""
    stabilization_parameter(ctx)

Return the stabilization parameter `αₕ` stored in the context. This can be
`nothing` when stabilization is disabled.
"""
stabilization_parameter(ctx::AbstractAssemblyContext) = ctx.αₕ

"""
    has_stabilization(ctx)

Return `true` when the context carries a non-`nothing` stabilization parameter.
"""
has_stabilization(ctx::AbstractAssemblyContext) = !isnothing(stabilization_parameter(ctx))

"""
    with_time(ctx, t)

Create a new `TimeAssemblyContext` from `ctx` with updated time `t`, preserving
domains and stabilization settings.
"""
with_time(ctx::TimeAssemblyContext, t) = TimeAssemblyContext(domains(ctx), t, stabilization_parameter(ctx))

export AbstractAssemblyContext
export FrequencyAssemblyContext, TimeAssemblyContext
export domains, frequency, current_time
export stabilization_parameter, has_stabilization, with_time

end
