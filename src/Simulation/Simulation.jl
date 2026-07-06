"""
    module Simulation

High-level simulation orchestrator for HydroElasticFEM.

Given physics entities, triangulations and integration domains, builds
FE spaces, assembles the FE operator (frequency or time-domain), solves,
and returns a `SimResult`.

# Main entry point
- `simulate(config, entities_trians...; dom, ...)` — frequency-domain
- `simulate(config, tconfig, entities_trians...; dom, ...)` — time-domain
"""
module Simulation

using Parameters
using Gridap
using Gridap.ODEs

import ..Geometry as G
import ..ParameterHandler as PH
using ..ParameterHandler  # brings SimConfig, TimeConfig into scope
import ..Physics as P
import ..AssemblyContexts as AC


# FESpaceAssembly (build_fe_spaces, build_test/trial_fe_space)
include("FESpaceAssembly.jl")
using .FESpaceAssembly
const FA = FESpaceAssembly

# FEOperators (FieldMap, assemble_*, detect_couplings)
include("FEOperators.jl")
using .FEOperators

"""
    HEFEM_Problem

Container for all built entities in a simulation.

# Fields
- `geometry` — model and domain configuration
- `triangulations` — mesh/triangulation objects
- `integration_domains` — integration domains (measures, normals, etc.)
- `entities` — list of physics entities
- `fe_spaces` — finite element spaces (test/trial, field map, etc.)
- `fe_operator` — assembled FE operator
- `sim_config` — simulation configuration (`FreqDomainConfig` or `TimeDomainConfig`)
"""
struct HEFEM_Problem{T<:PH.SimulationConfig,C<:AC.AbstractAssemblyContext}
    model::DiscreteModel
    triangulations::G.TankTriangulations
    integration_domains::G.IntegrationDomains
    assembly_context::C
    entities::Vector{P.PhysicsParameters}  # physics entities (e.g. Membrane, FreeSurface, etc.)
    field_map::Dict{Symbol,Int}  # maps physics entity symbols to FE space indices
    test_fe_space::MultiFieldFESpace
    trial_fe_space::MultiFieldFESpace
    fe_operator::Union{FEOperator,TransientFEOperator}
    sim_config::T
end

"""
    get_model(prob::HEFEM_Problem) -> DiscreteModel

Return the discrete mesh model from a built problem.
"""
get_model(prob::HEFEM_Problem) = prob.model

"""
    get_triangulations(prob::HEFEM_Problem) -> TankTriangulations

Return the triangulations container from a built problem.
"""
get_triangulations(prob::HEFEM_Problem) = prob.triangulations

"""
    get_integration_domains(prob::HEFEM_Problem) -> IntegrationDomains

Return the integration domains (measures and normals) from a built problem.
"""
get_integration_domains(prob::HEFEM_Problem) = prob.integration_domains
"""
    get_assembly_context(prob::HEFEM_Problem)

Return the immutable assembly context attached to a built problem.
This is either a `FrequencyAssemblyContext` or `TimeAssemblyContext`.
"""
get_assembly_context(prob::HEFEM_Problem) = prob.assembly_context
"""
    get_entities(prob::HEFEM_Problem) -> Vector{PhysicsParameters}

Return the ordered list of physics entities from a built problem.
"""
get_entities(prob::HEFEM_Problem) = prob.entities

"""
    get_field_map(prob::HEFEM_Problem) -> Dict{Symbol,Int}

Return the symbol-to-field-index map from a built problem.
"""
get_field_map(prob::HEFEM_Problem) = prob.field_map

"""
    get_test_fe_space(prob::HEFEM_Problem) -> MultiFieldFESpace

Return the test `MultiFieldFESpace` from a built problem.
"""
get_test_fe_space(prob::HEFEM_Problem) = prob.test_fe_space

"""
    get_trial_fe_space(prob::HEFEM_Problem) -> MultiFieldFESpace

Return the trial `MultiFieldFESpace` from a built problem.
"""
get_trial_fe_space(prob::HEFEM_Problem) = prob.trial_fe_space

"""
    get_fe_operator(prob::HEFEM_Problem) -> Union{FEOperator,TransientFEOperator}

Return the assembled FE operator from a built problem.
"""
get_fe_operator(prob::HEFEM_Problem) = prob.fe_operator

"""
    get_sim_config(prob::HEFEM_Problem) -> SimulationConfig

Return the simulation configuration (`FreqDomainConfig` or `TimeDomainConfig`)
from a built problem.
"""
get_sim_config(prob::HEFEM_Problem) = prob.sim_config

function _find_free_surface(physics::Vector{P.PhysicsParameters})
    idx = findfirst(p -> p isa P.FreeSurface, physics)
    isnothing(idx) ? nothing : physics[idx]
end

function _has_damping_zone_bc(physics::Vector{P.PhysicsParameters})
    any(
        p -> p isa P.PotentialFlow && any(bc -> bc isa P.DampingZoneBC && bc.enabled, p.boundary_conditions),
        physics,
    )
end

function _has_radiation_bc(physics::Vector{P.PhysicsParameters})
    any(
        p -> p isa P.PotentialFlow && any(bc -> bc isa P.RadiationBC && bc.enabled, p.boundary_conditions),
        physics,
    )
end

"""
    build_frequency_context(domains, physics, config) -> FrequencyAssemblyContext

Build a `FrequencyAssemblyContext` from integration domains and
frequency-domain simulation inputs.

Also validates that any enabled `RadiationBC` entries have a compatible
frequency value via `_radiation_frequency`.

# Arguments
- `domains::IntegrationDomains` — integration domains (measures and normals)
- `physics::Vector{PhysicsParameters}` — physics entity list
- `config::FreqDomainConfig` — frequency-domain configuration (provides `ω`)
"""
function build_frequency_context(domains::G.IntegrationDomains,
                                 physics::Vector{P.PhysicsParameters},
                                 config::PH.FreqDomainConfig)
    for entity in physics
        if entity isa P.PotentialFlow
            for bc in entity.boundary_conditions
                if bc isa P.RadiationBC && bc.enabled
                    P._radiation_frequency(entity)
                end
            end
        end
    end

    free_surface = _find_free_surface(physics)
    αₕ = isnothing(free_surface) ? nothing :
        -im * config.ω / free_surface.g * (1.0 - free_surface.βₕ) / free_surface.βₕ
    AC.FrequencyAssemblyContext(domains, config.ω, αₕ)
end

"""
    build_time_context(domains, physics, config, tconfig) -> TimeAssemblyContext

Build a `TimeAssemblyContext` from integration domains and time-domain
simulation inputs.

Errors if a `DampingZoneBC` is present but `tconfig` is `nothing` or
`tconfig.αₕ` is `nothing`.

# Arguments
- `domains::IntegrationDomains` — integration domains (measures and normals)
- `physics::Vector{PhysicsParameters}` — physics entity list
- `config::TimeDomainConfig` — time-domain configuration (provides `t₀`)
- `tconfig` — `TimeConfig` (or `nothing`); provides `t₀` and `αₕ` overrides
"""
function build_time_context(domains::G.IntegrationDomains,
                            physics::Vector{P.PhysicsParameters},
                            config::PH.TimeDomainConfig,
                            tconfig)
    if _has_damping_zone_bc(physics)
        isnothing(tconfig) && error("Time-domain damping-zone problems require `tconfig` to be passed to `build_problem`.")
        isnothing(tconfig.αₕ) && error("Time-domain damping-zone problems require `TimeConfig.αₕ`.")
    end

    if _has_radiation_bc(physics)
        for entity in physics
            if entity isa P.PotentialFlow
                for bc in entity.boundary_conditions
                    if bc isa P.RadiationBC && bc.enabled
                        P._radiation_frequency(entity)
                    end
                end
            end
        end
    end

    t₀ = isnothing(tconfig) ? config.t₀ : tconfig.t₀
    αₕ = isnothing(tconfig) ? nothing : tconfig.αₕ
    AC.TimeAssemblyContext(domains, t₀, αₕ)
end

function _entity_ambient_dimension(entity::P.PhysicsParameters)
    if isdefined(P, :ambient_dimension) && applicable(P.ambient_dimension, entity)
        return P.ambient_dimension(entity)
    end
    if applicable(G.ambient_dimension, entity)
        return G.ambient_dimension(entity)
    end
    return nothing
end

function _check_ambient_dimension_consistency(domain,
                                              physics::Vector{P.PhysicsParameters})
    ddom = G.ambient_dimension(domain)
    for entity in physics
        dent = _entity_ambient_dimension(entity)
        isnothing(dent) && continue
        if dent != ddom
            error(
                "Ambient-dimension mismatch: domain has dim=$(ddom), " *
                "but $(typeof(entity)) has dim=$(dent).",
            )
        end
    end
    nothing
end

function _build_problem_parts(domain, physics::Vector{P.PhysicsParameters}, config::PH.SimulationConfig)
    _check_ambient_dimension_consistency(domain, physics)
    model = G.build_model(domain)
    trians = G.build_triangulations(domain, model)
    degrees = get_integration_degrees(trians, physics)
    measures = G.get_integration_domains(trians, degree=degrees)
    X, Y, fmap = FA.build_fe_spaces(physics, trians, config)
    return model, trians, measures, X, Y, fmap
end

"""
    build_problem(domain, physics, config; rhs_fn=nothing) -> HEFEM_Problem
    build_problem(domain, physics, config; tconfig=nothing, rhs_fn=nothing) -> HEFEM_Problem

Construct an immutable `HEFEM_Problem` from a domain descriptor, a list of
physics entities, and a simulation configuration.

This is the primary entry point for building a HydroElasticFEM simulation.
It performs, in order:
1. Ambient-dimension consistency check between `domain` and all entities.
2. `build_model(domain)` — create the discrete mesh.
3. `build_triangulations(domain, model)` — extract all triangulations.
4. `get_integration_domains(trians)` — build `Measure`s and normals.
5. `build_fe_spaces(physics, trians, config)` — build `MultiFieldFESpace`.
6. `build_frequency_context` or `build_time_context` — set up the assembly context.
7. `build_frequency_fe_operator` or `build_time_fe_operator` — assemble the operator.

# Arguments
- `domain` — geometry descriptor (`TankDomain`, `CartesianDomain`, or `GmshDomain`)
- `physics::Vector{PhysicsParameters}` — ordered list of physics entities
  (e.g. `[fluid, fs, beam]`).  Order determines the `MultiFieldFESpace` layout.
- `config::FreqDomainConfig` or `config::TimeDomainConfig` — simulation config
- `tconfig::TimeConfig` — (time-domain only) time integration configuration
  including `t₀` and optional `αₕ` for damping-zone stabilisation
- `rhs_fn` — optional user-supplied forcing function

# Returns
`HEFEM_Problem` containing the model, triangulations, integration domains,
assembly context, FE spaces, FE operator, and simulation config.

# Example
```julia
tank = TankDomain(L=20.0, H=5.0, nx=40, ny=4)
config = FreqDomainConfig(ω=1.0)
fluid = PotentialFlow()
fs    = FreeSurface()
prob  = build_problem(tank, [fluid, fs], config)
```

See also: [`simulate`](@ref), [`HEFEM_Problem`](@ref)
"""
function build_problem(domain, physics::Vector{P.PhysicsParameters}, config::PH.FreqDomainConfig; rhs_fn=nothing)
    model, trians, measures, X, Y, fmap = _build_problem_parts(domain, physics, config)
    ctx = build_frequency_context(measures, physics, config)
    op = build_frequency_fe_operator(physics, ctx, fmap, X, Y; rhs_fn=rhs_fn)
    HEFEM_Problem(model, trians, measures, ctx, physics, fmap, Y, X, op, config)
end

function build_problem(domain, physics::Vector{P.PhysicsParameters}, config::PH.TimeDomainConfig; tconfig=nothing, rhs_fn=nothing)
    model, trians, measures, X, Y, fmap = _build_problem_parts(domain, physics, config)
    ctx = build_time_context(measures, physics, config, tconfig)
    op = build_time_fe_operator(physics, ctx, fmap, X, Y; rhs_fn=rhs_fn)
    HEFEM_Problem(model, trians, measures, ctx, physics, fmap, Y, X, op, config)
end


"""
    get_integration_degrees(trians, physics) -> Dict{Symbol,Int}

Determine the required Gaussian integration degree for each triangulation
domain, based on the polynomial order declared in each entity's `fe` field.

For each entity, the degree for its `space_domain_symbol` domain is set to
`2 * fe.order` on first encounter, and to the maximum of existing and new
degree on subsequent encounters.  Any domain present in `trians` but not
covered by a physics entity defaults to degree `2`.

# Arguments
- `trians::TankTriangulations` — triangulations container (used to enumerate domains)
- `physics::Vector{PhysicsParameters}` — ordered list of physics entities

# Returns
`Dict{Symbol,Int}` mapping triangulation-domain symbols to integration degrees.
"""
function get_integration_degrees(trians::G.TankTriangulations, physics::Vector{P.PhysicsParameters})
    
    # Get max FE order across all entities for each domain
    degrees = Dict{Symbol, Int}()
    for p in physics
        fe = p.fe
        if fe !== nothing
            domain_symbol = p.space_domain_symbol
            if haskey(degrees, domain_symbol)
                degrees[domain_symbol] = max(degrees[domain_symbol], 2*fe.order)
            else
                degrees[domain_symbol] = 2*fe.order
            end
        end
    end

    # Add default degree for any domains not covered by physics entities
    for (name, trian) in pairs(trians.data)
        if !haskey(degrees, name)
            degrees[name] = 2  # default degree
        end
    end

    return degrees
    
end


include("SimResult.jl")
include("simulate.jl")

export SimResult
export simulate
export build_problem
export build_frequency_context, build_time_context
export HEFEM_Problem
export FESpaceAssembly, FEOperators
export build_fe_spaces, build_test_fe_space, build_trial_fe_space
export FieldMap, detect_couplings, build_fe_operator
export build_frequency_fe_operator, build_time_fe_operator
export get_assembly_context
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module Simulation
