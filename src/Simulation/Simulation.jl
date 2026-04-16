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
    entities::Vector{P.PhysicsParameters}  # physics entities (e.g. Membrane2D, FreeSurface, etc.)
    field_map::Dict{Symbol,Int}  # maps physics entity symbols to FE space indices
    test_fe_space::MultiFieldFESpace
    trial_fe_space::MultiFieldFESpace
    fe_operator::Union{FEOperator,TransientFEOperator}
    sim_config::T
end

# Getter functions for HEFEM_Problem fields
get_model(prob::HEFEM_Problem) = prob.model
get_triangulations(prob::HEFEM_Problem) = prob.triangulations
get_integration_domains(prob::HEFEM_Problem) = prob.integration_domains
"""
    get_assembly_context(prob::HEFEM_Problem)

Return the immutable assembly context attached to a built problem.
This is either a `FrequencyAssemblyContext` or `TimeAssemblyContext`.
"""
get_assembly_context(prob::HEFEM_Problem) = prob.assembly_context
get_entities(prob::HEFEM_Problem) = prob.entities
get_field_map(prob::HEFEM_Problem) = prob.field_map
get_test_fe_space(prob::HEFEM_Problem) = prob.test_fe_space
get_trial_fe_space(prob::HEFEM_Problem) = prob.trial_fe_space
get_fe_operator(prob::HEFEM_Problem) = prob.fe_operator
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

"""
    build_frequency_context(domains, physics, config)

Build a `FrequencyAssemblyContext` from integration domains and
frequency-domain simulation inputs. The function also triggers
radiation-boundary validation for enabled `RadiationBC` entries.
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
    build_time_context(domains, physics, config, tconfig)

Build a `TimeAssemblyContext` from integration domains and time-domain
simulation inputs. For damping-zone cases, a non-`nothing` `tconfig`
with `αₕ` is required.
"""
function build_time_context(domains::G.IntegrationDomains,
                            physics::Vector{P.PhysicsParameters},
                            config::PH.TimeDomainConfig,
                            tconfig)
    if _has_damping_zone_bc(physics)
        isnothing(tconfig) && error("Time-domain damping-zone problems require `tconfig` to be passed to `build_problem`.")
        isnothing(tconfig.αₕ) && error("Time-domain damping-zone problems require `TimeConfig.αₕ`.")
    end
    t₀ = isnothing(tconfig) ? config.t₀ : tconfig.t₀
    αₕ = isnothing(tconfig) ? nothing : tconfig.αₕ
    AC.TimeAssemblyContext(domains, t₀, αₕ)
end

function _build_problem_parts(domain, physics::Vector{P.PhysicsParameters}, config::PH.SimulationConfig)
    model = G.build_model(domain)
    trians = G.build_triangulations(domain, model)
    degrees = get_integration_degrees(trians, physics)
    measures = G.get_integration_domains(trians, degree=degrees)
    X, Y, fmap = FA.build_fe_spaces(physics, trians, config)
    return model, trians, measures, X, Y, fmap
end

"""
    build_problem(config, entities_trians...; dom, ...)

Given a simulation config and physics entities, construct the immutable
assembly context, FE spaces, and FE operator.
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
    get_integration_degrees(trians::G.TankTriangulations, physics::Vector{P.PhysicsParameters})

Determines the required integration degree for each domain based on the FE order of the physics entities 
defined on that domain. Returns a dictionary mapping domain symbols to integration degrees.

# Arguments
- `trians` — triangulations for each domain (used to get domain symbols)
- `physics` — vector of physics entities, each with an optional `fe` field containing
"""
function get_integration_degrees(trians::G.TankTriangulations, physics::Vector{P.PhysicsParameters})
    
    # Get max FE order across all entities for each domain
    degrees = Dict{Symbol, Int}()
    for p in physics
        fe = p.fe
        if fe !== nothing
            domain_symbol = p.space_domain_symbol
            if haskey(degrees, domain_symbol)
                degrees[domain_symbol] = max(degrees[domain_symbol], fe.order)
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
