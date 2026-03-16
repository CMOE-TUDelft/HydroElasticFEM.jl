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
struct HEFEM_Problem{T<:PH.SimulationConfig}
    model::DiscreteModel
    triangulations::G.TankTriangulations
    integration_domains::G.IntegrationDomains
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
get_entities(prob::HEFEM_Problem) = prob.entities
get_field_map(prob::HEFEM_Problem) = prob.field_map
get_test_fe_space(prob::HEFEM_Problem) = prob.test_fe_space
get_trial_fe_space(prob::HEFEM_Problem) = prob.trial_fe_space
get_fe_operator(prob::HEFEM_Problem) = prob.fe_operator
get_sim_config(prob::HEFEM_Problem) = prob.sim_config

"""
    build_problem(config, entities_trians...; dom, ...)

Given a `SimConfig` and tuples of physics entities + triangulations,
constructs the `HEFEM_Problem` by building FE spaces, assembling the FE operator, and building the solver.
"""
function build_problem(domain, physics::Vector{P.PhysicsParameters}, config::PH.SimulationConfig; tconfig=nothing, rhs_fn=nothing)
    
    # Build discrete model and triangulations from geometry
    model = G.build_model(domain)
    trians = G.build_triangulations(domain, model)
    degrees = get_integration_degrees(trians, physics)
    measures = G.get_integration_domains(trians, degree=degrees)

    # Build FE spaces
    X, Y, fmap = FA.build_fe_spaces(physics, trians; transient=!isnothing(tconfig))

    # Build FE Operator
    if isa(config, PH.FreqDomainConfig)
        op = build_fe_operator(physics, measures, config.ω, fmap, X, Y; rhs_fn=rhs_fn)
    elseif isa(config, PH.TimeDomainConfig)
        op = build_fe_operator(physics, measures, fmap, X, Y; rhs_fn=rhs_fn)
    end

    HEFEM_Problem(model, trians, measures, physics, fmap, Y, X, op, config)

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
export HEFEM_Problem
export FESpaceAssembly, FEOperators
export build_fe_spaces, build_test_fe_space, build_trial_fe_space
export FieldMap, detect_couplings, build_fe_operator
export assemble_weakform
export assemble_mass, assemble_damping, assemble_stiffness, assemble_rhs
export assemble_residual, assemble_jacobian, assemble_jacobian_t, assemble_jacobian_tt

end # module Simulation
