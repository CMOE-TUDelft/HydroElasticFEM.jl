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
- `solver` — solver object
"""
struct HEFEM_Problem
    model::DiscreteModel
    triangulations::G.TankTriangulations
    integration_domains::G.IntegrationDomains
    entities::Vector{P.PhysicsParameters}  # physics entities (e.g. Membrane2D, FreeSurface, etc.)
    test_fe_space::MultiFieldFESpace
    trial_fe_space::MultiFieldFESpace
    fe_operator::FEOperator
    solver::FESolver
end

# Getter functions for HEFEM_Problem fields
get_model(prob::HEFEM_Problem) = prob.model
get_triangulations(prob::HEFEM_Problem) = prob.triangulations
get_integration_domains(prob::HEFEM_Problem) = prob.integration_domains
get_entities(prob::HEFEM_Problem) = prob.entities
get_test_fe_space(prob::HEFEM_Problem) = prob.test_fe_space
get_trial_fe_space(prob::HEFEM_Problem) = prob.trial_fe_space
get_fe_operator(prob::HEFEM_Problem) = prob.fe_operator
get_solver(prob::HEFEM_Problem) = prob.solver

"""
    build_problem(config, entities_trians...; dom, ...)

Given a `SimConfig` and tuples of physics entities + triangulations,
constructs the `HEFEM_Problem` by building FE spaces, assembling the FE operator, and building the solver.
"""
function build_problem(domain, physics::Vector{P.PhysicsParameters}; tconfig=nothing)
    
    # Build discrete model and triangulations from geometry
    model = G.build_model(domain)
    trians = G.build_triangulations(domain, model)
    measures = G.get_integration_domains(trians, degree=2)

    # Build FE spaces
    X, Y, fmap = FA.build_fe_spaces(physics, trians; transient=!isnothing(tconfig))

    # Build FE Operator
    op = build_fe_operator(physics, measures, fmap, X, Y)

    # Build solver
    solver = LUSolver()

    HEFEM_Problem(model, trians, measures, physics, Y, X, op, solver)

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
