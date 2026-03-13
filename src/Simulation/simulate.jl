# ─────────────────────────────────────────────────────────────
# Frequency-domain simulate
# ─────────────────────────────────────────────────────────────


"""
    simulate(config::SimConfig, entities, trians;
             dom::IntegrationDomains, couplings=nothing, rhs_fn=nothing)

Run a **frequency-domain** simulation.

# Arguments
- `problem::HEFEM_Problem{PH.FreqDomainConfig}` — problem container with FE spaces, operator, and other entities

# Returns
`SimResult` containing FE spaces, operator, and solution.
"""
function simulate(problem::HEFEM_Problem{PH.FreqDomainConfig})

    X = get_test_fe_space(problem)
    Y = get_trial_fe_space(problem)
    fmap = get_field_map(problem)
    op = get_fe_operator(problem)
    config = get_sim_config(problem)

    solver = isnothing(config.solver) ? LUSolver() : config.solver
    solution = solve(solver, op)

    SimResult(X, Y, fmap, op, solution)
end

# ─────────────────────────────────────────────────────────────
# Time-domain simulate
# ─────────────────────────────────────────────────────────────


"""
    simulate(config::SimConfig, tconfig::TimeConfig, problem::HEFEM_Problem)

Run a **time-domain** simulation.

# Arguments
- `problem::HEFEM_Problem` — problem container with FE spaces, operator, and other entities
- `tconfig::TimeConfig` — time-stepping and initial-condition parameters

# Returns
`SimResult` containing FE spaces, transient operator, and ODE solution.
"""
function simulate(problem::HEFEM_Problem{PH.TimeDomainConfig}, tconfig::TimeConfig)

    X = get_test_fe_space(problem)
    Y = get_trial_fe_space(problem)
    fmap = get_field_map(problem)
    op = get_fe_operator(problem)
    config = get_sim_config(problem)

    ls = isnothing(config.solver) ? LUSolver() : config.solver
    ode_solver = GeneralizedAlpha2(ls, tconfig.Δt, tconfig.ρ∞)

    # Interpolate initial conditions
    X0 = X(tconfig.t₀)
    u0   = interpolate_everywhere(tconfig.u0, X0)
    u0t  = isnothing(tconfig.u0t) ? interpolate_everywhere(tconfig.u0, X0) : interpolate_everywhere(tconfig.u0t, X0)
    u0tt = isnothing(tconfig.u0tt) ? interpolate_everywhere(tconfig.u0, X0) : interpolate_everywhere(tconfig.u0tt, X0)

    solution = solve(ode_solver, op, tconfig.t₀, tconfig.tf, (u0, u0t, u0tt))

    SimResult(X, Y, fmap, op, solution)
end
