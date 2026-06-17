# ─────────────────────────────────────────────────────────────
# Frequency-domain simulate
# ─────────────────────────────────────────────────────────────


"""
    simulate(problem::HEFEM_Problem{FreqDomainConfig}) -> SimResult

Solve a fully assembled **frequency-domain** hydroelastic problem.

Builds or re-uses the `AffineFEOperator` stored in `problem`, solves the complex
linear system with the configured solver (default: `LUSolver()`), and returns the
solution wrapped in a `SimResult`.

# Arguments
- `problem::HEFEM_Problem{FreqDomainConfig}`: assembled problem container returned
  by `build_problem` with a `FreqDomainConfig` simulation config

# Returns
- `result::SimResult`: contains FE spaces, operator, and the complex-valued FE solution

# Example
```julia
import HydroElasticFEM.Simulation as S
problem = S.build_problem(domain, entities, config)
result  = S.simulate(problem)
phi_h, kappa_h = result.solution  # unpack multi-field solution
```

# Reference
[C23] Colomés et al. (2023), Int. J. Numer. Methods Eng., 124(3), 714-751.
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
    simulate(problem::HEFEM_Problem{TimeDomainConfig}, tconfig::TimeConfig) -> SimResult

Solve a fully assembled **time-domain** hydroelastic problem.

Builds a `GeneralizedAlpha2` ODE solver from `tconfig`, interpolates initial
conditions, and integrates the system over `[t₀, tf]`.  Returns the ODE solution
wrapped in a `SimResult`.

# Arguments
- `problem::HEFEM_Problem{TimeDomainConfig}`: assembled problem container returned
  by `build_problem` with a `TimeDomainConfig` simulation config
- `tconfig::TimeConfig`: time-stepping parameters (time step `Δt`, spectral radius
  `ρ∞`, initial conditions `u0`, `u0t`, `u0tt`, time window `t₀` to `tf`)

# Returns
- `result::SimResult`: contains FE spaces, transient operator, and ODE solution iterator

# Example
```julia
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.ParameterHandler as PH
tconfig  = PH.TimeConfig(Δt=0.01, tf=10.0)
problem  = S.build_problem(domain, entities, config)
result   = S.simulate(problem, tconfig)
for (t, u_h) in result.solution
    # process time step
end
```

# Reference
[C23] Colomés et al. (2023), Int. J. Numer. Methods Eng., 124(3), 714-751.
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
