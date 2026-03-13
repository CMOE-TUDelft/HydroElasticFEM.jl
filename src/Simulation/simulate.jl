# ─────────────────────────────────────────────────────────────
# Frequency-domain simulate
# ─────────────────────────────────────────────────────────────


"""
    simulate(config::SimConfig, entities, trians;
             dom::IntegrationDomains, couplings=nothing, rhs_fn=nothing)

Run a **frequency-domain** simulation.

# Arguments
- `config::SimConfig` — must have `domain == :frequency` and `ω` set
- `entities` — vector of entities (PhysicsParameters or Vector{ResonatorSingle})
- `trians` — `G.TankTriangulations`
- `dom::IntegrationDomains` — integration measures/normals
- `couplings` — optional explicit coupling pairs; auto-detected if `nothing`
- `rhs_fn` — optional callable `rhs_fn(y::FieldMap) -> DomainContribution`;
  if `nothing`, uses zero right-hand side (requires `:dΩ` in `dom`)

# Returns
`SimResult` containing FE spaces, operator, and solution.
"""
function simulate(config::SimConfig, 
                  entities::Vector{<:P.PhysicsParameters}, 
                  trians::G.TankTriangulations;
                  dom::G.IntegrationDomains,
                  couplings=nothing, rhs_fn=nothing)
    config.domain == :frequency ||
        throw(ArgumentError("Use the (config, tconfig, ...) method for time-domain"))

    X, Y, fmap = FA.build_fe_spaces(entities, trians)

    ω = config.ω

    op = build_fe_operator(entities, dom, ω, fmap, X, Y;
                           rhs_fn=rhs_fn)

    solver = isnothing(config.solver) ? LUSolver() : config.solver
    solution = solve(solver, op)

    SimResult(X, Y, fmap, op, solution)
end

# ─────────────────────────────────────────────────────────────
# Time-domain simulate
# ─────────────────────────────────────────────────────────────


"""
    simulate(config::SimConfig, tconfig::TimeConfig, entities, trians;
             dom::IntegrationDomains, couplings=nothing, rhs_fn=nothing)

Run a **time-domain** simulation.

# Arguments
- `config::SimConfig` — must have `domain == :time`
- `tconfig::TimeConfig` — time-stepping and initial-condition parameters
- `entities` — vector of entities (PhysicsParameters or Vector{ResonatorSingle})
- `trians` — `G.TankTriangulations`
- `dom::IntegrationDomains` — integration measures/normals
- `couplings` — optional explicit coupling pairs; auto-detected if `nothing`
- `rhs_fn` — optional callable `rhs_fn(t, y::FieldMap) -> DomainContribution`;
  if `nothing`, uses zero right-hand side (requires `:dΩ` in `dom`)

# Returns
`SimResult` containing FE spaces, transient operator, and ODE solution.
"""
function simulate(config::SimConfig, tconfig::TimeConfig,
                  entities::Vector{<:P.PhysicsParameters}, 
                  trians::G.TankTriangulations;
                  dom::G.IntegrationDomains,
                  couplings=nothing, rhs_fn=nothing)
    config.domain == :time ||
        throw(ArgumentError("Use the (config, ...) method for frequency-domain"))

    X, Y, fmap = FA.build_fe_spaces(entities, trians; transient=true)

    op = build_fe_operator(entities, dom, fmap, X, Y;
                           rhs_fn=rhs_fn)

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
