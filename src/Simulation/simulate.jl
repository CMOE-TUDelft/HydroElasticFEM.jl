# ─────────────────────────────────────────────────────────────
# Coupling detection
# ─────────────────────────────────────────────────────────────

"""
    detect_couplings(entities) -> Vector{Tuple}

Auto-detect active coupling pairs among `entities` by checking all
ordered pairs `(a, b)` against the trait functions `has_mass_form`,
`has_damping_form`, `has_stiffness_form`, and `has_rhs_form`.

Returns a vector of `(entity_a, entity_b)` tuples that have at
least one active coupling form.
"""
function detect_couplings(entities)
    pairs = Tuple[]
    for (i, ea) in enumerate(entities)
        for (j, eb) in enumerate(entities)
            i == j && continue
            if (E.has_mass_form(ea, eb) || E.has_damping_form(ea, eb) ||
                E.has_stiffness_form(ea, eb) || E.has_rhs_form(ea, eb))
                push!(pairs, (ea, eb))
            end
        end
    end
    return pairs
end

# ─────────────────────────────────────────────────────────────
# Internal assembly helpers
# ─────────────────────────────────────────────────────────────

_add(a, b) = isnothing(a) ? b : a + b

"""Sum single-entity + coupling contributions for a given form."""
function _assemble_form(f, has_f, entities, coupling_pairs, dom, fmap, x, y)
    xd = FO.FieldMap(x, fmap)
    yd = FO.FieldMap(y, fmap)
    val = nothing
    for e in entities
        if has_f(e)
            val = _add(val, f(e, dom, xd, yd))
        end
    end
    for (ea, eb) in coupling_pairs
        if has_f(ea, eb)
            val = _add(val, f(ea, eb, dom, xd, yd))
        end
    end
    isnothing(val) && error("No active $(nameof(f)) contributions found")
    return val
end

_has_weakform(e) =
    E.has_mass_form(e) || E.has_damping_form(e) || E.has_stiffness_form(e)
_has_weakform(a, b) =
    E.has_mass_form(a, b) || E.has_damping_form(a, b) || E.has_stiffness_form(a, b)

"""Assemble frequency-domain bilinear form (single + coupling entities)."""
function _assemble_bilinear(entities, coupling_pairs, dom, ω, fmap, x, y)
    xd = FO.FieldMap(x, fmap)
    yd = FO.FieldMap(y, fmap)
    val = nothing
    for e in entities
        if _has_weakform(e)
            val = _add(val, E.weakform(e, dom, ω, xd, yd))
        end
    end
    for (ea, eb) in coupling_pairs
        if _has_weakform(ea, eb)
            val = _add(val, E.weakform(ea, eb, dom, ω, xd, yd))
        end
    end
    isnothing(val) && error("No active weak form contributions found")
    return val
end

# ─────────────────────────────────────────────────────────────
# Frequency-domain simulate
# ─────────────────────────────────────────────────────────────

"""
    simulate(config::SimConfig, entities_trians::Pair...;
             dom::IntegrationDomains, couplings=nothing, rhs_fn=nothing)

Run a **frequency-domain** simulation.

# Arguments
- `config::SimConfig` — must have `domain == :frequency` and `ω` set
- `entities_trians::Pair...` — sequence of `entity => triangulation`
- `dom::IntegrationDomains` — integration measures/normals
- `couplings` — optional explicit coupling pairs; auto-detected if `nothing`
- `rhs_fn` — optional callable `rhs_fn(y::FieldMap) -> DomainContribution`;
  if `nothing`, uses zero right-hand side (requires `:dΩ` in `dom`)

# Returns
`SimResult` containing FE spaces, operator, and solution.
"""
function simulate(config::SimConfig, entities_trians::Pair...;
                  dom::G.IntegrationDomains,
                  couplings=nothing, rhs_fn=nothing)
    config.domain == :frequency ||
        throw(ArgumentError("Use the (config, tconfig, ...) method for time-domain"))

    entities = [first(p) for p in entities_trians]
    X, Y, fmap = FA.build_fe_spaces(entities_trians...)

    coupling_pairs = isnothing(couplings) ? detect_couplings(entities) : couplings
    ω = config.ω

    a(x, y) = _assemble_bilinear(entities, coupling_pairs, dom, ω, fmap, x, y)

    l = if rhs_fn !== nothing
        y -> rhs_fn(FO.FieldMap(y, fmap))
    else
        first_sym = first(keys(fmap))
        y -> ∫(0 * FO.FieldMap(y, fmap)[first_sym])dom[:dΩ]
    end

    op = AffineFEOperator(a, l, X, Y)
    solver = isnothing(config.solver) ? LUSolver() : config.solver
    solution = solve(solver, op)

    SimResult(X, Y, fmap, op, solution)
end

# ─────────────────────────────────────────────────────────────
# Time-domain simulate
# ─────────────────────────────────────────────────────────────

"""
    simulate(config::SimConfig, tconfig::TimeConfig, entities_trians::Pair...;
             dom::IntegrationDomains, couplings=nothing, rhs_fn=nothing)

Run a **time-domain** simulation.

# Arguments
- `config::SimConfig` — must have `domain == :time`
- `tconfig::TimeConfig` — time-stepping and initial-condition parameters
- `entities_trians::Pair...` — sequence of `entity => triangulation`
- `dom::IntegrationDomains` — integration measures/normals
- `couplings` — optional explicit coupling pairs; auto-detected if `nothing`
- `rhs_fn` — optional callable `rhs_fn(t, y::FieldMap) -> DomainContribution`;
  if `nothing`, uses zero right-hand side (requires `:dΩ` in `dom`)

# Returns
`SimResult` containing FE spaces, transient operator, and ODE solution.
"""
function simulate(config::SimConfig, tconfig::TimeConfig,
                  entities_trians::Pair...;
                  dom::G.IntegrationDomains,
                  couplings=nothing, rhs_fn=nothing)
    config.domain == :time ||
        throw(ArgumentError("Use the (config, ...) method for frequency-domain"))

    entities = [first(p) for p in entities_trians]
    X, Y, fmap = FA.build_fe_spaces(entities_trians...; transient=true)

    coupling_pairs = isnothing(couplings) ? detect_couplings(entities) : couplings

    # Time-domain decomposed forms: stiffness, damping, mass
    a(t, x, y)    = _assemble_form(E.stiffness, E.has_stiffness_form,
                                    entities, coupling_pairs, dom, fmap, x, y)
    c(t, x_t, y)  = _assemble_form(E.damping, E.has_damping_form,
                                    entities, coupling_pairs, dom, fmap, x_t, y)
    m(t, x_tt, y) = _assemble_form(E.mass, E.has_mass_form,
                                    entities, coupling_pairs, dom, fmap, x_tt, y)

    l = if rhs_fn !== nothing
        (t, y) -> rhs_fn(t, FO.FieldMap(y, fmap))
    else
        first_sym = first(keys(fmap))
        (t, y) -> ∫(0 * FO.FieldMap(y, fmap)[first_sym])dom[:dΩ]
    end

    op = TransientLinearFEOperator((a, c, m), l, X, Y;
        constant_forms=(true, true, true))

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
