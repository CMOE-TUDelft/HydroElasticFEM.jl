
"""
    materialize(param, trian)

Convert `param` to something Gridap can use in a weak form on `trian`.

- If `param` is a `Float64` (or any scalar), it is returned as-is.
- If `param` is a `Function`, it is wrapped in a `CellField` on `trian`.

The triangulation is typically obtained from the test function via
`get_triangulation(v)`, so no geometry objects need to be stored outside
the physics module.
"""
materialize(param::Float64, ::Any) = param
materialize(param::Function, trian) = CellField(param, trian)

"""
   _add_contribution(a, b)
   
Utility function to sum contributions from multiple forms, handling `nothing` values.
"""
function _add_contribution(a, b)
    isnothing(a) && return b
    isnothing(b) && return a
    return a + b
end

"""
    _forcing_contribution(f, v, dΩ)

Body-forcing term `∫ v f dΩ`, or `nothing` when `f` is an exact numeric
zero (the default forcing when no `rhs_fn` is given), so that no zero
integral is assembled at every step.
"""
_forcing_contribution(f, v, dΩ) = (f isa Number && iszero(f)) ? nothing : ∫(v * f)dΩ

_space_measure_key(s) = Symbol("d", getfield(s, :space_domain_symbol))
_space_measure(dom::IntegrationDomains, s) = dom[_space_measure_key(s)]

function _space_measure(dom::IntegrationDomains, entities::AbstractVector)
    isempty(entities) && throw(ArgumentError("Cannot resolve a space-domain measure for an empty entity vector."))
    symbols = unique(getfield.(entities, :space_domain_symbol))
    length(symbols) == 1 || throw(ArgumentError("Entity vector has inconsistent `space_domain_symbol` values: $symbols"))
    return dom[Symbol("d", first(symbols))]
end


"""
    _as_space_function(v)

Utility function to convert a value `v` to a function if it is not already one.
This allows for flexible specification of parameters that can be either constants 
or spatially varying functions.
"""
_as_space_function(v) = v isa Function ? v : (x -> v)

"""
    _resolve_space_function(v, ctx)

Resolve a user input into a pure space function using any transient metadata
present in the assembly context.

Supported input shapes:
- constant values
- space functions: `x -> ...`
- time-indexed space functions: `t -> (x -> ...)`
- space-time functions: `(x, t) -> ...`
"""
function _resolve_space_function(v, t)
    !(v isa Function) && return (x -> v)

    if isnothing(t)
        return _as_space_function(v)
    end

    # The shape of `v` is resolved once here, so the returned closure is
    # evaluated at every quadrature point without any try/catch.
    if _has_arity(v, 1)
        # `t -> (x -> ...)` or a space function `x -> ...`: only calling it tells.
        vt = try
            v(t)
        catch
            nothing
        end
        vt isa Function && return vt
    end
    _has_arity(v, 2) && return x -> v(x, t)
    return v
end

# True if `f` has a method taking exactly `n` positional arguments.
_has_arity(f, n::Int) = any(m -> m.nargs == n + 1 && !m.isva, methods(f))

# Frequency-domain: time is not meaningful, so fall back to pure space resolution.
_resolve_space_function(v, ::AC.FrequencyAssemblyContext) = _resolve_space_function(v, nothing)
# Time-domain: extract the current time from the context and resolve.
_resolve_space_function(v, ctx::AC.TimeAssemblyContext) = _resolve_space_function(v, AC.current_time(ctx))
# IntegrationDomains has no time, behave the same as the time=nothing case.
_resolve_space_function(v, ::IntegrationDomains) = _resolve_space_function(v, nothing)
