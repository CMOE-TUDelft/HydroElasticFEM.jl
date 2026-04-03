
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

    try
        vt = v(t)
        if vt isa Function
            return vt
        end
    catch
    end

    return x -> begin
        try
            v(x, t)
        catch
            v(x)
        end
    end
end

_resolve_space_function(v, ::AC.FrequencyAssemblyContext) = _resolve_space_function(v, nothing)
_resolve_space_function(v, ctx::AC.TimeAssemblyContext) = _resolve_space_function(v, AC.current_time(ctx))
_resolve_space_function(v, ::IntegrationDomains) = _resolve_space_function(v, nothing)
