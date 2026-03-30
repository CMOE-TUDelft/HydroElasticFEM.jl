
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