module Physics

using Parameters

abstract type PhysicsParameters end

function print_parameters(params::PhysicsParameters)
    error("print_parameters not implemented for $(typeof(params))")
end



end