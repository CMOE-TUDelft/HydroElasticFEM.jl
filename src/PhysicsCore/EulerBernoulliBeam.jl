"""
EBBeamParameters <: PhysicsParameters

This struct defines the parameters for an Euler-Bernoulli beam model. It includes the following fields:
  - `E: Float64 = 210e9`: Young's modulus of the beam material. 
  - `I: Float64 = 1e-6`: Second moment of area of the beam cross-section.
"""
@with_kw struct EBBeamParameters <: PhysicsParameters
    E::Float64 = 210e9  # Young's modulus in Pascals (default value for steel)
    I::Float64 = 1e-6    # Second moment of area in m^4 (default value for a small beam)
end

"""
print_parameters(params::EBBeamParameters)

Print the parameters of the Euler-Bernoulli beam model in a human-readable format.
"""
function print_parameters(params::EBBeamParameters)
    @printf("\n[MSG] Beam Properties:\n")
    @printf("Young's Modulus (E): %.2e Pa\n", params.E)
    @printf("Second Moment of Area (I): %.2e m^4\n", params.I)
end
