module BeamNoJoints

using Gridap
using Gridap.CellData
using Printf


"""
Custom Structs
=============

"""
# ---------------------Start---------------------
abstract type BeamBndType end
struct Free <: BeamBndType end
struct Fixed <: BeamBndType end


struct Beam2D
	
  L::Real # Length of membrane
  m::Real # Mass per unit length per unit width
  E::Real # Young's Modulus
  I::Real # Second Moment of Area 
  τ::Real # Stiffness Proportional Structural Damping coefficient 
  bndType::BeamBndType # Boundary Type (:free or :fixed)
  
  # Derived quantities
  EI::Real # Flexural Rigidity
  τEI::Real # Damping Rigidity
  MTotal::Real # Total Mass per unit width
  ωn1::Real # Dry Analytical Natural frequency

  function Beam2D( L, m, E, I, τ, bndType::BeamBndType )
    MTotal = m * L
    EI = E * I
    τEI = τ * EI
    
    ωn1 = 22.3733 * sqrt( EI / (m * L^4) )    
    """
    See free-free vibration frequency formula from

    [Wiki](https://en.wikipedia.org/wiki/Euler%E2%80%93Bernoulli_beam_theory)
    """

    new( L, m, E, I, τ, bndType, EI, τEI, MTotal, ωn1 )
  end   

end

function Beam2D(bndType::BeamBndType = Free())
  Beam2D( 0.0, 0.0, 0.0, 0.0, 0.0, bndType )
end



# ----------------------End----------------------



"""
Functions
=============

"""
# ----------------------Start--------------------
function print_properties( beam2D::Beam2D, ρw::Real = 1025 )
  
  mᵨ = beam2D.m / ρw
  EIᵨ = beam2D.EI / ρw

  @printf("\n[MSG] Beam Properties:\n")
  @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", ρw)
  @printf("[VAL] L = %.4f m\n", beam2D.L)
  @printf("[VAL] m, mᵨ = %.4f kg/m2, %.4f m\n", beam2D.m, mᵨ)
  @printf("[VAL] E = %.4f Pa\n", beam2D.E)
  @printf("[VAL] I = %.4f m4/m\n", beam2D.I)
  @printf("[VAL] τ = %.4f \n", beam2D.τ)
  @printf("[VAL] EI, EIᵨ = %.4f Nm2/m, %.4f m5/s2\n", beam2D.EI, EIᵨ)
  @printf("[VAL] τEI = %.4f \n", beam2D.τEI)
  @printf("[VAL] beamBndType = %s \n", string(beam2D.bndType))
  @printf("[VAL] MTotal = %.4f kg \n", beam2D.MTotal)
  @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s\n", beam2D.ωn1)
  @printf("[MSG] See free-free vibration frequency formula involving 22.3733 from Wiki.\n")
  println()
  
end
# ----------------------End----------------------

end 