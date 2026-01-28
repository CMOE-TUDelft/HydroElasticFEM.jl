module Membrane

using Parameters
using Gridap
using Gridap.CellData
using Printf


"""
Custom Structs
=============

"""
# ---------------------Start---------------------
abstract type MembraneBndType end
struct Free <: MembraneBndType end
struct Fixed <: MembraneBndType end


@with_kw struct Membrane2D
	
  L::Real # Length of membrane
  m::Real # Mass per unit length per unit width
  T::Real # Pre-Tension
  τ::Real # Proportional Structural Damping coefficient
  bndType::MembraneBndType # Boundary Type (:free or :fixed)

  # Derived quantities
  MTotal::Real # Total Mass per unit width
  ωn1::Real # Dry Analytical Natural frequency

  function Membrane2D( L, m, T, τ, bndType::MembraneBndType )
    MTotal = m * L
    ωn1 = (π / L) * sqrt( T / m )
    new( L, m, T, τ, bndType, MTotal, ωn1 )
  end

end

function Membrane2D(bndType::MembraneBndType)
  Membrane2D( 0.0, 0.0, 0.0, 0.0, bndType )
end

# ----------------------End----------------------



"""
Functions
=============

"""
# ----------------------Start--------------------
function print_membrane_props( memb2D::Membrane2D, ρw::Real = 1025 )
  
  mᵨ = memb2D.m / ρw
  Tᵨ = memb2D.T / ρw

  @printf("\n[MSG] Membrane Properties:\n")
  @printf("[VAL] Density of water, ρw = %.2f kg/m3\n", ρw)
  @printf("[VAL] Lm = %.4f m\n", memb2D.L)
  @printf("[VAL] m, mᵨ = %.4f kg/m2, %.4f m\n", memb2D.m, mᵨ)
  @printf("[VAL] T, Tᵨ = %.4f MN/m, %.4f m3/s2\n", memb2D.T/1e6, Tᵨ)
  @printf("[VAL] τ = %.4f \n", memb2D.τ)
  @printf("[VAL] memBndType = %s \n", string(memb2D.bndType))
  @printf("[VAL] MTotal = %.4f kg/m \n", memb2D.MTotal)
  @printf("[VAL] 1st Dry Analytical Natural Freq, ωn1 = %.4f rad/s \n", memb2D.ωn1)
  println()
  
end
# ----------------------End----------------------

end