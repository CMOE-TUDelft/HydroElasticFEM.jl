module Membrane

using Parameters
using Gridap
using Gridap.CellData


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
# ----------------------End----------------------



"""
Functions
=============

"""
# ----------------------Start--------------------

# ----------------------End----------------------

end