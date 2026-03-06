module Resonator

using Parameters
using Gridap
using Gridap.CellData
using Printf


"""
Custom Structs
=============

"""
# ---------------------Start---------------------
struct Single
	
	M::Real # Mass
  K::Real # Stiffness
  C::Real # Damping
  XZ::VectorValue{2, Float64} # Position
  ωn1::Real # Natural frequencies
  
  function Single( M::Real, K::Real, C::Real,
    XZ::VectorValue{2, Float64} )
    
    ωn1 = sqrt(K / M)
    new( M, K, C, XZ, ωn1 )
  end

end

function Single()
  Single( 0.0, 0.0, 0.0, VectorValue(0.0, 0.0) )
end
# ----------------------End----------------------



"""
Functions
=============

"""
# ----------------------Start--------------------

function print_properties( resn::Single )
  
  @printf("\n[MSG] Resonator Properties:\n")
  @printf("[VAL] M = %.4f kg\n", resn.M)
  @printf("[VAL] K = %.4f N/m\n", resn.K)
  @printf("[VAL] C = %.4f Ns/m\n", resn.C)
  @printf("[VAL] XZ = (%.4f, %.4f) m\n", resn.XZ[1], resn.XZ[2])
  @printf("[VAL] ωn1 = %.4f rad/s\n", resn.ωn1)
  println()

end

function print_properties( resn::Vector{Single} )
  print_properties.(resn)
end

# @with_kw struct Array1D
	
#   N::Int
# 	M::Vector{<:Real} # Mass
#   K::Vector{<:Real} # Stiffness
#   C::Vector{<:Real} # Damping
#   X::Vector{<:Real} # Position
#   ωn1::Vector{<:Real} # Natural frequencies
  
# end


# function Array1D( 
#   N::Int, M::Vector{<:Real}, K::Vector{<:Real}, 
#   C::Vector{<:Real}, 
#   XZ::Vector{VectorValue{2,Float64}} )

#   if (length(M) != N || length(K) != N || 
#     length(C) != N || length(XZ) != N)

#     throw(ArgumentError("M, K, C, and XZ must be of length N"))
#   end

#   Array1D(N, M, K, C, XZ)
# end


function Array1D( 
  N::Int, M::Real, K::Real, C::Real, 
  XZ::Vector{VectorValue{2,Float64}} )

  if( length(XZ) != N)
    throw(ArgumentError("XZ must be of length N"))
  end

  M = fill(M, N)
  K = fill(K, N)
  C = fill(C, N)

  Array1D(N, M, K, C, XZ)
end


# function Array1D( 
#   N::Int, M::Vector{<:Real}, K::Vector{<:Real}, 
#   C::Vector{<:Real}, X::Vector{<:Real} )

#   ωn1 = sqrt.(K ./ M)

#   Array1D(N, M, K, C, X, ωn1)

# end
function Array1D( 
  N::Int, M::Vector{<:Real}, K::Vector{<:Real}, 
  C::Vector{<:Real}, XZ::Vector{VectorValue{2,Float64}} )

  if (length(M) != N || length(K) != N || 
    length(C) != N || length(XZ) != N)

    throw(ArgumentError("M, K, C, and XZ must be of length N"))
  end

  rA = Vector{Single}()

  for (m,k,c,xz) in zip(M,K,C,XZ)
    push!(rA, Single(m,k,c,xz))
  end

  return rA
end
# ----------------------End----------------------

end