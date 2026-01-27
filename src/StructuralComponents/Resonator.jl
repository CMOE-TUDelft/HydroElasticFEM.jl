module Resonator

using Parameters
using Gridap
using Gridap.CellData


"""
Custom Structs
=============

"""
# ---------------------Start---------------------
@with_kw struct Single
	
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
# ----------------------End----------------------



"""
Functions
=============

"""
# ----------------------Start--------------------
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

  rA = []

  for (m,k,c,xz) in zip(M,K,C,XZ)
    push!(rA, Single(m,k,c,xz))
  end

  return rA
end
# ----------------------End----------------------

end