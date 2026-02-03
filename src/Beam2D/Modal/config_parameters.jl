"""
Beam2D Modal Analysis Configuration Parameters.

Parameters for the Beam2D LRHS modal analysis.
"""
@with_kw struct Beam_params

  resDir::String = "data/sims_202508/beam_modes_free/"
  fileName::String = "lrhs"

  # Constants
  ρw = 1025 #kg/m3 water

  order::Int = 2
  vtk_output::Bool = true  

  H0 = 10 #m #still-water depth

  # Beam2D parameters
  beam2D = BeamNoJoints.Beam2D(
    2*H0,   # L
    192.956,   # m
    500e6,    # E
    6.667e-4,    # I 
    0.0,    # τ
    BeamNoJoints.Free()
  )

  # Domain 
  nx = 120
  ny = 10
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  LΩ = 6*H0 
  x₀ = 0.0
  xm₀ = x₀ + 2*H0
  xm₁ = xm₀ + beam2D.L  


  # Number of natural frequencies
  nωₙ = 6

  # Iterative solution for wet natural frequencies
  αRelax = 0.8
  maxIter = 20
  
  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, nωₙ) .+ 1.0
  
  
end