"""
Memb_LRHS_params

Parameters for the Memb2D LRHS modal analysis.
"""
@with_kw struct Memb_LRHS_params

  resDir::String = "data/sims_202508/mem_modes_free/"
  fileName::String = "lrhs"

  # Constants
  ρw = 1025 #kg/m3 water

  order::Int = 2
  vtk_output::Bool = true  

  H0 = 10 #m #still-water depth

  # Membrane parameters
  memb2D = Membrane.Membrane2D(
    2*H0,   # L
    0.9*ρw,   # m
    0.1/4*g* (2*H0)^2 *ρw,    # T 
    0.0,    # τ
    Membrane.Free()
  )

  # Domain 
  nx = 120
  ny = 10
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  LΩ = 6*H0 
  x₀ = 0.0
  xm₀ = x₀ + 2*H0
  xm₁ = xm₀ + memb2D.L  

  # Resonator parameters divided by ρw
  resn = Resonator.Single( 1000, 0.0, 0.0, Point(3*H0, 0.0) )

  # Number of natural frequencies
  nωₙ = 6

  # Iterative solution for wet natural frequencies
  αRelax = 0.8
  maxIter = 20
  
  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, nωₙ) .+ 1.0
  
  
end