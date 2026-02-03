using HydroElasticFEM: Membrane, BeamNoJoints

# ---------------------Start---------------------
"""
Memb_params_warmup

Parameters for the Memb2D module.
"""
@with_kw struct Memb_LRHS_warmup
  name::String = "data/sims_202306/run/spec_free"
  order::Int = 2
  vtk_output::Bool = true

  # Constants
  ρw = 1025 #kg/m3 water

  H0 = 10 #m #still-water depth

  # Wave parameters
  # ω, S, η₀ = jonswap(0.4, 2.5; 
  #     plotflag=true, plotloc=filename, nω=145)
  # println(ω[1], "\t", ω[2], "\t", ω[end])
  # ω = ω[2:end]
  # S = S[2:end]
  # η₀ = η₀[2:end]
  # ω = [2*π/2.53079486745378, 2*π/2.0]
  # η₀ = [0.25, 0.25]
  ω = 1:1:2
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)

  # Membrane parameters
  memb2D = Membrane.Membrane2D(
    1*H0,         # L
    0.9*ρw,      # m
    0.1/4*g* (1*H0)^2 *ρw,  # T
    0.0,         # τ
    Membrane.Free()  # bndType
  )  

  # Domain 
  nx = 50
  ny = 8
  mesh_ry = 1.05 #Ratio for Geometric progression of eleSize
  Ld = 0*H0 #damping zone length
  LΩ = 5*H0 + Ld #2*Ld
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xdᵢₙ + 2*H0
  xm₁ = xm₀ + memb2D.L


  # Resonator parameters
  resn = Resonator.Array1D(
    1, 
    1e3, 
    5.9e3, 
    0.0,
    [Point(xm₀ + memb2D.L/2.0,0.0)]
  )


  # Probes
  prbx=[  0.0, 0.0, 2.0, 5.0, 10.0, 
          12.7, 13.7, 15, 18.0, 20.0, 
          22.0, 25.0, 28.0, 30.0, 35.0, 
          40.0, 45.0, 50.0, 50.0 ]
  prbPowx=[ 15.0, 35.0 ]

end
# ----------------------End----------------------



# ---------------------Start---------------------
"""
Beam_LRHS_warmup

Parameters for the Beam2D module.
"""
@with_kw struct Beam_LRHS_warmup
  name::String = "data/sims_202306/run/spec_free"
  order::Int = 4
  vtk_output::Bool = true

  # Constants
  ρw = 1025 #kg/m3 water

  H0 = 10 #m #still-water depth

  # Wave parameters
  # ω, S, η₀ = jonswap(0.4, 2.5; 
  #     plotflag=true, plotloc=filename, nω=145)
  # println(ω[1], "\t", ω[2], "\t", ω[end])
  # ω = ω[2:end]
  # S = S[2:end]
  # η₀ = η₀[2:end]
  # ω = [2*π/2.53079486745378, 2*π/2.0]
  # η₀ = [0.25, 0.25]
  ω = 1:1:2
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)

  # Beam parameters
  beam2D = BeamNoJoints.Beam2D(
    H0,   # L
    192.956,   # m
    500e6,    # E
    6.667e-4,    # I 
    0.0,    # τ
    BeamNoJoints.Free()
  )

  # Domain 
  nx = 25
  ny = 4
  mesh_ry = 1.05 #Ratio for Geometric progression of eleSize
  Ld = 0*H0 #damping zone length
  LΩ = 5*H0 + Ld #2*Ld
  x0 = -Ld
  domain =  (x0, x0+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdin = 0.0
  xb0 = xdin + 2*H0
  xb1 = xb0 + beam2D.L


  # Resonator parameters
  resn = Resonator.Array1D(
    1, 
    1e3, 
    5.9e3, 
    0.0,
    [Point(xb0 + beam2D.L/2.0,0.0)]
  )


  # Probes
  prbx=[  0.0, 0.0, 2.0, 5.0, 10.0, 
          12.7, 13.7, 15, 18.0, 20.0, 
          22.0, 25.0, 28.0, 30.0, 35.0, 
          40.0, 45.0, 50.0, 50.0 ]
  prbPowx=[ 15.0, 35.0 ]

end
# ----------------------End----------------------