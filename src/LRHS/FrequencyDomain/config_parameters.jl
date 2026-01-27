"""
Memb_params_warmup

Parameters for the Memb2D module.
"""
# ---------------------Start---------------------
@with_kw struct Memb_LRHS_warmup
  name::String = "data/sims_202306/run/spec_free"
  order::Int = 2
  vtk_output::Bool = true

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
  memBndType::String = "free" # "free" or "fixed"
  Lm = 1*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.0#damping coeff


  # Domain 
  nx = 50
  ny = 8
  mesh_ry = 1.05 #Ratio for Geometric progression of eleSize
  Ld = 0*H0 #damping zone length
  LΩ = 50*H0 + Ld #2*Ld
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xdᵢₙ + 2*H0
  xm₁ = xm₀ + Lm


  # Resonator parameters
  rS = Resonator.Array1D(
    1, 
    1e3, 
    5.9e3, 
    0.0,
    [Point(xm₀ + Lm/2.0,0.0)]
  )


  # Probes
  prbx=[  0.0, 0.0, 2.0, 5.0, 10.0, 
          12.7, 13.7, 15, 18.0, 20.0, 
          22.0, 25.0, 28.0, 30.0, 35.0, 
          40.0, 45.0, 50.0, 50.0 ]
  prbPowx=[ 15.0, 35.0 ]

end
# ----------------------End----------------------



"""
Memb_params

Parameters for the Memb2D module.
"""
# ---------------------Start---------------------
@with_kw struct Memb_LRHS_params
  name::String = "data/sims_202306/run/spec_free"
  order::Int = 2
  vtk_output::Bool = true

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
  ω = 0.7:0.5:5
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)

  # Membrane parameters
  memBndType::String = "free" # "free" or "fixed"
  Lm = 2*H0 #m
  Wm = Lm  
  mᵨ = 0.9 #mass per unit area of membrane / ρw
  Tᵨ = 0.1/4*g*Lm*Lm #T/ρw
  τ = 0.0#damping coeff  

  # Domain 
  nx = 1650
  ny = 20
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  Ld = 15*H0 #damping zone length
  LΩ = 18*H0 + Ld #2*Ld
  x₀ = -Ld
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xdᵢₙ + 8*H0
  xm₁ = xm₀ + Lm

  # Resonator parameters
  rS = Resonator.Array1D(
    1, 
    1e3, 
    5.9e3, 
    0.0,
    [Point(xm₀ + Lm/2.0,0.0)]
  )

  # Probes
  prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
          52.7, 53.7, 55, 60.0, 80.0, 
          85.0, 90.0, 95.0, 100.0, 120.0, 
          125.0, 140.0, 160.0, 180.0 ]
  prbPowx=[ 55.0, 125.0 ]

end
# ----------------------End----------------------