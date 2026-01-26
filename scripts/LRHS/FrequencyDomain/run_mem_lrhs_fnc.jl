using Parameters
using Gridap
using WaveSpec
using .Constants
using HydroElasticFEM.Resonator
using HydroElasticFEM: PKG_ROOT

# Here you may include files from the source directory
# include( joinpath(PKG_ROOT,"src","LRHS",
#   "FrequencyDomain","mem_freq_lrmm_free_fnc.jl") )
include( joinpath(PKG_ROOT,"src","LRHS",
  "FrequencyDomain","mem_freq_lrhs_rad_fnc.jl") )

ρw = 1025 #kg/m3 water

resDir::String = "data/sims_202601/runlrmm"

# Remove all contents in the folder resDir
if isdir(resDir)
  rm(resDir; recursive=true, force=true)
end
mkpath(resDir)

# Warm-up run
params = Memb2D.Memb_params_warmup(name = resDir)
Memb2D.main(params)

# Remove all contents in the folder resDir
if isdir(resDir)
  rm(resDir; recursive=true, force=true)
end
mkpath(resDir)

# Production run
@with_kw struct run_params
  name = resDir
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
  ω = [1:0.1:5;]
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)
  # k = dispersionRelAng.(H0, ω; msg=false)

  # Membrane parameters
  memBndType="free"  #"free" or "fixed"
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
    [1*ρw], 
    [2.4*2.4*ρw], 
    [0.0],
    [Point(xm₀ + Lm/2.0,0.0)]    
  )
  
  # # Resonator parameters
  # rS = Resonator.Array1D(
  #   2, 
  #   [1e3, 1e3], 
  #   [2.4*2.4*1e3, 3.4*3.4*1e3], 
  #   [0.0, 0.0],
  #   # [Point(xm₀ + Lm/2.0,0.0)]
  #   [Point(90.0, 0.0), Point(85.0, 0.0)]
  # )

  # # Resonator parameters
  # rS = Resonator.Array1D(
  #   1, 
  #   1.0, 
  #   0.0, 
  #   0.0,
  #   [Point(xm₀ + Lm/2.0,0.0)]
  # )
  

  # Probes
  prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
          52.7, 53.7, 55, 60.0, 80.0, 
          85.0, 90.0, 95.0, 100.0, 120.0, 
          125.0, 140.0, 160.0, 180.0 ]
  prbPowx=[ 55.0, 125.0 ]

end
params = run_params()
Memb2D.main(params)
