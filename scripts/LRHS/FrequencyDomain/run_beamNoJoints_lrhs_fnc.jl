using Parameters
using Printf
using Gridap
using WaveSpec
using .Constants
using HydroElasticFEM: print_properties,  Resonator, BeamNoJoints
using HydroElasticFEM: PKG_ROOT

# Here you may include files from the source directory
# include( joinpath(PKG_ROOT,"src","LRHS",
#   "FrequencyDomain","mem_freq_lrmm_free_fnc.jl") )
include( joinpath(PKG_ROOT,"src","LRHS",
  "FrequencyDomain","beamNoJoints_freq_lrhs_rad_fnc.jl") )


resDir::String = "data/sims_202601/runlrhs"

# # Remove all contents in the folder resDir
# if isdir(resDir)
#   rm(resDir; recursive=true, force=true)
# end
# mkpath(resDir)

# Warm-up run
# ---------------------Start--------------------- 
isdir(resDir*"/warmup") || mkpath(resDir*"/warmup")
params = BeamLRHS2D.Beam_LRHS_warmup(name = resDir*"/warmup")
BeamLRHS2D.main(params)
rm(resDir*"/warmup"; recursive=true, force=true)
# ----------------------End----------------------

ρw = params.ρw

# Production run
# ---------------------Start---------------------

H0 = 10
Lb = 2*H0

# Beam parameters
beam2D = BeamNoJoints.Beam2D(
  Lb,   # L
  192.956,   # m
  500e6,    # E
  6.667e-4,    # I 
  0.0,    # τ
  BeamNoJoints.Free()
)

@with_kw struct run_params
  name = resDir
  order::Int = 4
  vtk_output::Bool = true

  # Constants
  ρw = ρw

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
  ω = [0.6:0.1:3.5;]
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)


  # Beam parameters
  beam2D = beam2D
  
  # Domain 
  nx = 1300
  ny = 20
  mesh_ry = 1.15 #Ratio for Geometric progression of eleSize
  Ld = 0*H0 #damping zone length
  LΩ = 2*30*H0+beam2D.L #2*Ld
  x0 = -30*H0
  domain =  (x0, x0+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdin = 0.0
  xb0 = xdin 
  xb1 = xb0 + beam2D.L  

  # Resonator parameters
  resn = Resonator.Array1D(
    1, 
    [0.1*beam2D.MTotal], 
    [0.1*beam2D.MTotal*(1.0^2)], 
    [0.0],
    [Point(xb0 + beam2D.L/2.0,0.0)]   
  )

  # Probes
  prbx=[  -300.0, -250.0, -200.0, -150.0, -100.0,
    -60, -57, -55, 
    -50, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0,
    150, 200, 250, 300 ]
  prbPowx=[ -50, 100 ]

end
params = run_params()

beamName = "beam_m=" * @sprintf("%0.2f", beam2D.m) *
  "_EI=" * @sprintf("%0.2f", beam2D.EI)

resonatorName = "resnM=" * @sprintf("%0.2f", params.resn[1].M) *
  "_resnK=" * @sprintf("%0.2f", params.resn[1].K)

BeamLRHS2D.main(params)
# ----------------------End----------------------
nothing

# ============================================================================
# Resonator Options
# ============================================================================

# # Resonator parameters
# resn = Resonator.Array1D(
#   1, 
#   [1*ρw], 
#   [2.4*2.4*ρw], 
#   [0.0],
#   [Point(xb0 + beam2D.L/2.0,0.0)]   
# )

# # Resonator parameters
# resn = Resonator.Array1D(
#   2, 
#   [1e3, 1e3], 
#   [2.4*2.4*1e3, 3.4*3.4*1e3], 
#   [0.0, 0.0],
#   # [Point(xb0 + beam2D.L/2.0,0.0)]
#   [Point(90.0, 0.0), Point(85.0, 0.0)]
# )

# # Resonator parameters
# resn = Resonator.Array1D(
#   1, 
#   1.0, 
#   0.0, 
#   0.0,
#   [Point(xb0 + beam2D.L/2.0,0.0)]
# )