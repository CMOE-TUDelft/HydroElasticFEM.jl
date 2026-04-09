using Parameters
using Printf
using Gridap
using WaveSpec
using .Constants
using HydroElasticFEM: print_properties,  Resonator, BeamNoJoints
using HydroElasticFEM: PKG_ROOT

# Here you may include files from the source directory
include( joinpath(PKG_ROOT,"src","LRHS",
  "FrequencyDomain","beamNoJoints_freq_lrhs_rad_fnc.jl") )

include("./_plotting_utilities.jl")


# Directory for results
resDir::String = "data/runlrhs"

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

H0 = 25
ρw = params.ρw

# Define beam parameters
# ---------------------Start---------------------
beam2D = BeamNoJoints.Beam2D()
let
  global beam2D
  E = 500e6 #Pa
  ρ = 910 #kg/m3

  Lb = 40 #m
  h_outer = 0.25 #m
  h_inner = 0.05 #m 
  I = 1/12*(h_outer^3 - h_inner^3) #m4/m

  m = ρ*(h_outer - h_inner) #Mass per unit length unit width

  beam2D = BeamNoJoints.Beam2D(
    Lb,   # L
    m,   # m
    E,    # E
    I,    # I 
    0.0,    # τ
    BeamNoJoints.Free()
  )
  print_properties(beam2D)

  # Sanity check
  @printf(
    "[IMP] Total possible displaced water mass: %0.2f kg per unit width \n", 
    ρw * beam2D.L * h_outer)
end

beamName = "beam_mrho=" * @sprintf("%0.2f", beam2D.m/ρw) *
  "_EIrho=" * @sprintf("%0.0f", beam2D.EI/ρw)

# Domain positions
x0 = -300
xb0 = 0.0
xb1 = xb0 + beam2D.L
LΩ = 300 + beam2D.L + 300
nx = 1280
ny = 20
mesh_ry = 1.15 #Ratio for Geometric progression of eleSize
# ----------------------End----------------------


# Define resonator parameters
# ---------------------Start---------------------
resn = Resonator.Single()
let
  global resn
  resnM = 0.10 * beam2D.MTotal
  ωn = 0.50 #rad/s
  resnK = resnM * ωn^2
  resnC = 0.0 * resnK
  resn = Resonator.Single(
    resnM, resnK, resnC,
    Point(xb0 + beam2D.L/2.0, 0.0)
  )
end

resonatorName = "resnM=" * @sprintf("%0.2f", resn.M) *
  "_resnK=" * @sprintf("%0.2f", resn.K)

# ----------------------End----------------------

# Run parameters
# ---------------------Start---------------------
ω_run = [0.3:0.02:3.5;]
ω_nat = [
  # Wet LRHS natural frequencies
  0.499384, 0.550596, 1.07248, 1.37029, 1.67608, 2.12124, 2.8505
]
@show ω_nat
push!(ω_run, ω_nat...)
ω_run = unique(sort(ω_run))


run_params = BeamLRHS2D.Beam_LRHS_warmup(
  name = resDir,
  order = 4,
  vtk_output = true,
  
  H0 = H0,

  # Wave Parameters
  ω = ω_run,
  η₀ = 0.10*ones(length(ω_run)),
  α = randomPhase(ω_run; seed=100),

  beam2D = beam2D,
  resn = resn,

  # Domain
  nx = nx,
  ny = ny,
  mesh_ry = mesh_ry, #Ratio for Geometric progression of eleSize
  LΩ = LΩ,
  x0 = x0,
  xb0 = xb0,
  xb1 = xb1,

  # Probes
  prbx=[  -300.0, -250.0, -200.0, -150.0, -100.0,
    -60, -57, -55, 
    -50, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0,
    150, 200, 250, 300 ],
  prbPowx=[ -50, 100 ]

)

cache_for_plots = BeamLRHS2D.main(run_params)

# Add additional variable to cache for plotting
cache_for_plots = (; cache_for_plots..., 
  ω_run, ω_nat )

plot_power_coefficient(cache_for_plots)
plot_resonator_RAO(cache_for_plots)
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