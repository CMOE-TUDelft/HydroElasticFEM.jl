using Revise
using Parameters
using Printf
using Gridap
using WaveSpec
using .Constants
using HydroElasticFEM: print_properties, Resonator, Membrane, print_properties
using HydroElasticFEM: PKG_ROOT

# Here you may include files from the source directory
# include( joinpath(PKG_ROOT,"src","LRHS",
#   "FrequencyDomain","mem_freq_lrmm_free_fnc.jl") )
include( joinpath(PKG_ROOT,"src","LRHS",
  "FrequencyDomain","mem_freq_lrhs_rad_fnc.jl") )

include( joinpath(PKG_ROOT,"Paper2026","scripts","LRHS",
  "FrequencyDomain","_plotting_utilities.jl") )

resDir::String = "data/paper2026_01/runlrhs"


# Define the test case
# ---------------------Start---------------------
ρw = MembLRHS2D.Memb_LRHS_warmup().ρw
H0 = 25 #m #still-water depth

# Membrane parameters
memb2D = Membrane.Membrane2D(Membrane.Free())
let
  global memb2D
  
  E = 100e6 #Pa
  ρ = 910 #kg/m3

  Lm = 50.0 #m
  h_memb = 0.05*16 #m
  ϵx = 0.002 #Pre Strain 

  σx = E*ϵx #Pre stress
  T = σx*h_memb #Pre tension

  m = ρ*h_memb #Mass per unit area

  @printf("Pre-Stress σx = %.2f MPa\n", σx/1e6)
  
  memb2D = Membrane.Membrane2D(
    Lm,         # L
    m,      # m
    T,  # T
    0.0,         # τ
    Membrane.Free()  # bndType
  )
end

# Position of the membrane in the domain
xm0 = 0.0

# Resonator parameters
resn = Resonator.Array1D(
  1, 
  [0.10*memb2D.MTotal], 
  [0.10*memb2D.MTotal*(0.5^2)], 
  [0.0],
  [Point(xm0 + memb2D.L/2.0,0.0)]   
)

print_properties(memb2D)
print_properties(resn)

membName = "memb_mrho=" * @sprintf("%0.2f", memb2D.m/ρw) *
  "_Trho=" * @sprintf("%0.2f", memb2D.T/ρw)

resonatorName = "resnM=" * @sprintf("%0.2f", resn[1].M) *
  "_resnK=" * @sprintf("%0.2f", resn[1].K)

# ----------------------End----------------------


# # Remove all contents in the folder resDir
# if isdir(resDir)
#   rm(resDir; recursive=true, force=true)
# end
# mkpath(resDir)

# Warm-up run
# ---------------------Start--------------------- 
isdir(resDir*"/warmup") || mkpath(resDir*"/warmup")
params = MembLRHS2D.Memb_LRHS_warmup(name = resDir*"/warmup")
MembLRHS2D.main(params)
rm(resDir*"/warmup"; recursive=true, force=true)
# ----------------------End----------------------



# Production run
# ---------------------Start---------------------
# ω_run = [0.6:0.02:4;]
ω_run = [0.3:0.02:2.5;]
ω_nat = [
  0.51797, 0.931481, 1.87729, #2.79444, 3.72934, 4.65741,
  0.497066, 0.960528, 1.31806, 1.71873, 2.17691
]

push!(ω_run, ω_nat...)
ω_run = unique(sort(ω_run))

@with_kw struct run_params
  name = resDir
  order::Int = 2
  vtk_output::Bool = true

  ρw = ρw

  H0 = H0 #m #still-water depth

  # Wave parameters
  # ω, S, η₀ = jonswap(0.4, 2.5; 
  #     plotflag=true, plotloc=filename, nω=145)
  # println(ω[1], "\t", ω[2], "\t", ω[end])
  # ω = ω[2:end]
  # S = S[2:end]
  # η₀ = η₀[2:end]
  # ω = [0.6:0.1:3.5;]
  ω = ω_run
  T = 2*π./ω
  η₀ = 0.10*ones(length(ω))
  α = randomPhase(ω; seed=100)


  # Membrane parameters
  memb2D = memb2D
  
  # Domain 
  nx = 1300
  ny = 20
  mesh_ry = 1.15 #Ratio for Geometric progression of eleSize
  Ld = 0*H0 #damping zone length
  LΩ = 2*300+memb2D.L #2*Ld
  x₀ = -300
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)
  xdᵢₙ = 0.0
  xm₀ = xm0
  xm₁ = xm0 + memb2D.L  
  
  # Resonator parameters
  resn = resn  

  # Probes
  prbx=[  -300.0, -250.0, -200.0, -150.0, -100.0,
    -60, -57, -55, 
    -50, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 100.0,
    150, 200, 250, 300 ]
  prbPowx=[ -50, 100 ]

end
params = run_params()

cache_for_plots = MembLRHS2D.main(params)
cache_for_plots = (; cache_for_plots..., 
  ω_nat
)
@show cache_for_plots

# plot_power_coefficient(cache_for_plots)
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
#   [Point(xm₀ + memb2D.L/2.0,0.0)]    
# )

# # Resonator parameters
# resn = Resonator.Array1D(
#   2, 
#   [1e3, 1e3], 
#   [2.4*2.4*1e3, 3.4*3.4*1e3], 
#   [0.0, 0.0],
#   # [Point(xm₀ + memb2D.L/2.0,0.0)]
#   [Point(90.0, 0.0), Point(85.0, 0.0)]
# )

# # Resonator parameters
# resn = Resonator.Array1D(
#   1, 
#   1.0, 
#   0.0, 
#   0.0,
#   [Point(xm₀ + memb2D.L/2.0,0.0)]
# )