using Parameters
using Gridap
using Printf
using WaveSpec.Constants
using HydroElasticFEM: print_properties,  Resonator, BeamNoJoints
using HydroElasticFEM: PKG_ROOT

include(joinpath(PKG_ROOT,
  "src","Beam2D","Modal","_config_parameters.jl"))

bndType = "free" # "free" or "fixed"
analysisType = "complexMass" 
include(joinpath(PKG_ROOT,
  "src","Beam2D","Modal","modes_beamNoJoints_complexMass.jl"))

# bndType = "free" # "free" or "fixed"
# analysisType = "dampedSys" 
# include(joinpath(PKG_ROOT,
#   "src","Beam2D","Modal","modes_beamNoJoints_dampedSys.jl"))


# Directory for results
resDir::String = "data/sims_202601/beam_modes_free/"


H0 = 10
ρw = Beam_params().ρw
Lb = 2*H0


# Common parameters for all runs
paramsBase = Beam_params(

  vtk_output = true,

  H0 = H0, #m #still-water depth

  # nx = 60,
  # ny = 5,
  
  # beam2D = beam2D

  nωₙ = 6, #number of natural frequencies to compute

  # Iterative solution for wet natural frequencies
  αRelax = 0.5,
  maxIter = 20,

  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, 1) .+ 1.0
)


# Loop over all combinations of m_rho, T_rho and resonator parameters

beam2D = BeamNoJoints.Beam2D(
  Lb,   # L
  192.956,   # m
  500e6,    # E
  6.667e-4,    # I 
  0.0,    # τ
  BeamNoJoints.Free()
)

beamName = "beam_m=" * @sprintf("%0.2f", beam2D.m) *
  "_EI=" * @sprintf("%0.2f", beam2D.EI)

caseDir = resDir*beamName
isdir(caseDir) || mkpath(caseDir)
# if( isdir(caseDir) )
#   rm(caseDir, recursive=true) #remove old data
#   @printf("Removed old data in %s\n", caseDir)
#   # return
# end
# mkdir(caseDir)    


# Update paramsBase for each run
params = Beam_params(
  paramsBase;

  resDir = caseDir,
  fileName = "beam",

  beam2D = beam2D,  
)

# Run case
BeamNoJointsModes.run_case(params)