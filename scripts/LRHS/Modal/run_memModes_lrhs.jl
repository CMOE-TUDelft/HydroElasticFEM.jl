using Parameters
using Gridap
using Printf
using WaveSpec.Constants
using HydroElasticFEM: Resonator, Membrane
using HydroElasticFEM: PKG_ROOT

include(joinpath(PKG_ROOT,
  "src","LRHS","Modal","config_parameters.jl"))

memBndType = "free" # "free" or "fixed"
analysisType = "complexMass" 
include(joinpath(PKG_ROOT,
  "src","LRHS","Modal","memModes_lrhs_complexMass.jl"))

# memBndType = "free" # "free" or "fixed"
# analysisType = "dampedSys" 
# include(joinpath(PKG_ROOT,
#   "src","LRHS","Modal","memModes_lrhs_dampedSys.jl"))


# Directory for results
resDir::String = "data/sims_202601/lrhs_modes_free/"


# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.05, 0.10, 0.25, 0.50, 0.75, 0.8]

# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.8]

m_rho = [0.9] 
T_rho = [98.1] 

H0 = 10
ρw = Memb_LRHS_params().ρw
Lm = 2*H0


# Common parameters for all runs
paramsBase = Memb_LRHS_params(

  vtk_output = false,

  H0 = H0, #m #still-water depth

  # nx = 60,
  # ny = 5,
  
  # memb2D = memb2D
  # resn_ρw = rS,

  nωₙ = 5, #number of natural frequencies to compute

  # Iterative solution for wet natural frequencies
  αRelax = 0.5,
  maxIter = 20,

  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, 1) .+ 1.0
)


# Loop over all combinations of m_rho, T_rho and resonator parameters
for imrho in m_rho
  for iTrho in T_rho

    memb2D = Membrane.Membrane2D(
      2*H0,         # L
      imrho*ρw,      # m
      iTrho*ρw,  # T
      0.0,         # τ
      Membrane.Free()  # bndType
    )

    membName = "memb_mrho=" * @sprintf("%0.2f", memb2D.m/ρw) *
      "_Trho=" * @sprintf("%0.2f", memb2D.T/ρw)

    caseDir = resDir*membName
    isdir(caseDir) || mkpath(caseDir)
    # if( isdir(caseDir) )
    #   rm(caseDir, recursive=true) #remove old data
    #   @printf("Removed old data in %s\n", caseDir)
    #   # return
    # end
    # mkdir(caseDir)    

    iresnMᵨ = 0.1*memb2D.MTotal
    iresnKᵨ = iresnMᵨ*1.0
    resn = Resonator.Single( iresnMᵨ, iresnKᵨ, 0.0, Point(30.0, 0.0) )

    resonatorName = "resnM=" * @sprintf("%0.2f", iresnMᵨ) *
      "_resnK=" * @sprintf("%0.2f", iresnKᵨ)
    

    # Update paramsBase for each run
    params = Memb_LRHS_params(
      paramsBase;

      resDir = caseDir,
      fileName = "lrhs_"*analysisType*"_"*memBndType*"_"*resonatorName,

      memb2D = memb2D,
      resn = resn
    )

    # Run case
    MembraneModes.run_case(params)
  end
end