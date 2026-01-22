using Parameters
using Gridap
using Printf
using WaveSpec.Constants
# using HydroElasticFEM.Resonator
using HydroElasticFEM: PKG_ROOT


memBndType = "free" # "free" or "fixed"
include(joinpath(PKG_ROOT,
  "src","Membrane2D","Modal","memModes_complexMass.jl"))

# caseTypeName = "memb_free"
# include(joinpath(PKG_ROOT,
#   "src","Membrane2D","Modal","memModes_dampedSys_free.jl"))


# Directory for results
resDir::String = "data/sims_202512/mem_modes_free/"


# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.05, 0.10, 0.25, 0.50, 0.75, 0.8]

# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.8]

mfac = [0.9] 
tfac = [0.1] 

H0 = 10


# Common parameters for all runs
paramsBase = MembraneModes.Memb_params(

  vtk_output = false,

  memBndType = memBndType,

  H0 = H0, #m #still-water depth
  Lm = 2*H0, #m

  # nx = 60,
  # ny = 5,
  
  # mfac = imfac,
  # tfac = itfac,
  # rS_by_ρw = rS,

  nωₙ = 4, #number of natural frequencies to compute

  # Iterative solution for wet natural frequencies
  αRelax = 0.5,
  maxIter = 20,

  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, 1) .+ 1.0
)


# Loop over all combinations of mfac, tfac and resonator parameters
for imfac in mfac
  for itfac in tfac

    membName = "mfac=" * @sprintf("%0.2f", imfac) *
      "_tfac=" * @sprintf("%0.2f", itfac)

    # Case directory
    # caseName = "mem_modes_ten" * @sprintf("%0.2f", itfac) *
    #   "_mass" * @sprintf("%0.2f", imfac)
    caseName = "run"

    caseDir = resDir*caseName
    # if( isdir(caseDir) )
    #   rm(caseDir, recursive=true) #remove old data
    #   @printf("Removed old data in %s\n", caseDir)
    #   # return
    # end
    # mkdir(caseDir)
    

    # Update paramsBase for each run
    params = MembraneModes.Memb_params(
      paramsBase;

      resDir = caseDir,
      fileName = "mem_"*memBndType*"_"*membName,

      mfac = imfac,
      tfac = itfac,      
    )

    # Run case
    MembraneModes.run_case(params)
  end
end