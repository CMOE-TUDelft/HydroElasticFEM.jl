using Parameters
using Gridap
using Printf
using WaveSpec.Constants
using HydroElasticFEM: print_properties,  Resonator, Membrane
using HydroElasticFEM: PKG_ROOT

include(joinpath(PKG_ROOT,
  "src","LRHS","Modal","_config_parameters.jl"))

memBndType = "free" # "free" or "fixed"
analysisType = "complexMass" 
include(joinpath(PKG_ROOT,
  "src","LRHS","Modal","modes_mem_lrhs_complexMass.jl"))

# memBndType = "free" # "free" or "fixed"
# analysisType = "dampedSys" 
# include(joinpath(PKG_ROOT,
#   "src","LRHS","Modal","modes_mem_lrhs_dampedSys.jl"))


# Directory for results
resDir::String = "data/paper2026_01/modes_LRHS/"


# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.05, 0.10, 0.25, 0.50, 0.75, 0.8]

# mfac = [0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0 ]
# tfac = [0.8]

ρw = Memb_LRHS_params().ρw

memb_m = [728] 
memb_T = [160000] 
Lm = 50 #m

H0 = 25.0 #m


# Common parameters for all runs
paramsBase = Memb_LRHS_params(

  vtk_output = true,

  H0 = H0, #m #still-water depth

  # # Domain Mesh1
  # nx = 120,
  # ny = 10,
  # mesh_ry = 1.2, 

  # # Domain Mesh2 
  # nx = 375,
  # ny = 15,
  # mesh_ry = 1.15, 

  LΩ = 3*Lm, 
  x₀ = 0.0,
  xm₀ = Lm,
  xm₁ = 2*Lm,  
  
  # memb2D = memb2D
  # resn_ρw = rS,

  nωₙ = 7, #number of natural frequencies to compute

  # Iterative solution for wet natural frequencies
  αRelax = 0.5,
  maxIter = 20,

  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, 1) .+ 1.0
)


# Loop over all combinations of m_rho, T_rho and resonator parameters
for imemb_m in memb_m
  for imemb_T in memb_T

    memb2D = Membrane.Membrane2D(
      Lm,         # L
      imemb_m,     # m
      imemb_T,  # T
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

    iresnM = 0.10*memb2D.MTotal
    iresnK = 0.0*iresnM*0.5^2 #resonator natural frequency = 0.5 rad/s
    resn = Resonator.Single( iresnM, iresnK, 0.0, Point(1.5*Lm, 0.0) )

    resonatorName = "resnM=" * @sprintf("%0.2f", iresnM) *
      "_resnK=" * @sprintf("%0.2f", iresnK)
    

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