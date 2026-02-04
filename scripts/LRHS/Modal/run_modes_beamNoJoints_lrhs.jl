using Parameters
using Gridap
using Printf
using WaveSpec.Constants
using HydroElasticFEM: print_properties,  Resonator, BeamNoJoints
using HydroElasticFEM: PKG_ROOT


bndType = "free" # "free" or "fixed"
analysisType = "complexMass" 
include(joinpath(PKG_ROOT,
  "src","LRHS","Modal","modes_beamNoJoints_lrhs_complexMass.jl"))


# Directory for results
resDir::String = "data/beam_modes_free/"


H0 = 25
ρw = BeamModes.Beam_LRHS_params().ρw

# Define beam parameters
# ---------------------Start---------------------
beam2D = BeamNoJoints.Beam2D()
let
  global beam2D
  E = 250e6 #Pa
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
x0 = 0.0
xb0 = 0.0 + beam2D.L
xb1 = xb0 + beam2D.L
LΩ = beam2D.L + beam2D.L + beam2D.L
nx = 120
ny = 10
mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
# ----------------------End----------------------


# Define resonator parameters
# ---------------------Start---------------------
resn = Resonator.Single()
let
  global resn
  resnM = 0.05 * beam2D.MTotal
  ωn = 0.50 #rad/s
  resnK = resnM * ωn^2
  resnC = 0.05 * resnK
  resn = Resonator.Single(
    resnM, resnK, resnC,
    Point(xb0 + beam2D.L/2.0, 0.0)
  )
end

resonatorName = "resnM=" * @sprintf("%0.2f", resn.M) *
  "_resnK=" * @sprintf("%0.2f", resn.K)

# ----------------------End----------------------

# Sanity check
@printf(
  "[IMP] Total mass of the structure: %0.2f kg per unit width \n", 
  beam2D.MTotal + resn.M )

# Common parameters for all runs
paramsBase = BeamModes.Beam_LRHS_params(

  order = 4,
  vtk_output = true,

  H0 = H0, #m #still-water depth
  
  nωₙ = 6, #number of natural frequencies to compute

  # Iterative solution for wet natural frequencies
  αRelax = 0.5,
  maxIter = 20  
)


caseDir = resDir*beamName
isdir(caseDir) || mkpath(caseDir)
# if( isdir(caseDir) )
#   rm(caseDir, recursive=true) #remove old data
#   @printf("Removed old data in %s\n", caseDir)
#   # return
# end
# mkdir(caseDir)    


# Update paramsBase for each run
params = BeamModes.Beam_LRHS_params(
  paramsBase;

  resDir = caseDir,
  fileName = "lrhs_"*analysisType*"_"*bndType*"_"*resonatorName,

  beam2D = beam2D,  
  resn = resn,

  # Domain 
  nx = nx,
  ny = ny,
  mesh_ry = mesh_ry, #Ratio for Geometric progression of eleSize
  LΩ = LΩ, 
  x0 = x0,
  xb0 = xb0,
  xb1 = xb1  
)

# Run case
BeamModes.run_case(params)