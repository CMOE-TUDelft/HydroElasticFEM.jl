# --------------------Start--------------------
"""
  map_vertical_GP_for_const_dep(y, r, n, H0; dbgmsg = false)

Scale vertical coordinates using a geometric progression for constant depth meshes.

# Arguments
- `y`: Vertical coordinate to transform (expected range: -H0 to 0)
- `r::Real`: Geometric ratio for mesh spacing. If r ≈ 1.0, uniform spacing is used
- `n::Integer`: Number of mesh divisions along depth
- `H0::Real`: Total depth (positive value, mesh extends from -H0 to 0)
- `dbgmsg`: Optional debug flag to print mesh spacing information (default: false)

# Returns
- Transformed vertical coordinate based on geometric progression

# Description
This function transforms vertical coordinates from a uniform distribution to a 
geometrically progressed distribution. The mesh is scaled along the depth using 
a geometric progression (GP) where the depth ranges from -H0 to 0.

# Special Cases
- If r ≈ 1.0, returns y (uniform mesh)
"""
function map_vertical_GP_for_const_dep(
  y, r::Real, n::Integer, H0::Real; dbgmsg = false )    

  # Mesh along depth as a GP
  # Depth is 0 to -H0    
  if(r ≈ 1.0)
    return y  
  else
    a0 = H0 * (r-1) / (r^n - 1)    
    if(dbgmsg)
      ln = 0:n
      ly = -a0 / (r-1) * (r.^ln .- 1)         
      @show hcat( ly, [ 0; ly[1:end-1] - ly[2:end] ] )
    end
    
    if y ≈ 0
      return 0.0
    end
    j = abs(y) / H0 * n  
    return -a0 / (r-1) * (r^j - 1)
  end
end
# ----------------------End---------------------


# --------------------Start--------------------
"""
  print_properties()

Print properties of various objects in a formatted manner.

This is a generic function that serves as a placeholder for method dispatch.
Specific implementations should be defined for different types to display
their relevant properties and parameters in a human-readable format.
"""
function print_properties() end

function print_properties(ele::BeamNoJoints.Beam2D)
    BeamNoJoints.print_properties(ele)
end

function print_properties(ele::Membrane.Membrane2D)
    Membrane.print_properties(ele)
end

function print_properties(ele::Union{Resonator.Single, Vector{Resonator.Single}})
    Resonator.print_properties(ele)
end
# ----------------------End---------------------

