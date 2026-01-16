module Membrane_modes

using Parameters
using JLD2
using Gridap
using Plots
using WaveSpec.Constants
using LinearAlgebra
using TickTock
using DataFrames
using Printf


function run_case( params )
    
  @unpack mfac, tfac = params
  @unpack order, vtk_output = params
  @unpack H0, Lm = params
  @unpack nωₙ = params

  @unpack αRelax, maxIter, ωn_guess = params

  @unpack resDir, fileName = params
  
  fileName = resDir*"/"*fileName

  # Constants
  ρw = 1025 #kg/m3 water    
  @show g #defined in .Constants
  
  # Membrane parameters
  mᵨ = mfac  #mass per unit area of membrane / ρw
  Tᵨ = tfac *g*H0*H0 #T/ρw

  # Stabilisation
  βh = 1.0#0.5

  # Domain 
  @unpack nx, ny, mesh_ry, LΩ = params
  @unpack x₀, xm₀, xm₁ = params
  domain =  (x₀, x₀+LΩ, -H0, 0.0)
  partition = (nx, ny)  
  @show Lm
  @show LΩ
  @show domain
  @show partition
  @show mesh_ry
  @show (xm₀, xm₁)
  @show isinteger(Lm/LΩ*nx)
  @show LΩ/nx
  @show H0/ny
  #@show Ld*k/2/π
  #@show cosh.(k*H0*0.5)./cosh.(k*H0)
  println()


  # Mesh
  function f_y(y, r, n, H0; dbgmsg = false)
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
  map(x) = VectorValue( x[1], f_y(x[2], mesh_ry, ny, H0) )
  model = CartesianDiscreteModel(domain,partition,map=map)


  # Labelling
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω,"surface",[3,4,6])   # assign the label "surface" to the entity 3,4 and 6 (top corners and top side)
  add_tag_from_tags!(labels_Ω,"bottom",[1,2,5])    # assign the label "bottom" to the entity 1,2 and 5 (bottom corners and bottom side)
  add_tag_from_tags!(labels_Ω,"inlet",[7])         # assign the label "inlet" to the entity 7 (left side)
  add_tag_from_tags!(labels_Ω,"outlet",[8])        # assign the label "outlet" to the entity 8 (right side)
  add_tag_from_tags!(labels_Ω, "water", [9])       # assign the label "water" to the entity 9 (interior)


  # Triangulations
  Ω = Interior(model) #same as Triangulation()
  Γ = Boundary(model,tags="surface") #same as BoundaryTriangulation()
  Γin = Boundary(model,tags="inlet")
  Γot = Boundary(model,tags="outlet")


  # Auxiliar functions
  function is_mem(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xm₀ <= x[1] <= xm₁ ) * ( x[2] ≈ 0.0)
  end


  # Masking and Beam Triangulation
  xΓ = get_cell_coordinates(Γ)
  Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
  Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
  Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask))
  Γη = Triangulation(Γ, findall(Γm_to_Γ_mask))
  Γκ = Triangulation(Γ, findall(!,Γm_to_Γ_mask))


  # Construct the tag for membrane boundary
  Λmb = Boundary(Γm)
  xΛmb = get_cell_coordinates(Λmb)
  xΛmb_n1 = findall(model.grid_topology.vertex_coordinates .== xΛmb[1])
  xΛmb_n2 = findall(model.grid_topology.vertex_coordinates .== xΛmb[2])
  new_entity = num_entities(labels_Ω) + 1
  labels_Ω.d_to_dface_to_entity[1][xΛmb_n1[1]] = new_entity
  labels_Ω.d_to_dface_to_entity[1][xΛmb_n2[1]] = new_entity
  add_tag!(labels_Ω, "mem_bnd", [new_entity])

  
  if vtk_output == true
    writevtk(model, fileName*"_model")
    writevtk(Ω,fileName*"_O")
    writevtk(Γ,fileName*"_G")
    writevtk(Γm,fileName*"_Gm")  
    writevtk(Γfs,fileName*"_Gfs")
    writevtk(Λmb,fileName*"_Lmb")  
  end


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓm = Measure(Γm,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  dΛmb = Measure(Λmb,degree)


  # Normals
  @show nΛmb = get_normal_vector(Λmb)


  # Dirichlet Fnc
  gη(x) = ComplexF64(0.0)

  # FE spaces
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64},
    dirichlet_tags=["mem_bnd"]) #diri
#   V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
#     vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  U_Γη = TrialFESpace(V_Γη, gη) #diri
#   U_Γη = TrialFESpace(V_Γη)

    
  # Weak form: Constant matrices
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  
  m11(η,v) = ∫( mᵨ*v*η )dΓm
  
  k11(η,v) = 
    ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓm +  
    ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri    
  
  k11Dry(η,v) = 
    ∫( Tᵨ*∇(v)⋅∇(η) )dΓm +
    ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri


  c12(ϕ,v) = ∫( v*ϕ )dΓm

  c21(η,w) = ∫( w*η )dΓm  

  k33(κ,u) = ∫( βh*u*g*κ )dΓfs  

  l1(v) = ∫( 0*v )dΓm
  l2(w) = ∫( 0*w )dΩ
  l3(u) = ∫( 0*u )dΓfs 
  l4(ξ) = ∫( 0*(ξ⋅î1) )dΩ 
  println("[MSG] Done Weak form")

  # Global matrices: constant matrices
  M11 = get_matrix(AffineFEOperator( m11, l1, U_Γη, V_Γη ))
  K11 = get_matrix(AffineFEOperator( k11, l1, U_Γη, V_Γη ))
  K11Dry = get_matrix(AffineFEOperator( k11Dry, l1, U_Γη, V_Γη ))
  C12 = get_matrix(AffineFEOperator( c12, l1, U_Ω, V_Γη ))

  C21 = get_matrix(AffineFEOperator( c21, l2, U_Γη, V_Ω ))

  K33 = get_matrix(AffineFEOperator( k33, l3, U_Γκ, V_Γκ ))


  println("[MSG] Done Global matrices")
  println(K11 == transpose(K11))

  #xp = range(xm₀, xm₁, size(V,2)+2)  

  ## Dry Natural Frequencies
  # --------------------Start--------------------
  λDry = LinearAlgebra.eigvals(M11\Matrix(K11Dry))
  VDry = LinearAlgebra.eigvecs(M11\Matrix(K11Dry))
  ωn_dry = sqrt.(λDry[1:nωₙ])
  da_V_dry = [ VDry[:,i] for i in [1:nωₙ;] ]
  meffDry = 
    diag( transpose(VDry[:,1:nωₙ]) * M11 * VDry[:,1:nωₙ] )
  
  dfDry = DataFrame(
    ωn = ωn_dry,
    V = da_V_dry,
    meff = meffDry
  )
  # ---------------------End---------------------


  ## Function to run the frequency analysis 
  #  for a given membrane mass and tension factor
  # --------------------Start--------------------
  function run_freq(ω, MassType = :Complex)

    ωreal = real(ω)
    αh = -im*ω/g*(1-βh)/βh #Disabled if βh=1.0
    k = dispersionRelAng(H0, ωreal; msg=false)
    @show ω, k, αh
  
    # Weak form: ω dependent
    k22(ϕ,w) = ∫( ∇(w)⋅∇(ϕ) )dΩ + ∫(-im*ω*ϕ * βh*αh*w)dΓfs +
        ∫( -w * im * k * ϕ )dΓin + ∫( -w * im * k * ϕ )dΓot
  
    c23(κ,w) = ∫( im*ω*w*κ )dΓfs + ∫(αh*βh*w*g*κ)dΓfs
  
    c32(ϕ,u) = ∫( -im*ω*u*ϕ*βh )dΓfs 
    
    # Global matrices: ω dependent
    K22 = get_matrix(AffineFEOperator( k22, l2, U_Ω, V_Ω ))
    C23 = get_matrix(AffineFEOperator( c23, l2, U_Γκ, V_Ω ))
    C32 = get_matrix(AffineFEOperator( c32, l3, U_Ω, V_Γκ ))
  
    # Solution
    tick()
    Mϕ = K22 - ( C23 * (Matrix(K33) \ C32) )
    Mhat = C12 * (Mϕ \ C21)
    Mtot = M11 + Mhat
    tock()    

    
    # # Eigen values memb only
    # λ = LinearAlgebra.eigvals(Mtot\Matrix(K11 + DTot))
    # V = LinearAlgebra.eigvecs(Mtot\Matrix(K11 + DTot))      
    # # @show real.(λ[1:nωₙ])
    # # ωₙ = sqrt.(real.(λ))

    # Eigen values System Complex
    if( MassType == :Complex )          
      println("[MSG] Using Complex Mass Matrix for Eigenvalue Problem")
      # Keeping the complex valued matrices
      Sol = Mtot \ K11 
    elseif( MassType == :Real )
      println("[MSG] Using Real Mass Matrix for Eigenvalue Problem")
      # Eigen values System Real only
      MtotReal = real.(Mtot)
      Sol = MtotReal \ K11
    end

    # Sol = MFull \ KFull
    λ = LinearAlgebra.eigvals(Sol)
    V = LinearAlgebra.eigvecs(Sol)
    
    # ω_all = sqrt.(λ)
    # @show λ_idx = sortperm(ω_all, by=real)
    # λ = λ[λ_idx]
    # V = V[:,λ_idx]

    meff = diag(transpose(V[:,1:nωₙ]) * Mtot * V[:,1:nωₙ])
    # meff = diag(transpose(V[:,1:nωₙ]) * MFullReal * V[:,1:nωₙ])

    # Wrong
    # Ur, S, Vr = svd(Mtot\Matrix(K11))    
    # λ = reverse(S)
    # V = reverse(Vr, dims=2)
    # rλ = real.(λ[1:nωₙ])
    # @show rλ    
    # return(rλ[1:nωₙ], V[:,1:nωₙ])

    @show sqrt.(λ[1:nωₙ])
    println()
    return(λ[1:nωₙ], V[:,1:nωₙ], meff)
  end
  # ---------------------End---------------------


  ## Wet Natural Frequencies Real
  # --------------------Start--------------------  
  da_ωₙ = zeros(Float64, nωₙ) 
  println("[MSG] Starting iterative solution for wet natural frequencies with Real Mass Matrix")
  println(ωn_guess)
  ωₙ = zeros(Float64, nωₙ) .+ real.(ωn_guess)
  da_V = []
  da_meff=[]
  da_iter = zeros(Int, nωₙ)

  
  for i in 1:nωₙ
    # global da_ωₙ, da_V  
    # global ωₙ, ω
    local V, lIter, meff
    lIter = 0    
    Δω = 1 
    ω = ωₙ[i]   
    ωₒ = ω 
    while ((Δω > 1e-3) && (lIter < maxIter))
      
      ωᵣ = αRelax * ω + (1 - αRelax) * ωₒ      
      # ωᵣ = real(ωᵣ)

      λ, V, meff = run_freq(ωᵣ, :Real)
      ωₒ = ω      
      ω = sqrt(λ[i])

      Δω = abs((ω - ωₒ)/ωₒ)
      lIter += 1
      # @show ωₙ
      @show i, ω, Δω, lIter
    end    

    da_ωₙ[i] = ω    
    da_iter[i] = lIter
    push!(da_V, V[:,i]) 
    push!(da_meff, meff[i])
  end

  dfWetReal = DataFrame(
    ωn = da_ωₙ,
    V = da_V,
    meff = da_meff,
    iter = da_iter
  )
  # ---------------------End---------------------



  ## Wet Natural Frequencies Complex
  # --------------------Start--------------------  
  da_ωₙ = zeros(ComplexF64, nωₙ)
  println("\n[MSG] Starting iterative solution for wet natural frequencies Complex")  
  ωₙ = dfWetReal[:,"ωn"] .+ 0.0*im
  da_V = []
  da_meff=[]
  da_iter = zeros(Int, nωₙ)

  
  for i in 1:nωₙ
    # global da_ωₙ, da_V  
    # global ωₙ, ω
    local V, lIter, meff
    lIter = 0    
    Δω = 1 
    ω = ωₙ[i]   
    ωₒ = ω 
    while ((Δω > 1e-3) && (lIter < maxIter))
      
      ωᵣ = αRelax * ω + (1 - αRelax) * ωₒ
      # ωᵣ = real(ωᵣ)

      λ, V, meff = run_freq(ωᵣ, :Complex)
      ωₒ = ω      
      ω = sqrt(λ[i])

      Δω = abs((ω - ωₒ)/ωₒ)
      lIter += 1
      # @show ωₙ
      @show i, ω, Δω, lIter
    end    

    da_ωₙ[i] = ω    
    da_iter[i] = lIter
    push!(da_V, V[:,i]) 
    push!(da_meff, meff[i])
  end

  dfWetComplex = DataFrame(
    ωn = da_ωₙ,
    V = da_V,
    meff = da_meff,
    iter = da_iter
  )
  # ---------------------End---------------------


  @show dfDry  
  @show dfWetReal
  @show dfWetComplex

  data = Dict(
    "params" => params,
    "dfDry" => dfDry,
    "dfWetReal" => dfWetReal,
    "dfWetComplex" => dfWetComplex
  )

  save(fileName*"_modesdata.jld2", data)
end


"""
Memb_params

Parameters for the VIV.jl module.
"""
@with_kw struct MembLR_params

  resDir::String = "data/sims_202508/mem_modes_free/"
  fileName::String = "mem"

  order::Int = 2
  vtk_output::Bool = true  

  H0 = 10 #m #still-water depth
  Lm = 2*H0 #m

  # Domain 
  nx = 120
  ny = 10
  mesh_ry = 1.2 #Ratio for Geometric progression of eleSize
  LΩ = 6*H0 
  x₀ = 0.0
  xm₀ = x₀ + 2*H0
  xm₁ = xm₀ + Lm

  mfac = 0.9
  tfac = 0.1

  # Number of natural frequencies
  nωₙ = 6

  # Iterative solution for wet natural frequencies
  αRelax = 0.8
  maxIter = 20
  
  # initial guess for natural frequencies
  ωn_guess = zeros(ComplexF64, nωₙ) .+ 1.0
  
  
end




end