module MembraneModes

using Parameters
using JLD2
using Gridap
using Plots
using WaveSpec.Constants
using LinearAlgebra
using TickTock
using DataFrames
using Printf

abstract type MemBndType end
struct Free <: MemBndType end
struct Fixed <: MemBndType end

function run_case( params )

  @printf("\n[MSG] Method 2: Damped System\n\n")
    
  @unpack mfac, tfac = params
  @unpack memBndType = params
  @unpack order, vtk_output = params
  @unpack H0, Lm = params
  @unpack nωₙ = params

  @unpack αRelax, maxIter, ωn_guess = params

  @unpack resDir, fileName = params
  
  fileName = resDir*"/"*fileName

  # Validate and convert memBndType to symbol
  memBndType = if memBndType == "free"
    Free()
  elseif memBndType == "fixed"
    Fixed()
  else
    error("memBndType should be either 'free' or 'fixed', got: ", memBndType)
  end
  @show memBndType  

  # Constants
  ρw = 1025 #kg/m3 water    
  @show g #defined in .Constants
  
  # Membrane parameters
  mᵨ = mfac  #mass per unit area of membrane / ρw
  Tᵨ = tfac *g*H0*H0 #T/ρw

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
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)
  
  if(memBndType == Free())
    V_Γη = TestFESpace(Γη, reffe, conformity=:H1,
      vector_type=Vector{ComplexF64})
    U_Γη = TrialFESpace(V_Γη)      
  elseif(memBndType == Fixed())
    V_Γη = TestFESpace(Γη, reffe, conformity=:H1,
      vector_type=Vector{ComplexF64},
      dirichlet_tags=["mem_bnd"]) #diri
    U_Γη = TrialFESpace(V_Γη, gη) #diri
  else
    error("memBndType should be either 'free' or 'fixed', got: ", memBndType)
  end


  ## Weak form: Constant matrices
  # --------------------Start--------------------
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  
  m11(η,v) = ∫( mᵨ*v*η )dΓm  


  c12(ϕ,v) = ∫( v*ϕ )dΓm

  c21(η,w) = ∫( -w*η )dΓm  

  c22_tmp(ϕ,w) =  ∫( w * ϕ )dΓin + ∫( w * ϕ )dΓot

  c23(κ,w) = ∫( -w*κ )dΓfs 

  c32(ϕ,u) = ∫( u*ϕ )dΓfs 


  k11(η,v) = k11(η,v,memBndType)
  k11Dry(η,v) = k11Dry(η,v,memBndType)
    
  k11(η,v,::Free) = 
    ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓm #+  
    # ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri    
  
  k11Dry(η,v,::Free) = 
    ∫( Tᵨ*∇(v)⋅∇(η) )dΓm #+
    # ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri

  k11(η,v,::Fixed) = 
    ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓm +  
    ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri    
  
  k11Dry(η,v,::Fixed) = 
    ∫( Tᵨ*∇(v)⋅∇(η) )dΓm +
    ∫(- Tᵨ*v*∇(η)⋅nΛmb )dΛmb #diri

  
  k22(ϕ,w) = ∫( ∇(w)⋅∇(ϕ) )dΩ

  k33(κ,u) = ∫( u*g*κ )dΓfs  

  l1(v) = ∫( 0*v )dΓm
  l2(w) = ∫( 0*w )dΩ
  l3(u) = ∫( 0*u )dΓfs 
  println("[MSG] Done Weak form")

  # Global matrices: constant matrices
  M11 = get_matrix(AffineFEOperator( m11, l1, U_Γη, V_Γη ))
  
  C12 = get_matrix(AffineFEOperator( c12, l1, U_Ω, V_Γη ))
  C21 = get_matrix(AffineFEOperator( c21, l2, U_Γη, V_Ω ))
  C22_tmp = get_matrix(AffineFEOperator( c22_tmp, l2, U_Ω, V_Ω ))
  C23 = get_matrix(AffineFEOperator( c23, l2, U_Γκ, V_Ω ))
  C32 = get_matrix(AffineFEOperator( c32, l3, U_Ω, V_Γκ ))

  K11 = get_matrix(AffineFEOperator( k11, l1, U_Γη, V_Γη ))
  K11Dry = get_matrix(AffineFEOperator( k11Dry, l1, U_Γη, V_Γη ))
  K22 = get_matrix(AffineFEOperator( k22, l2, U_Ω, V_Ω ))
  K33 = get_matrix(AffineFEOperator( k33, l3, U_Γκ, V_Γκ ))  

  # Eliminate η equation
  M22 = - C23 * (Matrix(K33) \ C32)

  println("[MSG] Done Global matrices")  
  # ----------------------End---------------------

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


  ## WET NATURAL FREQUENCIES  
  ## Nonlinear System
  # --------------------Start--------------------
  function run_freq(ω)
    """
    Run one iteration of the frequency calculation
    Here ω is Real
    """
    
    if(ω isa Complex)
      error("ω should be Real, but got Complex ω = ", ω)    
    end

    k = dispersionRelAng(H0, ω; msg=false)
    
    # # Easy implementation : Very slow 3s per call
    # Mϕ = -ω^2 * M22 - im*k*C22_tmp + K22
    # Sol = C12 * (Mϕ \ C21) # This is the slow step

    # Faster implementation : 0.5s per call
    α = ComplexF64(-ω^2)
    β = ComplexF64(-im*k)
    Mϕ = (α .* M22) .+ (β.*C22_tmp) .+ K22
    Sol = C12 * (Mϕ \ Matrix(C21))

    ## Added mass matrix
    A = - real( Sol )
    ## Radiation damping matrix
    C = - imag( Sol )*ω

    # # Version 1: Works, Complex Valued
    # MTot = M11 - Sol
    # KTot = K11
    # AFull = MTot \ KTot

    # λ = LinearAlgebra.eigvals(AFull)

    # return λ

    # Version 2: Damped system
    MTot = M11 + A
    CTot = C
    KTot = K11

    sz = size(MTot,1)
    AFull = [ zeros(sz, sz)           I(sz);
              -MTot\KTot        -MTot\CTot ]

    AFull = Matrix(AFull)    

    λ, V = LinearAlgebra.eigen(AFull)
    
    λ_idx = sortperm(abs.(imag.(λ)))    
    λ = λ[λ_idx]
    V = V[:, λ_idx]
    
    # @show λ
    return λ, V
    
  end
  # ---------------------End---------------------
  

  #  Iterative Solution for each mode
  # --------------------Start--------------------
  da_ωnneg = []
  da_ωnpos = []
  println("[MSG] Starting iterative solution for wet natural frequencies")    
  da_Vneg = []
  da_Vpos = []
  da_meff=[]
  da_iter = zeros(Int, nωₙ)

  startIndex = 1

  if(memBndType == Free())
    startIndex = 2 # First mode is zero frequency mode for free membrane
  
    # Push zeros in first mode for free membrane
    push!(da_ωnneg, 0.0im)
    push!(da_ωnpos, 0.0im)
    push!(da_Vneg, zeros(ComplexF64,2*size(M11,1)) )
    push!(da_Vpos, zeros(ComplexF64,2*size(M11,1)) )
  end

  for i in startIndex:nωₙ

    local V, lIter, meff, λ, cache, runTime
    local ω_save, VMode_neg, VMode_pos

    lIter = 0    
    Δω = 1 
    ω = dfDry.ωn[i]
    # If the dry freq is too close to zero, use the guess instead
    if( ω < 0.1 ) ω = ωn_guess[i] end 
    @printf("\nWet Modes %3d, Starting ω_guess =%04f \n",i, ω)
    ωₒ = ω     
    while ((Δω > 1e-5) && (lIter < maxIter))
      
      runTime = time_ns()
      ωᵣ = αRelax * ω + (1 - αRelax) * ωₒ      
      ωᵣ = real(ωᵣ)
      
      λ, V = run_freq(ωᵣ)
      ωₒ = ω      

      jneg = imag(λ[2*i]) < 0 ? 2i : 2i-1 # Mode with neg imag part
      jpos = jneg == 2i ? 2i-1 : 2i       # Mode with pos imag part

      # Version 2: Damped system
      ω = -imag(λ[jneg])

      ω_save = [λ[jneg], λ[jpos]]  # Complex conjugate pair      
      VMode_neg = V[:, jneg]
      VMode_pos = V[:, jpos]

      Δω = abs((ω - ωₒ)/ωₒ)
      lIter += 1
      runTime = (time_ns() - runTime)/1e9
      # @show ωₙ
      @show i, lIter, ω, Δω, runTime
    end        

    # Store results
    push!(da_ωnneg, ω_save[1])
    push!(da_ωnpos, ω_save[2])
    da_iter[i] = lIter
    push!(da_Vneg, VMode_neg) 
    push!(da_Vpos, VMode_pos)
    # push!(da_meff, meff)

    # Plot the eigenvalues for this mode iteration
    scatter(real.(λ), imag.(λ), label = "All Eigenvalues",
      title="Eigenvalues Mode $(i)", xlabel="Real(λ)", ylabel="Imag(λ)")
    scatter!(real.(ω_save), imag.(ω_save), label = "Mode $(i)", 
      markersize=5, color=:red)    
    # plot!(xlim = (-5,5))

    isdir(fileName*"_figs") || mkpath(fileName*"_figs")
    savefig(fileName*"_figs/eigenvalues_mode_$(lpad(i, 3, '0')).png")
  end
  closeall() #close plots

  dfWet = DataFrame(
    ωnneg = da_ωnneg,
    ωnpos = da_ωnpos,
    Vneg = da_Vneg,
    Vpos = da_Vpos,
    # meff = da_meff,
    iter = da_iter
  )
  # ---------------------End---------------------

  

  ## Print and Save Output
  # --------------------Start--------------------
  @show dfDry  
  @show dfWet

  data = Dict(
    "params" => params,
    "dfDry" => dfDry,
    "dfWet" => dfWet
  )

  save(fileName*"_modesdata.jld2", data)
  # ---------------------End---------------------
end


"""
Memb_params

Parameters for the VIV.jl module.
"""
@with_kw struct Memb_params

  resDir::String = "data/sims_202508/mem_modes_free/"
  fileName::String = "mem"

  order::Int = 2
  vtk_output::Bool = true  

  memBndType::String = "free"  #"fixed" or "free"
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