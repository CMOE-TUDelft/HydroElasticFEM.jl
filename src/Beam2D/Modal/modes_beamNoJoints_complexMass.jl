module BeamNoJointsModes

using Parameters
using JLD2
using Gridap
using Plots
using WaveSpec.Constants
using LinearAlgebra
using TickTock
using DataFrames
using Printf
using HydroElasticFEM: print_properties, Resonator, BeamNoJoints
using HydroElasticFEM: map_vertical_GP_for_const_dep


function run_case( params )  
  
  @printf("\n[MSG] Method 1: Complex Mass\n\n")

  # Constants
  @unpack ρw = params #kg/m3 water    
  @show g #defined in .Constants
    
  @unpack order, vtk_output = params
  @unpack H0 = params
  @unpack nωₙ = params

  @unpack αRelax, maxIter, ωn_guess = params

  @unpack resDir, fileName = params
  
  fileName = resDir*"/"*fileName

  # Beam parameters
  @unpack beam2D = params
  Lb = beam2D.L
  bndType = beam2D.bndType
  m_ρ = beam2D.m / ρw
  EI_ρ = beam2D.EI / ρw

  print_properties(beam2D)


  # Domain 
  @unpack nx, ny, mesh_ry, LΩ = params
  @unpack x0, xb0, xb1 = params
  domain =  (x0, x0+LΩ, -H0, 0.0)
  partition = (nx, ny)  
  @show Lb
  @show LΩ
  @show domain
  @show partition
  @show mesh_ry
  @show (xb0, xb1)
  @show isinteger(Lb/LΩ*nx)
  @show LΩ/nx
  @show H0/ny
  #@show Ld*k/2/π
  #@show cosh.(k*H0*0.5)./cosh.(k*H0)
  println()


  # Mesh
  map(x) = VectorValue(
    x[1],
    map_vertical_GP_for_const_dep(x[2], mesh_ry, ny, H0; dbgmsg=false)
  )
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
  function is_beam(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xb0 <= x[1] <= xb1 ) * ( x[2] ≈ 0.0)
  end
  function is_beam_node(coords) 
    n = length(coords)
    #println(n)
    x = (1/n)*sum(coords)
    #println(x)
    
    (xb0 < x[1] < xb1 ) * ( x[2] ≈ 0.0 )
  end


  # Masking and Beam Triangulation
  xΓ = get_cell_coordinates(Γ)
  Γb_to_Γ_mask = lazy_map(is_beam, xΓ)
  Γb = Triangulation(Γ, findall(Γb_to_Γ_mask))
  Γfs = Triangulation(Γ, findall(!, Γb_to_Γ_mask))
  Γη = Triangulation(Γ, findall(Γb_to_Γ_mask))
  Γκ = Triangulation(Γ, findall(!,Γb_to_Γ_mask))


  # Identify internal beam nodes for jump terms
  grid_dim_0_Γ = Skeleton(Γ) #Edges not included in the Skeleton
  xΓ_dim_0 = get_cell_coordinates(grid_dim_0_Γ)
  Λb_to_Γb_mask = lazy_map(is_beam_node,xΓ_dim_0)
  Λb = Triangulation(grid_dim_0_Γ,Λb_to_Γb_mask)
  
  if vtk_output == true
    isdir(fileName*"_vtk") || mkpath(fileName*"_vtk")
    writevtk(model, fileName*"_vtk/beam_model")
    writevtk(Ω,fileName*"_vtk/beam_O")
    writevtk(Γ,fileName*"_vtk/beam_G")
    writevtk(Γb,fileName*"_vtk/beam_Gb")  
    writevtk(Γfs,fileName*"_vtk/beam_Gfs")
    writevtk(Λb,fileName*"_vtk/beam_Lb")  
  end


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓb = Measure(Γb,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  dΛb = Measure(Λb,degree)


  # Normals
  @show nΛb = get_normal_vector(Λb)


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
  
  if(bndType == BeamNoJoints.Free())
    V_Γη = TestFESpace(Γη, reffe, conformity=:H1,
      vector_type=Vector{ComplexF64})
    U_Γη = TrialFESpace(V_Γη)      
  elseif(bndType == BeamNoJoints.Fixed())
    V_Γη = TestFESpace(Γη, reffe, conformity=:H1,
      vector_type=Vector{ComplexF64},
      dirichlet_tags=["mem_bnd"]) #diri
    U_Γη = TrialFESpace(V_Γη, gη) #diri
  else
    error("bndType should be either 'free' or 'fixed', got: ", bndType)
  end


  # Stabalisation
  h = LΩ/nx
  γ_m = 1.0*order*(order+1)
  println("\n[MSG] Stabalisation parameter γ_m = ", γ_m)
  println("[MSG] Stabalisation Element size h = ", h, "\n")


  ## Weak form: Constant matrices
  # --------------------Start--------------------
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  
  m11(η,v) = ∫( m_ρ*v*η )dΓb    


  c12(ϕ,v) = ∫( v*ϕ )dΓb

  c21(η,w) = ∫( -w*η )dΓb  

  c22_tmp(ϕ,w) =  ∫( w * ϕ )dΓin + ∫( w * ϕ )dΓot

  c23(κ,w) = ∫( -w*κ )dΓfs 

  c32(ϕ,u) = ∫( u*ϕ )dΓfs 


  k11(η,v) = k11(η,v,bndType)
  k11Dry(η,v) = k11Dry(η,v,bndType)
    
  k11(η,v,::BeamNoJoints.Free) = 
    ∫(  v*g*η + EI_ρ*Δ(v)*Δ(η) )dΓb +
    ∫(  EI_ρ * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) +
        -mean(Δ(v)) * jump(∇(η)⋅nΛb) + 
        γ_m/h*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) 
    )dΛb             
  
  k11Dry(η,v,::BeamNoJoints.Free) = 
    ∫(  EI_ρ*Δ(v)*Δ(η) )dΓb +
    ∫(  EI_ρ * ( - jump(∇(v)⋅nΛb) * mean(Δ(η)) +
        -mean(Δ(v)) * jump(∇(η)⋅nΛb) + 
        γ_m/h*( jump(∇(v)⋅nΛb) * jump(∇(η)⋅nΛb) ) ) 
    )dΛb             

  # k11(η,v,::BeamNoJoints.Fixed) = 
  #   ∫( v*g*η + Tᵨ*∇(v)⋅∇(η) )dΓb +  
  #   ∫(- Tᵨ*v*∇(η)⋅nΛb )dΛb #diri
  
  # k11Dry(η,v,::BeamNoJoints.Fixed) = 
  #   ∫( Tᵨ*∇(v)⋅∇(η) )dΓb +
  #   ∫(- Tᵨ*v*∇(η)⋅nΛb )dΛb #diri  
  

  k22(ϕ,w) = ∫( ∇(w)⋅∇(ϕ) )dΩ

  k33(κ,u) = ∫( u*g*κ )dΓfs  
  
  l1(v) = ∫( 0*v )dΓb
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

  ## DRY NATURAL FREQUENCIES
  # --------------------Start--------------------
  function run_dry_analysis()

    λDry, VDry = LinearAlgebra.eigen(M11\Matrix(K11Dry))
    ωn_dry = sqrt.(λDry[1:nωₙ])
    da_V_dry = [ VDry[:,i] for i in [1:nωₙ;] ]

    idx = sortperm(real.(ωn_dry))
    ωn_dry = ωn_dry[idx]
    da_V_dry = da_V_dry[idx]

    meffDry = 
      diag( transpose(VDry[:,1:nωₙ]) * M11 * VDry[:,1:nωₙ] )

    cache = (ωn_dry, da_V_dry, meffDry)
  end

  cache = run_dry_analysis()
  
  dfDry = DataFrame(
    ωn = cache[1],
    V = cache[2],
    meff = cache[3]
  )
  @show dfDry
  # ---------------------End---------------------


  ## WET NATURAL FREQUENCIES  
  # ===================================================

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
    # Sol = C12 * (Mϕ \ C21)  # This is the slow step

    # Faster implementation : 0.5s per call
    α = ComplexF64(-ω^2)
    β = ComplexF64(-im*k)
    Mϕ = (α .* M22) .+ (β .* C22_tmp) .+ K22
    Sol = C12 * (Mϕ \ Matrix(C21))

    # Version 1: Works, Complex Valued
    MTot = M11 - Sol
    KTot = K11

    AFull = MTot \ KTot
    
    λ, V = LinearAlgebra.eigen(AFull)

    ω_tmp = real.(sqrt.(λ))
    λ_idx = sortperm(abs.(ω_tmp))
    λ = λ[λ_idx]
    V = V[:,λ_idx]
    @show ω_tmp[1:nωₙ]

    cache = (MTot = MTot, K11 = K11)

    return λ, V, cache    
    
  end
  # ---------------------End---------------------

  ## Iterative Solution for each mode
  # --------------------Start--------------------
  da_ωₙ = zeros(ComplexF64, nωₙ) 
  println("[MSG] Starting iterative solution for wet natural frequencies")    
  da_V = []
  da_meff=[]
  da_iter = zeros(Int, nωₙ)

  startIndex = 1

  # if(bndType == BeamNoJoints.Free())
  #   startIndex = 2  # Skip first mode for free membrane
  
  #   # Push zeros in first mode for free membrane
  #   push!(da_V, zeros(ComplexF64,size(M11,1)))
  #   push!(da_meff, 0.0*im)
  # end

  for i in startIndex:nωₙ 
    
    local lIter, λ, VMode, ωc, meff, runTime

    lIter = 0    
    Δω = 1 
    ω = real.(dfDry.ωn[i])
    ωₒ = ω 
    while ((Δω > 1e-5) && (lIter < maxIter))
      
      runTime = time_ns()
      ωᵣ = αRelax * ω + (1 - αRelax) * ωₒ # Here ωᵣ is Real already
      ωᵣ = real(ωᵣ)

      λ, V, cache = run_freq(ωᵣ)
      
      ωₒ = ω      

      # Version 1: Works, Complex Valued
      ωc = sqrt(λ[i])
      ω = real(ωc)      

      VMode = V[:,i]
      meff = transpose(VMode) * cache.MTot * VMode      

      Δω = abs((ω - ωₒ)/ωₒ)
      lIter += 1
      runTime = (time_ns() - runTime)/1e9
      # @show ωₙ
      @show i, lIter, ω, Δω, runTime
      println()
    end        

    # Store results
    da_ωₙ[i] = ωc
    da_iter[i] = lIter
    push!(da_V, VMode) 
    push!(da_meff, meff)

    # Plot the eigenvalues for this mode iteration
    scatter(real.(λ), imag.(λ), label = "All Eigenvalues",
      title="Eigenvalues Mode $(i)", xlabel="Real(λ)", ylabel="Imag(λ)")
    scatter!([real.(λ[i])], [imag.(λ[i])], label = "Mode $(i)", 
      markersize=5, color=:red)    
    # plot!(xlim = (-5,5))
    
    isdir(fileName*"_figs") || mkpath(fileName*"_figs")
    savefig(fileName*"_figs/eigenvalues_mode_$(lpad(i, 3, '0')).png")

  end
  closeall() #close plots

  dfWet = DataFrame(
    ωn = da_ωₙ,
    V = da_V,
    meff = da_meff,
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



end