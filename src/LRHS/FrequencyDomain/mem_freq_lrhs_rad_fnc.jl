include("./config_parameters.jl")

module MembLRHS2D

using Gridap
using JLD2
using Plots
using WaveSpec
using WaveSpec.Constants
using WaveSpec.Jonswap
using DataFrames:DataFrame
using DataFrames:Matrix
using TickTock
using Parameters
using Printf
using HydroElasticFEM.Resonator
using HydroElasticFEM: PKG_ROOT
using HydroElasticFEM.MeshModifier: map_vertical_GP_for_const_dep
import HydroElasticFEM.WaveInput_FrequencyDomain as WI


abstract type MemBndType end
struct Free <: MemBndType end
struct Fixed <: MemBndType end


function powerDissipatedResonator(ω, resonator, q, η)
  return 0.5*resonator.C*ω*ω* (abs(q - η))^2
end


function main(params)

  ## Function to run each freq
  # ---------------------Start---------------------
  function run_freq(ω, η₀, α)

    tick()
    airyWave = WI.AiryWaveXZ(H0, ω, η₀, α)

    ηin(x) = WI.surface_elevation(airyWave, x)
    ϕin(x) = WI.velocity_potential(airyWave, x)
    ∇ϕin(x) = WI.potential_gradient(airyWave, x)

    @show WI.wave_properties(airyWave)
    k = airyWave.k

    # Numeric constants
    αₕ = -im*ω/g * (1-βₕ)/βₕ
    @show αₕ
    println()  

    # Damping
    # Ldw = min( 5.0*λ, Ld )
    μ₀ = 2.5#maximum([2.5, 5.24/(ω^0.922)])#2.5
    μ₁ᵢₙ(x) =  μ₀*(1.0 - sin(π/2 * (x[1]-x₀)/Ld ))
    # μ₁ᵢₙ(x) =  μ₀*(1.0 - sin(π/2 * min( (x[1]-x₀)/Ldw, 1.0 ) ))
    μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*k
    # μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*k
    ηd(x) = μ₂ᵢₙ(x)*ηᵢₙ(x)
    ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???

    # Weak form
    ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)

    res_membrane((ϕ,κ,η),(w,u,v)) = 
      res_membrane((ϕ,κ,η),(w,u,v), memBndType)

    res_membrane((ϕ,κ,η),(w,u,v), ::Free) = 
      ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η
        - mᵨ*v*ω^2*η + Tᵨ*(1-im*ω*τ)*∇(v)⋅∇(η) )dΓm  
      
    res_membrane((ϕ,κ,η),(w,u,v), ::Fixed) = 
      ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η
        - mᵨ*v*ω^2*η + Tᵨ*(1-im*ω*τ)*∇(v)⋅∇(η) )dΓm  + 
      ∫(- Tᵨ*(1-im*ω*τ)*v*∇(η)⋅nΛmb )dΛmb

    function a(trialVars, testVars)
      (ϕ,κ,η,q...) = trialVars
      (w,u,v,ξ...) = testVars

      # # Weak form: Damping zone formulation
      # val = ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
      #   ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
      #   ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
      #     - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1    +
      #   ∫( -w * im * k * ϕ )dΓot +
      #   # ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
      #   #   - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2    +
      #   ∫(  v*(g*η - im*ω*ϕ) +  im*ω*w*η
      #     - mᵨ*v*ω^2*η + Tᵨ*(1-im*ω*τ)*∇(v)⋅∇(η) )dΓm  #+ 
      #   #∫(- Tᵨ*(1-im*ω*τ)*v*∇(η)⋅nΛmb )dΛmb
      
      val = res_membrane((ϕ,κ,η),(w,u,v)) +
        ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
        ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
        ∫( -w * im * airyWave.k * ϕ )dΓot +
        ∫( -w * im * airyWave.k * ϕ )dΓin        

      for (qi, ξi, δi, iresn) in zip(q, ξ, δ_p_Arr, resn)
        val +=
          (+im * ω * iresn.C - iresn.K) / ρw * δi(v * ((qi ⋅ î1) - η)) +
          ∫((ξi ⋅ qi) * 0.0)dΩ +
          # ∫( -rM/cnstFEArea*ω^2*(q⋅ξ) + rK/cnstFEArea*(ξ⋅q) )dΩ +
          -iresn.M * ω^2 * δi(qi ⋅ ξi) +
          (-im * ω * iresn.C + iresn.K) * δi(qi ⋅ ξi - (ξi ⋅ î1) * η)
      end

      return val
    end

    l((w,u,v)) =  ∫( w*(∇ϕin⋅nΓin) )dΓin + ∫( -w * im * airyWave.k * ϕin )dΓin


    # Solution
    op = AffineFEOperator(a,l,X,Y)
    (ϕₕ, κₕ, ηₕ, qₕ...) = solve(op)

    # Function for inlet phase
    κin = interpolate_everywhere(ηin, 
      FESpace(Γκ, reffe, conformity=:H1, vector_type=Vector{ComplexF64}))
    κr = κₕ - κin

    # Energy flux (Power) calculation
    ηx = ∇(ηₕ)⋅VectorValue(1.0,0.0)
    Pd = sum(∫( abs(ηx)*abs(ηx) )dΓm)
    Pd = 0.5*Tᵨ*ρw*τ*ω*ω*Pd    

    Pd_r = [ 
      powerDissipatedResonator(ω, iresn, qi(iresn.XZ)⋅î1, ηₕ(iresn.XZ)) 
      for (iresn, qi) in zip(resn, qₕ)
    ]
    # Pd_r = [0.0, 0.0]
    

    # Wave energy flux
    ηrf = abs(κr(Point(prbPowx[1],0.0)))
    ηtr = abs(κₕ(Point(prbPowx[2],0.0)))
    kh = k*H0
    wave_n = 0.5*(1 + 2*kh/sinh(2*kh))
    Pin = (0.5*ρw*g*η₀*η₀)*(ω/k)*wave_n
    Prf = (0.5*ρw*g*ηrf*ηrf)*(ω/k)*wave_n
    Ptr = (0.5*ρw*g*ηtr*ηtr)*(ω/k)*wave_n
    PErr = Pin - Prf - Ptr - Pd - sum( Pd_r )
    println("Power In \t ",Pin,"  W/m")
    println("Power Ref \t ",Prf," W/m")
    println("Power Trans \t ",Ptr," W/m")
    println("Power Abs \t ",Pd," W/m")
    println("Power Abs Resonator \t ",Pd_r," W")
    println("Error \t ",PErr," W/m")    

    # Interpolation on prboes
    prb_κ = zeros(ComplexF64, 1, length(prbxy))
    prb_κ_x = zeros(ComplexF64, 1, length(prbxy))  
    

    prb_κ[prbfs] = κₕ(prbxy[prbfs])
    prb_κ[prbmem] = ηₕ(prbxy[prbmem])

    prb_κ_x[prbfs] = (∇(κₕ)⋅VectorValue(1.0,0.0))(prbxy[prbfs])
    prb_κ_x[prbmem] = (∇(ηₕ)⋅VectorValue(1.0,0.0))(prbxy[prbmem])
  
    push!(prbDa, prb_κ)  
    push!(prbDa_x, prb_κ_x)  

    push!(prbDaΓη, ηₕ(prxΓη))
    push!(prbDaΓκ, κₕ(prxΓκ))

    push!(prbPow, [ω, Pin, Prf, Ptr, Pd, PErr, 0.0, Pd_r...])

    # VTK Output
    # ---------------------Start---------------------
    if vtk_output == true
      
      freqName = filename*"_omg_" * @sprintf("%.3f", ω)

      if(isdir(freqName))
        rm(freqName; recursive=true, force=true)
      end
      mkpath(freqName)


      writevtk(Ω, freqName * "/mem_O_sol.vtu",
        cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ),
        "phi_abs" => abs(ϕₕ), "phi_ang" => angle∘(ϕₕ)])

      
      writevtk(Γη, freqName * "/mem_R_sol.vtu",
        cellfields = vcat(
          [ "q$i" => real(qhi⋅î1)*VectorValue(1.0,0.0) + imag(qhi⋅î1)*VectorValue(0.0,1.0)
            for (i, qhi) in enumerate(qₕ) ],
          [ "q$i"*"_XZ" => iresn.XZ 
            for (i, iresn) in enumerate(resn) ],
          [ "q$i"*"_MKC" => VectorValue(iresn.M,iresn.K,iresn.C)
            for (i, iresn) in enumerate(resn) ]
        ) )
      

      writevtk(Γκ, freqName * "/mem_Gk_sol.vtu",
        cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ),
        "eta_abs" => abs(κₕ), "eta_ang" => angle∘(κₕ),
        "etaR_re" => real(κr),"etaR_im" => imag(κr),
        "etaR_abs" => abs(κr), "etaR_ang" => angle∘(κr),
        "ηin_abs" => abs(κin), "ηin_ang" => angle∘(κin)])

      writevtk(Γη, freqName * "/mem_Ge_sol.vtu",
        cellfields = ["eta_re" => real(ηₕ),"eta_im" => imag(ηₕ),
        "eta_abs" => abs(ηₕ), "eta_ang" => angle∘(ηₕ)])
    end
    # ----------------------End----------------------
    
    tock()
    println()
    return 0
  end
  # ----------------------End----------------------  


  @unpack name, order, vtk_output = params
  @show name
  @show order
  @show vtk_output
  filename = name*"/mem"

  @unpack H0, ω, T, η₀, α = params 
  k = dispersionRelAng.(H0, ω; msg=false)
  
  ρw = 1025 #kg/m3 water
  @show H0  #m #still-water depth
  @show ω

  # Peak Wave
  ωₚ, indp = findmax(η₀)
  @show ωₚ = ω[indp]
  kₚ = dispersionRelAng(H0, ωₚ; msg=false)
  println("Peak Wave T, L ", 2*pi/ωₚ, " ", 2*pi/kₚ)


  # Membrane parameters
  @unpack Lm, mᵨ, Tᵨ, τ = params
  @unpack memBndType = params
  @show Lm  #m
  @show g #defined in .Constants
  @show mᵨ #mass per unit area of membrane / ρw
  @show Tᵨ #T/ρw
  @show τ #damping coeff

  # Validate and convert memBndType to symbol
  memBndType = if memBndType == "free"
    Free()
  elseif memBndType == "fixed"
    Fixed()
  else
    error("memBndType should be either 'free' or 'fixed', got: ", memBndType)
  end
  @show memBndType  


  # Domain 
  @unpack nx, ny, mesh_ry, Ld, LΩ, x₀ = params
  @unpack domain, partition, xdᵢₙ, xm₀, xm₁ = params
  @show Lm
  @show LΩ, Ld
  @show domain
  @show partition
  @show mesh_ry
  @show (xm₀, xm₁)
  @show isinteger(Lm/LΩ*nx)
  @show LΩ/nx
  @show H0/ny
  println()


  # Numeric constants
  h = LΩ / nx
  γ = 1.0*order*(order-1)/h
  βₕ = 0.5
  # αₕ = -im*ω/g * (1-βₕ)/βₕ
  @show h
  @show βₕ
  # @show αₕ
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
  function is_mem(xs) # Check if an element is inside the beam1
    n = length(xs)
    x = (1/n)*sum(xs)
    (xm₀ <= x[1] <= xm₁ ) * ( x[2] ≈ 0.0)
  end
  # function is_damping1(xs) # Check if an element is inside the damping zone 1
  #   n = length(xs)
  #   x = (1/n)*sum(xs)
  #   (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ 0.0)
  # end
  # function is_damping2(xs) # Check if an element is inside the damping zone 2
  #   n = length(xs)
  #   x = (1/n)*sum(xs)
  #   (xdₒₜ <= x[1] ) * ( x[2] ≈ 0.0)
  # end

  # Masking and Beam Triangulation
  xΓ = get_cell_coordinates(Γ)
  Γm_to_Γ_mask = lazy_map(is_mem, xΓ)
  # Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
  # Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
  Γm = Triangulation(Γ, findall(Γm_to_Γ_mask))
  # Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
  # Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
  # Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask .| 
  #   Γd1_to_Γ_mask ))# .| Γd2_to_Γ_mask))
  Γfs = Triangulation(Γ, findall(!, Γm_to_Γ_mask ))    
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


  writevtk(model, filename*"_model")
  if vtk_output == true
    writevtk(Ω,filename*"_O")
    writevtk(Γ,filename*"_G")
    writevtk(Γm,filename*"_Gm")  
    # writevtk(Γd1,filename*"_Gd1")
    # writevtk(Γd2,filename*"_Gd2")
    writevtk(Γfs,filename*"_Gfs")
    writevtk(Λmb,filename*"_Lmb")  
  end


  # Measures
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓm = Measure(Γm,degree)
  # dΓd1 = Measure(Γd1,degree)
  # dΓd2 = Measure(Γd2,degree)
  dΓfs = Measure(Γfs,degree)
  dΓin = Measure(Γin,degree)
  dΓot = Measure(Γot,degree)
  dΛmb = Measure(Λmb,degree)


  # Normals
  @show nΛmb = get_normal_vector(Λmb)
  nΓin = get_normal_vector(Γin)


  # Dirichlet Fnc
  gη(x) = ComplexF64(0.0)

  # FE spaces
  # ---------------------Start---------------------
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
    vector_type=Vector{ComplexF64})
  U_Ω = TrialFESpace(V_Ω)
  U_Γκ = TrialFESpace(V_Γκ)

  if(memBndType == Fixed())
    V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
      vector_type=Vector{ComplexF64},
      dirichlet_tags=["mem_bnd"]) #diri
    U_Γη = TrialFESpace(V_Γη, gη)
  elseif(memBndType == Free())
    V_Γη = TestFESpace(Γη, reffe, conformity=:H1, 
      vector_type=Vector{ComplexF64})
    U_Γη = TrialFESpace(V_Γη)
  else
    error("memBndType should be either 'free' or 'fixed', got: ", memBndType)
  end  
  # ----------------------End----------------------


  # Resonator FE Spaces
  # ---------------------Start---------------------
  @unpack resn = params

  [@show iresn, iresn.XZ for iresn in resn]

  V_Γq_Arr = [ ConstantFESpace( Ω, 
    vector_type=Vector{ComplexF64}, 
    field_type=VectorValue{1,ComplexF64} ) 
    for iresn in resn ]
  U_Γq_Arr = [ TrialFESpace(iV_Γq) for iV_Γq in V_Γq_Arr ]
  î1 = VectorValue(1.0)

  δ_p_Arr = [ DiracDelta(Γ, iresn.XZ) for iresn in resn ]
  # ----------------------End----------------------

  X = MultiFieldFESpace([U_Ω, U_Γκ, U_Γη, U_Γq_Arr...])
  Y = MultiFieldFESpace([V_Ω, V_Γκ, V_Γη, V_Γq_Arr...])

  # Probes
  @unpack prbx, prbPowx = params
  # prbx=[  -20.0, 0.0, 20.0, 40.0, 50.0, 
  #         52.7, 53.7, 55, 60.0, 80.0, 
  #         85.0, 90.0, 95.0, 100.0, 120.0, 
  #         125.0, 140.0, 160.0, 180.0 ]
  # prbPowx=[ 55.0, 125.0 ]
  @show prbxy = Point.(prbx, 0.0)
  prbmem = (xm₀ .<= prbx .<= xm₁ )
  @show prbfs = findall(!,prbmem)
  @show prbmem = findall(prbmem)

  lDa = zeros(ComplexF64, 1, length(prbxy))
  prbDa = DataFrame(lDa, :auto) # for η
  prbDa_x = DataFrame(lDa, :auto) #for ηₓ


  # Storing soln at Γ
  # Difficulties in doing it purely using get_cell_dof_valules()
  # instead going to do it using evaluate
  # maybe a bit slow, but shouldnt matter too much i guess.
  xΓη = get_cell_coordinates(Γη)
  xΓκ = get_cell_coordinates(Γκ)
  prxΓη = [val[1] for val in xΓη]
  tmp = [val[2] for val in xΓη]
  push!(prxΓη,tmp[end])
  prxΓκ = [val[1] for val in xΓκ]
  push!(prxΓκ,prxΓη[1])
  sort!(prxΓκ)
  lDa = zeros(ComplexF64, 1, length(prxΓη))
  prbDaΓη = DataFrame(lDa, :auto)
  lDa = zeros(ComplexF64, 1, length(prxΓκ))
  prbDaΓκ = DataFrame(lDa, :auto)

  prbPow = DataFrame(zeros(Float64, 1, 7+length(resn)), :auto)


  # # Remove old vtk files  
  # for entry in readdir(name)
  #   if startswith(entry, "mem_omg")
  #     rm(joinpath(name, entry); force=true, recursive=true)
  #   end
  # end    

  # Run weak-form for each freq
  # ---------------------Start---------------------
  run_freq.(ω, η₀, α)
  # ----------------------End----------------------

  prbDa = prbDa[2:end, :]
  prbDa_x = prbDa_x[2:end, :]
  prbDaΓη = prbDaΓη[2:end,:]
  prbDaΓκ = prbDaΓκ[2:end,:]
  println("\nω, Pin, Prf, Ptr, Pd, PErr, 0.0, Pd_r...")
  @show prbPow = prbPow[2:end,:]

  k = dispersionRelAng.(H0, ω; msg=false)


  # Plotting 
  # ---------------------Start---------------------
  if isdir(filename*"_plots")
    rm(filename*"_plots"; force=true, recursive=true)
  end
  mkpath(filename*"_plots")  

  plt1 = plot(ω, prbPow[:,3]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "K_R",
    title = "Reflection coefficient",
    ylims = (0,1.0))

  plt2 = plot(ω, prbPow[:,4]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "K_T",
    title = "Transmission coefficient",
    ylims = (0,1.0))

  plt3 = plot(ω, prbPow[:,5]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "K_A",
    title = "Absorption coefficient",
    ylims = (0,1.0))

  plt4 = plot(ω, 100 * prbPow[:,6]./prbPow[:,2], linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "Error %",
    title = "Power Relative Error")

  pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
    plot_title = "Power Balance")
  savefig(pltAll,filename*"_plots/mem_powerBalance"*".png")

  for lprb in 1:length(prbxy)    

    plt1 = plot(ω, abs.(prbDa[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "A (m)",
      title = "Amplitude")  

    plt2 = plot(ω, abs.(prbDa_x[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "dA/dx",
      title = "Slope Magnitude")
    
    plt3 = plot(ω, angle.(prbDa[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "α (rad)",
      title = "Phase")  

    plt4 = plot(ω, angle.(prbDa_x[:,lprb]), linewidth=3, 
      xlabel = "ω (rad/s)",
      ylabel = "α (rad)",
      title = "Slope Phase")
    
    xloc = prbx[lprb]
    pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
      plot_title = "x = $xloc")

    savefig(pltAll,filename*"_plots/mem_dxPrb_$lprb"*".png")
  end  
  closeall() #close plots
  # ----------------------End----------------------

  data = Dict("ω" => ω,
              "η₀" => η₀,
              "k" => k,
              "prbxy" => prbxy,
              "prbDa" => prbDa,            
              "prbDa_x" => prbDa_x,
              "prxΓκ" => prxΓκ,
              "prxΓη" => prxΓη,
              "prbDaΓκ" => prbDaΓκ,
              "prbDaΓη" => prbDaΓη,
              "prbPow" => prbPow,
              "resn" => resn,
              "params" => params )

  save(filename*"_data.jld2", data)

end

end