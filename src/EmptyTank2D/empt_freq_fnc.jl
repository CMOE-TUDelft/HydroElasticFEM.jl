module EmptTank2D_Freq

using Gridap
using Plots
using Printf
using JLD2
using WaveSpec
using .Constants
using .Jonswap
using DataFrames:DataFrame
using DataFrames:Matrix
import HydroElasticFEM.WaveInput_FrequencyDomain as WI

function run_freq(ω, η₀, α)
  
  airyWave = WI.AiryWaveXZ(H0, ω, η₀, α)

  ηin(x) = WI.surface_elevation(airyWave, x)
  ϕin(x) = WI.velocity_potential(airyWave, x)
  ∇ϕin(x) = WI.potential_gradient(airyWave, x)

  @show WI.wave_properties(airyWave)

  # Numeric constants
  αₕ = -im*ω/g * (1-βₕ)/βₕ
  @show αₕ
  println()  

  # Damping
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x) * airyWave.k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x) * airyWave.k
  ηd(x) = μ₂ᵢₙ(x)*ηin(x)
  ∇ₙϕd(x) = μ₁ᵢₙ(x)*vzfsᵢₙ(x) #???  

  # Weak form
  ∇ₙ(ϕ) = ∇(ϕ)⋅VectorValue(0.0,1.0)
  a((ϕ,κ),(w,u)) =      
    ∫(  ∇(w)⋅∇(ϕ) )dΩ   +
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓfs   +
    # ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
    #   - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd1  +
    # ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ 
    #   - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w) )dΓd2 
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓd1  +  
    ∫(  βₕ*(u + αₕ*w)*(g*κ - im*ω*ϕ) + im*ω*w*κ )dΓd2 +
    ∫( -w * im * airyWave.k * ϕ )dΓot +
    ∫( -w * im * airyWave.k * ϕ )dΓin 

  # l((w,u)) =  ∫( w*vxᵢₙ )dΓin - ∫( ηd*w - ∇ₙϕd*(u + αₕ*w) )dΓd1
  # l((w,u)) =  ∫( w*vxᵢₙ )dΓin - ∫( w * im * k * ϕin )dΓin
  l((w,u)) =  ∫( w*(∇ϕin⋅nΓin) )dΓin - ∫( w * im * airyWave.k * ϕin )dΓin
  
  # Solution
  op = AffineFEOperator(a,l,X,Y)
  (ϕₕ,κₕ) = solve(op)

  paraFolder = name*"/mem_omg_"*@sprintf("%.2f",ω)
  paraFile = paraFolder*"/mem"

  # VTK output for each frequency
  if vtk_output == true
    if( isdir(paraFolder) )
      rm(paraFolder; recursive=true)
    end
    mkdir(paraFolder)

    writevtk(Ω,paraFile * "_O_sol.vtu",
      cellfields = ["phi_re" => real(ϕₕ),"phi_im" => imag(ϕₕ),
      "phi_abs" => abs(ϕₕ), "phi_ang" => angle∘(ϕₕ)])
    writevtk(Γκ,paraFile * "_Gk_sol.vtu",
      cellfields = ["eta_re" => real(κₕ),"eta_im" => imag(κₕ),
      "eta_abs" => abs(κₕ), "eta_ang" => angle∘(κₕ)])
  end

  # Because unable to do interpolation on κₕ
  prb_κ = im*ω/g*ϕₕ(prbxy)     

  prb_κ_x = (im*ω/g * (∇(ϕₕ)⋅VectorValue(1.0,0.0)) )(prbxy)
 
  push!(prbDa, prb_κ)  
  push!(prbDa_x, prb_κ_x)  

  return 0
end


name::String = "data/sims_202601/empt_freq_spec"
order::Int = 2
vtk_output::Bool = true
filename = name*"/mem"

H0 = 10 #m #still-water depth

# Wave parameters
# ω, S, η₀ = jonswap(0.4, 2.5; 
#     plotflag=true, plotloc=filename, nω=145)
# println(ω[1], "\t", ω[2], "\t", ω[end])
# ω = ω[2:end]
# S = S[2:end]
# η₀ = η₀[2:end]
# ω = [2*π/2.53079486745378, 2*π/2.0]
# η₀ = [0.25, 0.25]
ω = 1:0.2:3
T = 2*π./ω
η₀ = 0.25*ones(length(ω))
α = randomPhase(ω; seed=100)


# Peak Wave
ωₚ, indp = findmax(η₀)
@show ωₚ = ω[indp]
kₚ = dispersionRelAng(H0, ωₚ; msg=false)
println("Peak Wave T, L ", 2*pi/ωₚ, " ", 2*pi/kₚ)


# Domain 
nx = 2600
ny = 20
mesh_ry = 1.1 #Ratio for Geometric progression of eleSize
Ld = 10*H0 #damping zone length
LΩ = 2*Ld + 3*2*H0
x₀ = -Ld
domain =  (x₀, x₀+LΩ, -H0, 0.0)
partition = (nx, ny)
xdᵢₙ = 0.0
xdₒₜ = x₀ + LΩ - Ld
@show LΩ
@show domain
@show partition
@show LΩ/nx
@show H0/ny
println()


# Numeric constants
h = LΩ / nx
γ = 1.0*order*(order-1)/h
βₕ = 0.5
#αₕ = -im*ω/g * (1-βₕ)/βₕ
@show h
@show βₕ
#@show αₕ
println()


# Damping
μ₀ = 2.5
μ₁ᵢₙ(x) = μ₀*(1.0 - sin(π/2*(x[1]-x₀)/Ld))
μ₁ₒᵤₜ(x) = μ₀*(1.0 - cos(π/2*(x[1]-xdₒₜ)/Ld))
# μ₂ᵢₙ(x) = μ₁ᵢₙ(x)*kₚ
# μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x)*kₚ


# Mesh
function f_y(y, r, n, H0)
  # Mesh along depth as a GP
  # Depth is 0 to -H0
  if(r ≈ 1.0)
    return y  
  else
    a0 = H0 * (r-1) / (r^n - 1)    
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
function is_damping1(xs) # Check if an element is inside the damping zone 1
  n = length(xs)
  x = (1/n)*sum(xs)
  (x₀ <= x[1] <= xdᵢₙ ) * ( x[2] ≈ 0.0)
end
function is_damping2(xs) # Check if an element is inside the damping zone 2
  n = length(xs)
  x = (1/n)*sum(xs)
  (xdₒₜ <= x[1] ) * ( x[2] ≈ 0.0)
end

# Masking damping zones
xΓ = get_cell_coordinates(Γ)
Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
Γfs = Triangulation(Γ, findall(!, Γd1_to_Γ_mask .| Γd2_to_Γ_mask))
Γκ = Γ

writevtk(model, filename*"_model")
if vtk_output == true
  writevtk(Ω,filename*"_O")
  writevtk(Γ,filename*"_G")
  writevtk(Γd1,filename*"_Gd1")
  writevtk(Γd2,filename*"_Gd2")
  writevtk(Γfs,filename*"_Gfs")
end


# Measures
degree = 2*order
dΩ = Measure(Ω,degree)
dΓd1 = Measure(Γd1,degree)
dΓd2 = Measure(Γd2,degree)
dΓfs = Measure(Γfs,degree)
dΓin = Measure(Γin,degree)
dΓot = Measure(Γot,degree)

# Normals
nΓin = get_normal_vector(Γin)

# FE spaces
reffe = ReferenceFE(lagrangian,Float64,order)
V_Ω = TestFESpace(Ω, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1, 
  vector_type=Vector{ComplexF64})
U_Ω = TrialFESpace(V_Ω)
U_Γκ = TrialFESpace(V_Γκ)
X = MultiFieldFESpace([U_Ω,U_Γκ])
Y = MultiFieldFESpace([V_Ω,V_Γκ])


# Probes
prbx=range(-80, 140, 23)
@show prbxy = Point.(prbx, 0.0)

lDa = zeros(ComplexF64, 1, length(prbxy))
prbDa = DataFrame(lDa, :auto) # for η
prbDa_x = DataFrame(lDa, :auto) #for ηₓ

run_freq.(ω, η₀, α)

@show prbDa = prbDa[2:end, :]
prbDa_x = prbDa_x[2:end, :]

## Plotting and saving data
# ---------------------Start---------------------
if isdir(filename*"_plots")
  rm(filename*"_plots"; force=true, recursive=true)
end
mkpath(filename*"_plots")  

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

  # plt4 = plot(ω, angle.(prbDa_x[:,lprb]), linewidth=3, 
  #   xlabel = "ω (rad/s)",
  #   ylabel = "α (rad)",
  #   title = "Slope Phase")

  k = dispersionRelAng.(Ref(H0), ω; msg=false)
  dByL = H0 ./ (2π ./ k)
  plt4 = plot(ω, dByL, linewidth=3, 
    xlabel = "ω (rad/s)",
    ylabel = "",
    title = "d / L")
  
  xloc = prbx[lprb]
  pltAll = plot(plt1, plt2, plt3, plt4, layout=4, dpi=330,
    plot_title = "x = $xloc")

  savefig(pltAll,filename*"_plots/mem_dxPrb_$lprb"*".png")
end
# ----------------------End----------------------

data = Dict("ω" => ω,
            "η₀" => η₀,
            "prbxy" => prbxy,
            "prbDa" => prbDa,
            "prbDa_x" => prbDa_x)

save(filename*"_data.jld2", data)

end