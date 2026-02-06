module EBBeamExample
using HydroElasticFEM
using Gridap

model = CartesianDiscreteModel((0,1),(10,))
Ω = Triangulation(model)
Λ = Skeleton(Ω)
dΩ = Measure(Ω,4)
dΛ = Measure(Λ,4)
nΛ = get_normal_vector(Λ)
reffe = ReferenceFE(lagrangian,Float64,2)
V = TestFESpace(model,reffe,dirichlet_tags="boundary")
U = TrialFESpace(V)
f(x) = 1.0
EI = 100.0
γ_m = 10.0
h = 0.1
a(u,v) = ∫(EI * Δ(v)*Δ(u))dΩ +
         ∫(  EI * ( - jump(∇(v)⋅nΛ) * mean(Δ(u)) +
          -mean(Δ(v)) * jump(∇(u)⋅nΛ) + 
          γ_m/h*( jump(∇(v)⋅nΛ) * jump(∇(u)⋅nΛ) ) ) )dΛ 
l(v) = ∫(f * v)dΩ
op = AffineFEOperator(a,l,U,V)
solver = LUSolver()
uh = solve(solver,op)
writevtk(Ω,"scripts/EBBeam_example",nsubcells=10,cellfields=["uh"=>uh])


end