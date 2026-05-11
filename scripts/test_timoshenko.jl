using Test, Gridap
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Geometry as D
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.ParameterHandler as FES

println("=== Struct traits ===")
beam = P.TimoshenkoBeam(E=210e9, ν=0.3, h_beam=0.1, b_beam=0.05, ρ_s=7800.0)
@assert P.variable_symbol(beam)  == :w
@assert P.variable_symbols(beam) == (:w, :θ)
@assert P.has_damping_form(beam) == false
cfgs = P.field_fe_configs(beam)
@assert cfgs[1].order == 2
@assert cfgs[2].order == 1
println("  OK")

println("=== Static solve (h/L = 0.1, <1% vs analytic) ===")
E, ν, κ = 210e9, 0.3, 5/6
L, q  = 1.0, 1e4
h, b  = 0.1*L, 0.05
Iv    = b*h^3/12
G     = E/(2*(1+ν))
EI    = E*Iv
κGA   = κ*G*b*h
w_ex  = 5*q*L^4/(384*EI) + q*L^2/(8*κGA)

model  = CartesianDiscreteModel((0.0,L),(20,))
Ω      = Triangulation(model)
dΩ     = Measure(Ω, 8)
dom    = D.IntegrationDomains(dΓη = dΩ)

rw = ReferenceFE(lagrangian, Float64, 2)
Vw = TestFESpace(model, rw; conformity=:H1, dirichlet_tags="boundary",
                 vector_type=Vector{Float64})
Uw = TrialFESpace(Vw, 0.0)
rθ = ReferenceFE(lagrangian, Float64, 1)
Vθ = TestFESpace(model, rθ; conformity=:H1, vector_type=Vector{Float64})
Uθ = TrialFESpace(Vθ)
X  = MultiFieldFESpace([Uw, Uθ])
Y  = MultiFieldFESpace([Vw, Vθ])
fmap = Dict(:w=>1, :θ=>2)

beam2 = P.TimoshenkoBeam(
  E=E, ν=ν, h_beam=h, b_beam=b, ρ_s=1.0, ρ_w=1.0, g=0.0, κ=κ,
  tangent=VectorValue(1.0),
  fe_w=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}),
  fe_θ=FES.FESpaceConfig(order=1, vector_type=Vector{Float64}),
)

a((w,θ),(vw,vθ)) = P.stiffness(beam2, dom,
  FO.FieldMap((w,θ),   fmap),
  FO.FieldMap((vw,vθ), fmap))

fmap_rhs = Dict(:w=>1)
src(x) = q
l((vw,vθ)) = P.rhs(beam2, dom,
  FO.FieldMap((src,),   fmap_rhs),
  FO.FieldMap((vw,vθ),  fmap))

op = AffineFEOperator(a, l, X, Y)
uh = solve(LUSolver(), op)
w_h = uh[1](Point(L/2))

err = abs(w_h - w_ex)/abs(w_ex)
println("  w_exact = $w_ex, w_h = $w_h, rel_err = $(round(err*100, sigdigits=3)) %")
@assert err < 0.01 "Error $(err) exceeds 1%"
println("  OK")

println("=== Thin limit (h/L=0.01, <2% from EB) ===")
h2   = 0.01*L
Iv2  = b*h2^3/12
EI2  = E*Iv2
w_EB = 5*q*L^4/(384*EI2)

model2 = CartesianDiscreteModel((0.0,L),(30,))
Ω2     = Triangulation(model2)
dΩ2    = Measure(Ω2, 8)
dom2   = D.IntegrationDomains(dΓη = dΩ2)

Vw2 = TestFESpace(model2, ReferenceFE(lagrangian,Float64,2);
      conformity=:H1, dirichlet_tags="boundary", vector_type=Vector{Float64})
Uw2 = TrialFESpace(Vw2, 0.0)
Vθ2 = TestFESpace(model2, ReferenceFE(lagrangian,Float64,1);
      conformity=:H1, vector_type=Vector{Float64})
Uθ2 = TrialFESpace(Vθ2)
X2  = MultiFieldFESpace([Uw2, Uθ2])
Y2  = MultiFieldFESpace([Vw2, Vθ2])

beam3 = P.TimoshenkoBeam(
  E=E, ν=ν, h_beam=h2, b_beam=b, ρ_s=1.0, ρ_w=1.0, g=0.0, κ=κ,
  tangent=VectorValue(1.0),
  fe_w=FES.FESpaceConfig(order=2, vector_type=Vector{Float64}),
  fe_θ=FES.FESpaceConfig(order=1, vector_type=Vector{Float64}),
)
a2((w,θ),(vw,vθ)) = P.stiffness(beam3, dom2,
  FO.FieldMap((w,θ),   fmap), FO.FieldMap((vw,vθ), fmap))
l2((vw,vθ)) = P.rhs(beam3, dom2,
  FO.FieldMap((src,), fmap_rhs), FO.FieldMap((vw,vθ), fmap))
op2 = AffineFEOperator(a2, l2, X2, Y2)
uh2 = solve(LUSolver(), op2)
w_h2 = uh2[1](Point(L/2))
err2 = abs(w_h2 - w_EB)/abs(w_EB)
println("  w_EB = $w_EB, w_h = $w_h2, rel_err = $(round(err2*100, sigdigits=3)) %")
@assert err2 < 0.02 "Error $(err2) exceeds 2%"
println("  OK")

println("\n=== All checks passed ===")
