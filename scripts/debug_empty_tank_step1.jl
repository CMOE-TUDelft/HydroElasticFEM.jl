using Printf
using Gridap
using WaveSpec

include("../examples/EmptyTankExample.jl")
using .EmptyTankExample
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.ParameterHandler as PH

p = EmptyTankExample.EmptyTankTutorialParams(nx=80, ny=8, order=2)
tp = EmptyTankExample.tank_parameters(H0=p.H0, nx=p.nx, ny=p.ny)

# Plain geometry
model = CartesianDiscreteModel(tp.domain, tp.partition, map=EmptyTankExample.gp_map(p.mesh_ry, p.ny, p.H0))
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "surface", [3, 4, 6])
add_tag_from_tags!(labels, "bottom", [1, 2, 5])
add_tag_from_tags!(labels, "inlet", [7])
add_tag_from_tags!(labels, "outlet", [8])
Γ = Boundary(model, tags="surface")
Γin = Boundary(model, tags="inlet")
Γout = Boundary(model, tags="outlet")
dΓ = Measure(Γ, 4)
dΓin = Measure(Γin, 4)
dΓout = Measure(Γout, 4)
nin_plain = get_normal_vector(Γin)
nout_plain = get_normal_vector(Γout)

len_fs_plain = sum(∫(1.0)dΓ)
len_in_plain = sum(∫(1.0)dΓin)
len_out_plain = sum(∫(1.0)dΓout)
np_in = nin_plain(Point(tp.x0, -p.H0 / 2))
np_out = nout_plain(Point(tp.x0 + tp.LΩ, -p.H0 / 2))

# Structured geometry via problem assembly
inc = EmptyTankExample.incident_wave(H0=p.H0, ω=p.ω, η0=p.η0, α=p.α)
sea_state = EmptyTankExample.build_regular_wave_state(H=2.0 * p.η0, T=2 * pi / p.ω, h=p.H0)
f_in(x) = (inc.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc.sea_state.k[1] * inc.ϕin(x)

tank = G.TankDomain(L=tp.LΩ, H=p.H0, nx=p.nx, ny=p.ny, map=EmptyTankExample.shifted_gp_map(tp.x0, p.mesh_ry, p.ny, p.H0))
pf = P.PotentialFlow(
  ρw=1025.0,
  g=WaveSpec.PhysicalConstants.g,
  sea_state=sea_state,
  boundary_conditions=[
    P.RadiationBC(domain=:dΓin),
    P.RadiationBC(domain=:dΓout),
    P.PrescribedInletPotentialBC(domain=:dΓin, forcing=f_in, quantity=:traction),
  ],
  fe=PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64}),
  space_domain_symbol=:Ω,
)
fs = P.FreeSurface(
  ρw=1025.0,
  g=WaveSpec.PhysicalConstants.g,
  βₕ=p.βₕ,
  fe=PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64}),
  space_domain_symbol=:Γκ,
)
problem = S.build_problem(tank, P.PhysicsParameters[pf, fs], S.FreqDomainConfig(ω=p.ω))
dom = S.get_integration_domains(problem)
trians = S.get_triangulations(problem)

len_fs_struct = sum(∫(1.0)dom[:dΓκ])
len_in_struct = sum(∫(1.0)dom[:dΓin])
len_out_struct = sum(∫(1.0)dom[:dΓout])
nin_struct = get_normal_vector(trians[:Γin])
nout_struct = get_normal_vector(trians[:Γout])
ns_in = nin_struct(Point(tp.x0, -p.H0 / 2))
ns_out = nout_struct(Point(tp.x0 + tp.LΩ, -p.H0 / 2))

@printf("plain     : |Γfs|=%.8f |Γin|=%.8f |Γout|=%.8f n_in=(%.4f,%.4f) n_out=(%.4f,%.4f)\n", len_fs_plain, len_in_plain, len_out_plain, np_in[1], np_in[2], np_out[1], np_out[2])
@printf("structured: |Γκ| =%.8f |Γin|=%.8f |Γout|=%.8f n_in=(%.4f,%.4f) n_out=(%.4f,%.4f)\n", len_fs_struct, len_in_struct, len_out_struct, ns_in[1], ns_in[2], ns_out[1], ns_out[2])
