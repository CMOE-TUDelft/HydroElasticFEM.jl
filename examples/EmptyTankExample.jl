module EmptyTankExample

using Gridap
using Parameters
using Printf
using WaveSpec

using HydroElasticFEM: PKG_ROOT, map_vertical_GP_for_const_dep
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.WaveInput_FrequencyDomain as WI

"""
Tutorial: Empty tank in the frequency domain
============================================

This example mirrors the legacy implementation in `src/EmptyTank2D/empt_freq_fnc.jl`
and then rewrites the same problem using the newer `HydroElasticFEM` code structure.

The problem is a 2D empty tank with:
- a linear free surface,
- inlet forcing from a monochromatic Airy wave,
- radiation boundary conditions on inlet and outlet,
- damping-zone partitions on the top boundary.

The two functions below solve the same monochromatic test case:
- `run_plain_implementation()` reproduces the explicit Gridap formulation.
- `run_structured_implementation()` uses `Geometry`, `Physics`, and `Simulation`.
"""

@with_kw struct EmptyTankTutorialParams
  H0::Float64 = 10.0
  ω::Float64 = 2.0
  η0::Float64 = 0.25
  α::Float64 = 0.0
  βₕ::Float64 = 0.5
  order::Int = 1
  nx::Int = 240
  ny::Int = 12
  mesh_ry::Float64 = 1.08
  probe_x::Vector{Float64} = collect(range(-40.0, 80.0, length=7))
end

function build_single_frequency_state(; H::Real, T::Real, h::Real, θ::Real = 0.0)
  spec = WaveSpec.ContinuousSpectrums.RegularWave(H, T)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(θ)
  ω = [2π / T]
  k = [WaveSpec.AiryWaves.solve_wavenumber(ω[1], h)]
  θ_vec = [θ]
  WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, ω, k, θ_vec, h, 1)
end

function tank_parameters(; H0::Real, nx::Integer, ny::Integer)
  Ld = 10.0 * H0
  LΩ = 2.0 * Ld + 3.0 * 2.0 * H0
  x0 = -Ld
  xd_in = 0.0
  xd_out = x0 + LΩ - Ld
  return (; Ld, LΩ, x0, xd_in, xd_out, domain=(x0, x0 + LΩ, -H0, 0.0), partition=(nx, ny))
end

function gp_map(mesh_ry, ny, H0)
  x -> VectorValue(x[1], map_vertical_GP_for_const_dep(x[2], mesh_ry, ny, H0; dbgmsg=false))
end

function shifted_gp_map(x0, mesh_ry, ny, H0)
  x -> VectorValue(
    x0 + x[1],
    map_vertical_GP_for_const_dep(x[2] - H0, mesh_ry, ny, H0; dbgmsg=false),
  )
end

function incident_wave(; H0::Real, ω::Real, η0::Real, α::Real)
  wave = WI.AiryWaveXZ(H0, ω, η0, α)
  ηin(x) = WI.surface_elevation(wave, x)
  ϕin(x) = WI.velocity_potential(wave, x)
  ∇ϕin(x) = WI.potential_gradient(wave, x)
  return (; wave, ηin, ϕin, ∇ϕin)
end

function probe_points(xs)
  Point.(xs, 0.0)
end

"""
    run_plain_implementation(; kwargs...)

Direct Gridap implementation following the same structure as
`src/EmptyTank2D/empt_freq_fnc.jl`.

This is the "everything spelled out explicitly" version:
- manual geometry tagging,
- manual top-boundary masking for free surface and damping zones,
- manual FE spaces and weak form definition,
- manual operator assembly and solve.
"""
function run_plain_implementation(; kwargs...)
  p = EmptyTankTutorialParams(; kwargs...)
  tp = tank_parameters(; H0=p.H0, nx=p.nx, ny=p.ny)
  inc = incident_wave(; H0=p.H0, ω=p.ω, η0=p.η0, α=p.α)
  probes = probe_points(p.probe_x)

  model = CartesianDiscreteModel(tp.domain, tp.partition, map=gp_map(p.mesh_ry, p.ny, p.H0))

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "surface", [3, 4, 6])
  add_tag_from_tags!(labels, "bottom", [1, 2, 5])
  add_tag_from_tags!(labels, "inlet", [7])
  add_tag_from_tags!(labels, "outlet", [8])
  add_tag_from_tags!(labels, "water", [9])

  Ω = Interior(model)
  Γ = Boundary(model, tags="surface")
  Γin = Boundary(model, tags="inlet")
  Γout = Boundary(model, tags="outlet")

  function is_damping_left(xs)
    c = sum(xs) / length(xs)
    (tp.x0 <= c[1] <= tp.xd_in) && (c[2] ≈ 0.0)
  end

  function is_damping_right(xs)
    c = sum(xs) / length(xs)
    (tp.xd_out <= c[1]) && (c[2] ≈ 0.0)
  end

  xΓ = get_cell_coordinates(Γ)
  Γd1_mask = lazy_map(is_damping_left, xΓ)
  Γd2_mask = lazy_map(is_damping_right, xΓ)
  Γd1 = Triangulation(Γ, findall(Γd1_mask))
  Γd2 = Triangulation(Γ, findall(Γd2_mask))
  Γfs = Triangulation(Γ, findall(!, Γd1_mask .| Γd2_mask))
  Γκ = Γ

  degree = 2 * p.order
  dΩ = Measure(Ω, degree)
  dΓfs = Measure(Γfs, degree)
  dΓd1 = Measure(Γd1, degree)
  dΓd2 = Measure(Γd2, degree)
  dΓin = Measure(Γin, degree)
  dΓout = Measure(Γout, degree)
  nΓin = get_normal_vector(Γin)

  reffe = ReferenceFE(lagrangian, Float64, p.order)
  VΩ = TestFESpace(Ω, reffe; conformity=:H1, vector_type=Vector{ComplexF64})
  VΓκ = TestFESpace(Γκ, reffe; conformity=:H1, vector_type=Vector{ComplexF64})
  UΩ = TrialFESpace(VΩ)
  UΓκ = TrialFESpace(VΓκ)
  X = MultiFieldFESpace([UΩ, UΓκ])
  Y = MultiFieldFESpace([VΩ, VΓκ])

  αₕ = -im * p.ω / WaveSpec.PhysicalConstants.g * (1.0 - p.βₕ) / p.βₕ
  k = inc.wave.k

  a((ϕ, κ), (w, u)) =
    ∫(∇(w) ⋅ ∇(ϕ))dΩ +
    ∫(p.βₕ * (u + αₕ * w) * (WaveSpec.PhysicalConstants.g * κ - im * p.ω * ϕ) + im * p.ω * w * κ)dΓfs +
    ∫(p.βₕ * (u + αₕ * w) * (WaveSpec.PhysicalConstants.g * κ - im * p.ω * ϕ) + im * p.ω * w * κ)dΓd1 +
    ∫(p.βₕ * (u + αₕ * w) * (WaveSpec.PhysicalConstants.g * κ - im * p.ω * ϕ) + im * p.ω * w * κ)dΓd2 +
    ∫(-im * k * w * ϕ)dΓin +
    ∫(-im * k * w * ϕ)dΓout

  l((w, u)) =
    ∫(w * (inc.∇ϕin ⋅ nΓin))dΓin -
    ∫(im * k * w * inc.ϕin)dΓin

  op = AffineFEOperator(a, l, X, Y)
  ϕₕ, κₕ = solve(op)

  κ_from_ϕ = (im * p.ω / WaveSpec.PhysicalConstants.g * ϕₕ)(probes)
  κ_direct = κₕ(probes)

  return (
    model=model,
    fields=(ϕ=ϕₕ, κ=κₕ),
    probes=probes,
    probe_surface_from_phi=κ_from_ϕ,
    probe_surface_direct=κ_direct,
    wave=inc.wave,
  )
end

"""
    run_structured_implementation(; kwargs...)

HydroElasticFEM version of the same empty-tank problem.

Compared with `run_plain_implementation`, this version delegates:
- mesh partitioning to `Geometry`,
- field definitions to `Physics`,
- operator assembly and solve to `Simulation`.

The damping zones are represented geometrically in `TankDomain2D`, while the
inlet forcing and radiation boundaries are attached to `PotentialFlow`.
"""
function run_structured_implementation(; kwargs...)
  p = EmptyTankTutorialParams(; kwargs...)
  tp = tank_parameters(; H0=p.H0, nx=p.nx, ny=p.ny)
  inc = incident_wave(; H0=p.H0, ω=p.ω, η0=p.η0, α=p.α)
  probes = probe_points(p.probe_x)

  tank = G.TankDomain2D(
    L=tp.LΩ,
    H=p.H0,
    nx=p.nx,
    ny=p.ny,
    map=shifted_gp_map(tp.x0, p.mesh_ry, p.ny, p.H0),
    damping_zones=[
      G.DampingZone1D(L=tp.Ld, x₀=[tp.x0, 0.0], domain_symbol=:Γd_1),
      G.DampingZone1D(L=tp.Ld, x₀=[tp.xd_out, 0.0], domain_symbol=:Γd_2),
    ],
  )

  sea_state = build_single_frequency_state(H=2.0 * p.η0, T=2π / p.ω, h=p.H0)
  inlet_traction(x) = (inc.∇ϕin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc.wave.k * inc.ϕin(x)

  potential = P.PotentialFlow(
    ρw=1025.0,
    g=WaveSpec.PhysicalConstants.g,
    sea_state=sea_state,
    boundary_conditions=[
      P.RadiationBC(domain=:dΓin),
      P.RadiationBC(domain=:dΓout),
      P.PrescribedInletPotentialBC(domain=:dΓin, forcing=inlet_traction, quantity=:traction),
    ],
    fe=PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64}),
    space_domain_symbol=:Ω,
  )

  freesurface = P.FreeSurface(
    ρw=1025.0,
    g=WaveSpec.PhysicalConstants.g,
    βₕ=p.βₕ,
    fe=PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64}),
    space_domain_symbol=:Γκ,
  )

  config = S.FreqDomainConfig(ω=p.ω)
  problem = S.build_problem(tank, P.PhysicsParameters[potential, freesurface], config)
  result = S.simulate(problem)

  ϕₕ, κₕ = result.solution

  return (
    problem=problem,
    result=result,
    fields=(ϕ=ϕₕ, κ=κₕ),
    probes=probes,
    probe_surface=κₕ(probes),
    wave=inc.wave,
  )
end

"""
    run_tutorial(; kwargs...)

Run both the plain and structured implementations and print a compact comparison
at the surface probes.
"""
function run_tutorial(; kwargs...)
  println("\n=== Empty Tank Tutorial ===")
  plain = run_plain_implementation(; kwargs...)
  structured = run_structured_implementation(; kwargs...)

  println("\nProbe comparison at y = 0:")
  for (pt, κ_plain, κ_struct) in zip(plain.probes, plain.probe_surface_direct, structured.probe_surface)
    @printf("x = %8.3f | plain |κ| = %10.4e | structured |κ| = %10.4e\n",
      pt[1], abs(κ_plain), abs(κ_struct))
  end

  return (; plain, structured)
end

if abspath(PROGRAM_FILE) == @__FILE__
  run_tutorial()
end

end
