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

"""
Tutorial: Empty tank in the frequency domain
============================================

This example mirrors the plain implementation using Gridap functions
and then rewrites the same problem using the `HydroElasticFEM` package.

The problem is a 2D empty tank with:
- a linear free surface,
- inlet forcing from a monochromatic Airy wave,
- radiation boundary conditions on inlet and outlet.

The two functions below solve the same monochromatic test case:
- `run_plain_implementation()` reproduces the explicit Gridap formulation.
- `run_structured_implementation()` uses `Geometry`, `Physics`, and `Simulation` modules from the `HydroElasticFEM` package.
"""

"""
    EmptyTankTutorialParams

Parameters for the empty tank tutorial.
- `H0` — still water depth
- `ω` — wave frequency
- `η0` — wave height (used to build the incident wave state)
- `α` — wave angle of incidence (0 = head-on)
- `βₕ` — free-surface damping parameter (0.5 = critical damping)
- `order` — FE order
- `nx`, `ny` — mesh resolution parameters
- `mesh_ry` — mesh grading parameter in the vertical direction (1.0 = uniform, >1 = finer near surface)
- `probe_x` — x-coordinates of surface probes for solution comparison between implementations
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
  probe_x::Vector{Float64} = collect(range(-10.0, 50.0, length=7))
end

"""
    build_regular_wave_state(; H, T, h, θ=0.0)

Builds a regular wave state using the `WaveSpec` package.
- `H` — wave height
- `T` — wave period
- `h` — water depth
- `θ` — wave angle of incidence (default = 0.0)
"""
function build_regular_wave_state(; H::Real, T::Real, h::Real, θ::Real = 0.0)
  spec = WaveSpec.ContinuousSpectrums.RegularWave(H, T)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(θ)
  ω = [2π / T]
  k = [WaveSpec.AiryWaves.solve_wavenumber(ω[1], h)]
  θ_vec = [θ]
  WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, ω, k, θ_vec, h, 1)
end

"""
    tank_parameters(; H0, nx, ny)

Helper function to compute tank geometry parameters from depth and mesh resolution.
- `H0` — still water depth
- `nx`, `ny` — mesh resolution parameters
Returns a named tuple with:
- `LΩ` — total tank length (set to 6 times the depth for this example)
- `x0` — x-coordinate of the left boundary
- `domain` — tuple defining the rectangular domain (x0, x0 + LΩ, -H0, 0)
- `partition` — tuple defining the mesh partition (nx, ny)
"""
function tank_parameters(; H0::Real, nx::Integer, ny::Integer)
  LΩ = 6.0 * H0
  x0 = -H0
  return (; LΩ, x0, domain=(x0, x0 + LΩ, -H0, 0.0), partition=(nx, ny))
end

"""
    gp_map(mesh_ry, ny, H0)

Mapping function for the mesh grading in the vertical direction.
- `mesh_ry` — grading parameter (1.0 = uniform, >1 = finer near surface)
- `ny` — number of vertical elements
- `H0` — still water depth
Returns a function that maps reference coordinates to physical coordinates with the specified grading.
"""
function gp_map(mesh_ry, ny, H0)
  x -> VectorValue(x[1], map_vertical_GP_for_const_dep(x[2], mesh_ry, ny, H0; dbgmsg=false))
end

"""
    shifted_gp_map(x0, mesh_ry, ny, H0)

Mapping function for the mesh grading in the vertical direction, with a horizontal shift.
- `x0` — horizontal shift (x-coordinate of the left boundary)
- `mesh_ry` — grading parameter (1.0 = uniform, >1 = finer near surface)
- `ny` — number of vertical elements
- `H0` — still water depth
Returns a function that maps reference coordinates to physical coordinates with the specified grading and horizontal shift.
"""
function shifted_gp_map(x0, mesh_ry, ny, H0)
  x -> VectorValue(
    x0 + x[1],
    map_vertical_GP_for_const_dep(x[2] - H0, mesh_ry, ny, H0; dbgmsg=false),
  )
end

"""
  incident_wave(; H0, ω, η0, α)

Helper function to build the incident wave state and associated functions for the inlet boundary condition.
- `H0` — still water depth
- `ω` — wave frequency
- `η0` — wave height (used to build the incident wave state)
- `α` — wave angle of incidence (0 = head-on)
Returns a named tuple with:
- `sea_state` — the `AiryState` object representing the incident wave state
- `ηin(x)` — function that returns the incident wave elevation at point `x`
- `ϕin(x)` — function that returns the incident wave potential at point `x`
- `vin(x)` — function that returns the incident wave velocity vector at point `x`
"""
function incident_wave(; H0::Real, ω::Real, η0::Real, α::Real)
  sea_state = build_regular_wave_state(H=2.0 * η0, T=2π / ω, h=H0, θ=α)
  wave(x) = WaveSpec.AiryWaves.generate_sea(sea_state, [x[1]], [0.], [x[2]], [0.0], vars=[:η, :ϕ, :u, :w])
  ηin(x) = wave(x)[:η][1]
  ϕin(x) = wave(x)[:ϕ][1]
  vin(x) = VectorValue(wave(x)[:u][1], wave(x)[:w][1])
  return (; sea_state, ηin, ϕin, vin)
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
- direct free-surface triangulation,
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
  Γfs = Γ
  Γκ = Γ

  degree = 2 * p.order
  dΩ = Measure(Ω, degree)
  dΓfs = Measure(Γfs, degree)
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
  k = inc.sea_state.k[1]

  a((ϕ, κ), (w, u)) =
    ∫(∇(w) ⋅ ∇(ϕ))dΩ +
    ∫(p.βₕ * (u + αₕ * w) * (WaveSpec.PhysicalConstants.g * κ - im * p.ω * ϕ) + im * p.ω * w * κ)dΓfs +
    ∫(-im * k * w * ϕ)dΓin +
    ∫(-im * k * w * ϕ)dΓout

  l((w, u)) =
    ∫(w * (inc.vin ⋅ nΓin))dΓin -
    ∫(im * k * w * inc.ϕin)dΓin

  op = AffineFEOperator(a, l, X, Y)
  ϕₕ, κₕ = solve(op)

  κ_from_ϕ = (im * p.ω / WaveSpec.PhysicalConstants.g * ϕₕ)(probes)
  κ_direct = κₕ(probes)

  filename = joinpath(PKG_ROOT, "data" , "VTK", "examples", "EmptyTankExample", "empty_tank_plain.vtu")
  isdir(dirname(filename)) || mkpath(dirname(filename))
  writevtk(Ω, filename, cellfields=["ϕ_re" => real(ϕₕ), "ϕ_im" => imag(ϕₕ)])

  return ( probe_surface=κₕ(probes), probe_potential=ϕₕ(probes) )
end

"""
    run_structured_implementation(; kwargs...)

HydroElasticFEM version of the same empty-tank problem.

Compared with `run_plain_implementation`, this version delegates:
- mesh partitioning to `Geometry`,
- field definitions to `Physics`,
- operator assembly and solve to `Simulation`.

The inlet forcing and radiation boundaries are attached to `PotentialFlow`.
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
  )

  sea_state = build_regular_wave_state(H=2.0 * p.η0, T=2π / p.ω, h=p.H0)
  f_in(x) = (inc.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc.sea_state.k[1] * inc.ϕin(x)

  p_flow = P.PotentialFlow(
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

  free_surface = P.FreeSurface(
    ρw=1025.0,
    g=WaveSpec.PhysicalConstants.g,
    βₕ=p.βₕ,
    fe=PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64}),
    space_domain_symbol=:Γκ,
  )

  config = S.FreqDomainConfig(ω=p.ω)
  problem = S.build_problem(tank, P.PhysicsParameters[p_flow, free_surface], config)
  result = S.simulate(problem)

  ϕₕ, κₕ = result.solution

  Ω = S.get_triangulations(problem)[:Ω]
  filename = joinpath(PKG_ROOT, "data" , "VTK", "examples", "EmptyTankExample", "empty_tank_structured.vtu")
  isdir(dirname(filename)) || mkpath(dirname(filename))
  writevtk(Ω, filename, cellfields=["ϕ_re" => real(ϕₕ), "ϕ_im" => imag(ϕₕ)])

  return ( probe_surface=κₕ(probes), probe_potential=ϕₕ(probes) )
  
end

"""
    compare_implementations(; atol=1e-8, rtol=1e-6, kwargs...)

Run both the plain and structured implementations and print a compact comparison
at the surface probes.
"""
function run_tutorial(; atol=1e-8, rtol=1e-6, kwargs...)
  plain = run_plain_implementation(; kwargs...)
  structured = run_structured_implementation(; kwargs...)

  plain_surface_vals = plain.probe_surface
  structured_surface_vals = structured.probe_surface
  err_surface = maximum(abs.(plain_surface_vals .- structured_surface_vals))

  plain_potential_vals = plain.probe_potential
  structured_potential_vals = structured.probe_potential
  err_potential = maximum(abs.(plain_potential_vals .- structured_potential_vals))  
  mismatch_count = 0

  @assert length(plain_surface_vals) == length(structured_surface_vals) "Probe arrays must have the same length."
  @assert length(plain_potential_vals) == length(structured_potential_vals) "Probe arrays must have the same length."
  
  for (κ_plain, κ_struct) in zip(plain_surface_vals, structured_surface_vals)
    if !isapprox(κ_plain, κ_struct; atol=atol, rtol=rtol)
      mismatch_count += 1
    end
  end

  for (ϕ_plain, ϕ_struct) in zip(plain_potential_vals, structured_potential_vals)
    if !isapprox(ϕ_plain, ϕ_struct; atol=atol, rtol=rtol)
      mismatch_count += 1
    end
  end

  if mismatch_count > 0
    @warn "Plain and structured implementations differ beyond tolerance." mismatch_count max_surface_probe_error=err_surface max_potential_probe_error=err_potential atol rtol
  end

  return nothing
end

end
