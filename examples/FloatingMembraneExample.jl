module FloatingMembraneExample

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
Tutorial: Floating membrane in the frequency domain
===================================================

This example implements the plain Gridap formulation for a floating membrane in
potential flow and then rewrites the same problem using the `HydroElasticFEM` 
package structure.

The model is a 2D tank with:
- a floating membrane on part of the free surface,
- a free-surface elevation field `κ` on the uncovered surface,
- a membrane displacement field `η` on the structure,
- radiation conditions on inlet and outlet,
- monochromatic incident-wave forcing at the inlet.

The two functions below solve the same monochromatic test case:
- `run_plain_implementation()` reproduces the explicit Gridap formulation.
- `run_structured_implementation()` uses the `Geometry`, `Physics`, and
  `Simulation` modules from `HydroElasticFEM`.
"""

@with_kw struct FloatingMembraneTutorialParams
  H0::Float64 = 10.0
  LΩ::Float64 = 60.0
  x0::Float64 = 0.0
  ω::Float64 = 2.0
  η0::Float64 = 0.25
  α::Float64 = 0.0
  βₕ::Float64 = 0.5
  order::Int = 1
  nx::Int = 240
  ny::Int = 12
  mesh_ry::Float64 = 1.08
  Lm::Float64 = 20.0
  xm0::Float64 = 20.0
  mᵨ::Float64 = 0.9
  Tᵨ::Float64 = 98.1
  τ::Float64 = 0.0
  mem_bnd_type::Symbol = :free
  probe_x::Vector{Float64} = collect(range(5.0, 55.0, length=9))
end

"""
    build_regular_wave_state(; H, T, h, θ=0.0)

Build a monochromatic regular-wave state using `WaveSpec`.
- `H` — wave height
- `T` — wave period
- `h` — water depth
- `θ` — incident-wave angle
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
    gp_map(mesh_ry, ny, H0)

Return the vertical mesh-grading map for the plain Cartesian tank model.
- `mesh_ry` — grading ratio
- `ny` — number of vertical elements
- `H0` — still-water depth
"""
function gp_map(mesh_ry, ny, H0)
  x -> VectorValue(x[1], map_vertical_GP_for_const_dep(x[2], mesh_ry, ny, H0; dbgmsg=false))
end

"""
    shifted_gp_map(x0, H0, mesh_ry, ny)

Return the vertical mesh-grading map for the structured tank, including the
horizontal shift `x0` and the `[0, H0]` vertical convention used by
`TankDomain2D`.
"""
function shifted_gp_map(x0, H0, mesh_ry, ny)
  x -> VectorValue(
    x0 + x[1],
    map_vertical_GP_for_const_dep(x[2] - H0, mesh_ry, ny, H0; dbgmsg=false),
  )
end

"""
    incident_wave(; H0, ω, η0, α)

Build the incident monochromatic wave and return helper functions for surface
elevation, potential, and velocity evaluation.
"""
function incident_wave(; H0::Real, ω::Real, η0::Real, α::Real)
  sea_state = build_regular_wave_state(H=2.0 * η0, T=2π / ω, h=H0, θ=α)
  wave(x) = WaveSpec.AiryWaves.generate_sea(sea_state, [x[1]], [0.0], [x[2]], [0.0], vars=[:η, :ϕ, :u, :w])
  ηin(x) = wave(x)[:η][1]
  ϕin(x) = wave(x)[:ϕ][1]
  vin(x) = VectorValue(wave(x)[:u][1], wave(x)[:w][1])
  return (; sea_state, ηin, ϕin, vin)
end

"""
    probe_points(xs)

Convert a vector of x-coordinates into free-surface probe points `(x, 0.0)`.
"""
function probe_points(xs)
  Point.(xs, 0.0)
end

"""
    membrane_indicator(xs, xm0, xm1)

Return a boolean mask identifying which probe coordinates lie on the membrane
span `[xm0, xm1]`.
"""
function membrane_indicator(xs, xm0, xm1)
  (xm0 .<= xs) .& (xs .<= xm1)
end

"""
    validate_params(p)

Validate that the membrane configuration encoded in
`FloatingMembraneTutorialParams` is physically contained inside the tank and
that the membrane boundary type is supported.
"""
function validate_params(p::FloatingMembraneTutorialParams)
  p.mem_bnd_type in (:free, :fixed) || error("`mem_bnd_type` must be either `:free` or `:fixed`.")
  xm1 = p.xm0 + p.Lm
  p.x0 <= p.xm0 || error("Membrane start must lie inside the tank.")
  xm1 <= p.x0 + p.LΩ || error("Membrane end must lie inside the tank.")
  return nothing
end

"""
    add_membrane_endpoint_tag!(model, Γη)

Create the `"mem_bnd"` face tag at the endpoints of the membrane
triangulation `Γη` so fixed-end membrane FE spaces can attach Dirichlet
conditions there.
"""
function add_membrane_endpoint_tag!(model, Γη)
  labels = get_face_labeling(model)
  Λmb = Boundary(Γη)
  xΛmb = get_cell_coordinates(Λmb)
  isempty(xΛmb) && return Λmb

  vertex_coords = model.grid_topology.vertex_coordinates
  vertex_ids = Int[]
  for x in xΛmb
    matches = findall(vertex_coords .== x)
    isempty(matches) && error("Could not locate membrane endpoint vertex in model topology.")
    push!(vertex_ids, matches[1])
  end

  new_entity = num_entities(labels) + 1
  for vid in unique(vertex_ids)
    labels.d_to_dface_to_entity[1][vid] = new_entity
  end
  add_tag!(labels, "mem_bnd", [new_entity])
  return Λmb
end

"""
    write_membrane_vtk(dirname_suffix, Ω, Γκ, Γη, ϕₕ, κₕ, ηₕ)

Write VTK output for the fluid field, uncovered free surface, and membrane
response into the example output folder identified by `dirname_suffix`.
"""
function write_membrane_vtk(dirname_suffix, Ω, Γκ, Γη, ϕₕ, κₕ, ηₕ)
  folder = joinpath(PKG_ROOT, "data", "VTK", "examples", "FloatingMembraneExample", dirname_suffix)
  isdir(folder) || mkpath(folder)

  vx = ∇(ϕₕ) ⋅ VectorValue(1.0, 0.0)
  vy = ∇(ϕₕ) ⋅ VectorValue(0.0, 1.0)
  ηx = ∇(ηₕ) ⋅ VectorValue(1.0, 0.0)

  writevtk(
    Ω,
    joinpath(folder, "floating_membrane_fluid"),
    cellfields=[
      "ϕ_re" => real(ϕₕ),
      "ϕ_im" => imag(ϕₕ),
      "vx_re" => real(vx),
      "vx_im" => imag(vx),
      "vy_re" => real(vy),
      "vy_im" => imag(vy),
    ],
  )
  writevtk(
    Γκ,
    joinpath(folder, "floating_membrane_kappa"),
    cellfields=["κ_re" => real(κₕ), "κ_im" => imag(κₕ)],
  )
  writevtk(
    Γη,
    joinpath(folder, "floating_membrane_eta"),
    cellfields=["η_re" => real(ηₕ), "η_im" => imag(ηₕ), "ηx_re" => real(ηx), "ηx_im" => imag(ηx)],
  )
end

"""
    run_plain_implementation(; kwargs...)

Solve the floating-membrane problem with the explicit Gridap formulation that
matches `src/Membrane2D/FrequencyDomain/mem_freq_rad_fnc.jl`.

Keyword arguments override `FloatingMembraneTutorialParams`.
Returns a named tuple with probe values for the combined surface response and
the velocity potential.
"""
function run_plain_implementation(; kwargs...)
  p = FloatingMembraneTutorialParams(; kwargs...)
  validate_params(p)

  inc = incident_wave(; H0=p.H0, ω=p.ω, η0=p.η0, α=p.α)
  probes = probe_points(p.probe_x)
  xm1 = p.xm0 + p.Lm
  probe_membrane = membrane_indicator(p.probe_x, p.xm0, xm1)

  model = CartesianDiscreteModel(
    (p.x0, p.x0 + p.LΩ, -p.H0, 0.0),
    (p.nx, p.ny),
    map=gp_map(p.mesh_ry, p.ny, p.H0),
  )

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

  xΓ = get_cell_coordinates(Γ)
  Γm_to_Γ_mask = lazy_map(xs -> begin
    n = length(xs)
    x = (1 / n) * sum(xs)
    (p.xm0 <= x[1] <= xm1) && (x[2] ≈ 0.0)
  end, xΓ)

  membrane_ids = findall(Γm_to_Γ_mask)
  isempty(membrane_ids) && error("The membrane region is not resolved by the current mesh. Increase `nx` or adjust `xm0`/`Lm`.")
  Γη = Triangulation(Γ, membrane_ids)
  Γκ = Triangulation(Γ, findall(!, Γm_to_Γ_mask))
  Λmb = add_membrane_endpoint_tag!(model, Γη)

  degree = 2 * p.order
  dΩ = Measure(Ω, degree)
  dΓm = Measure(Γη, degree)
  dΓfs = Measure(Γκ, degree)
  dΓin = Measure(Γin, degree)
  dΓout = Measure(Γout, degree)
  dΛmb = Measure(Λmb, degree)
  nΛmb = get_normal_vector(Λmb)
  nΓin = get_normal_vector(Γin)

  reffe = ReferenceFE(lagrangian, Float64, p.order)
  VΩ = TestFESpace(Ω, reffe; conformity=:H1, vector_type=Vector{ComplexF64})
  VΓκ = TestFESpace(Γκ, reffe; conformity=:H1, vector_type=Vector{ComplexF64})
  VΓη = if p.mem_bnd_type == :fixed
    TestFESpace(Γη, reffe; conformity=:H1, vector_type=Vector{ComplexF64}, dirichlet_tags=["mem_bnd"])
  else
    TestFESpace(Γη, reffe; conformity=:H1, vector_type=Vector{ComplexF64})
  end

  UΩ = TrialFESpace(VΩ)
  UΓκ = TrialFESpace(VΓκ)
  UΓη = if p.mem_bnd_type == :fixed
    TrialFESpace(VΓη, x -> ComplexF64(0.0))
  else
    TrialFESpace(VΓη)
  end

  X = MultiFieldFESpace([UΩ, UΓκ, UΓη])
  Y = MultiFieldFESpace([VΩ, VΓκ, VΓη])

  αₕ = -im * p.ω / WaveSpec.PhysicalConstants.g * (1.0 - p.βₕ) / p.βₕ
  k = inc.sea_state.k[1]

  function res_membrane((ϕ, κ, η), (w, u, v))
    common =
      ∫(v * (WaveSpec.PhysicalConstants.g * η - im * p.ω * ϕ) + im * p.ω * w * η - p.mᵨ * v * p.ω^2 * η + p.Tᵨ * (1 - im * p.ω * p.τ) * ∇(v) ⋅ ∇(η))dΓm
    if p.mem_bnd_type == :fixed
      return common + ∫(-p.Tᵨ * (1 - im * p.ω * p.τ) * v * ∇(η) ⋅ nΛmb)dΛmb
    end
    return common
  end

  a((ϕ, κ, η), (w, u, v)) =
    res_membrane((ϕ, κ, η), (w, u, v)) +
    ∫(∇(w) ⋅ ∇(ϕ))dΩ +
    ∫(p.βₕ * (u + αₕ * w) * (WaveSpec.PhysicalConstants.g * κ - im * p.ω * ϕ) + im * p.ω * w * κ)dΓfs +
    ∫(-im * k * w * ϕ)dΓin +
    ∫(-im * k * w * ϕ)dΓout

  l((w, u, v)) =
    ∫(w * (inc.vin ⋅ nΓin))dΓin -
    ∫(im * k * w * inc.ϕin)dΓin

  op = AffineFEOperator(a, l, X, Y)
  ϕₕ, κₕ, ηₕ = solve(op)

  probe_surface = similar(ϕₕ(probes))
  probe_surface[.!probe_membrane] = κₕ(probes[.!probe_membrane])
  probe_surface[probe_membrane] = ηₕ(probes[probe_membrane])

  write_membrane_vtk("plain", Ω, Γκ, Γη, ϕₕ, κₕ, ηₕ)

  return (
    probe_surface=probe_surface,
    probe_potential=ϕₕ(probes),
  )
end

"""
    run_structured_implementation(; kwargs...)

Solve the same floating-membrane problem using the structured
`HydroElasticFEM` geometry, physics, and simulation APIs.

Keyword arguments override `FloatingMembraneTutorialParams`.
Returns a named tuple with probe values for the combined surface response and
the velocity potential.
"""
function run_structured_implementation(; kwargs...)
  p = FloatingMembraneTutorialParams(; kwargs...)
  validate_params(p)

  probes = probe_points(p.probe_x)
  xm1 = p.xm0 + p.Lm
  probe_membrane = membrane_indicator(p.probe_x, p.xm0, xm1)

  tank = G.TankDomain2D(
    L=p.LΩ,
    H=p.H0,
    nx=p.nx,
    ny=p.ny,
    map=shifted_gp_map(p.x0, p.H0, p.mesh_ry, p.ny),
    structure_domains=[
      G.StructureDomain1D(L=p.Lm, x₀=[p.xm0, 0.0], domain_symbol=:Γm),
    ],
  )

  inc = incident_wave(; H0=p.H0, ω=p.ω, η0=p.η0, α=p.α)
  f_in(x) = (inc.vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * inc.sea_state.k[1] * inc.ϕin(x)

  potential = P.PotentialFlow(
    ρw=1025.0,
    g=WaveSpec.PhysicalConstants.g,
    sea_state=inc.sea_state,
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

  mem_fe = if p.mem_bnd_type == :fixed
    PH.FESpaceConfig(
      order=p.order,
      vector_type=Vector{ComplexF64},
      dirichlet_tags=["mem_bnd"],
      dirichlet_value=x -> ComplexF64(0.0),
    )
    @error "Fixed membrane boundary conditions are not yet implemented in the structured API. Please set `mem_bnd_type` to `:free`."
  else
    PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64})
  end

  membrane = P.Membrane2D(
    L=p.Lm,
    mᵨ=p.mᵨ,
    Tᵨ=p.Tᵨ,
    τ=p.τ,
    g=WaveSpec.PhysicalConstants.g,
    fe=mem_fe,
    space_domain_symbol=:Γη,
  )
  
  physics = P.PhysicsParameters[potential, free_surface, membrane]
  config = S.FreqDomainConfig(ω=p.ω)

  problem = S.build_problem(tank, physics, config)

  result = S.simulate(problem)
  ϕₕ, κₕ, ηₕ = result.solution

  trians = S.get_triangulations(problem)
  probe_surface = similar(ϕₕ(probes))
  probe_surface[.!probe_membrane] = κₕ(probes[.!probe_membrane])
  probe_surface[probe_membrane] = ηₕ(probes[probe_membrane])

  write_membrane_vtk("structured", trians[:Ω], trians[:Γκ], trians[:Γη], ϕₕ, κₕ, ηₕ)

  return (
    probe_surface=probe_surface,
    probe_potential=ϕₕ(probes)
  )
end

"""
    run_tutorial(; atol=1e-8, rtol=1e-6, kwargs...)

Run both the plain and structured floating-membrane implementations and report
their probe-level agreement.

Keyword arguments are forwarded to both implementations through
`FloatingMembraneTutorialParams`.
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

  println("Maximum surface probe error: ", err_surface)
  println("Maximum potential probe error: ", err_potential)

  mismatch_count = 0
  @assert length(plain_surface_vals) == length(structured_surface_vals) "Probe arrays must have the same length."
  @assert length(plain_potential_vals) == length(structured_potential_vals) "Probe arrays must have the same length."

  for (η_plain, η_struct) in zip(plain_surface_vals, structured_surface_vals)
    if !isapprox(η_plain, η_struct; atol=atol, rtol=rtol)
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
