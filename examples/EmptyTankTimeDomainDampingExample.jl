module EmptyTankTimeDomainDampingExample

using Gridap
using Gridap.ODEs
using Parameters
using WaveSpec

using HydroElasticFEM: PKG_ROOT, map_vertical_GP_for_const_dep
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S

"""
    EmptyTankTimeDomainDampingParams

Parameters for the time-domain empty-tank example with damping zones.
"""
@with_kw struct EmptyTankTimeDomainDampingParams
  H0::Float64 = 10.0
  LΩ::Float64 = 60.0
  x0::Float64 = -10.0
  nx::Int = 240
  ny::Int = 12
  mesh_ry::Float64 = 1.08
  order::Int = 1
  βₕ::Float64 = 0.5
  Δt::Float64 = 0.02
  tf::Float64 = 6.0
  ρ∞::Float64 = 1.0
  η0::Float64 = 0.05
  T::Float64 = 3.0
  λ::Float64 = 20.0
  Ld::Float64 = 10.0
  inlet_amp::Float64 = 0.15
  μ₁_in::Float64 = 2.5
  μ₂_in::Float64 = 1.0
  μ₁_out::Float64 = 2.5
  μ₂_out::Float64 = 1.0
  probe_x::Vector{Float64} = collect(range(0.0, 40.0, length=7))
  write_vtk::Bool = true
end

function _validate_time_domain_damping_params(p::EmptyTankTimeDomainDampingParams)
  p.Δt > 0.0 || error("`Δt` must be positive.")
  p.tf > 0.0 || error("`tf` must be positive.")
  p.T > 0.0 || error("`T` must be positive.")
  p.λ > 0.0 || error("`λ` must be positive.")
  p.Ld > 0.0 || error("`Ld` must be positive.")
  2.0 * p.Ld < p.LΩ || error("Damping zones overlap. Require `2*Ld < LΩ`.")
  p.βₕ > 0.0 || error("`βₕ` must be positive.")
  return nothing
end

function shifted_gp_map(x0, mesh_ry, ny, H0)
  x -> VectorValue(
    x0 + x[1],
    map_vertical_GP_for_const_dep(x[2] - H0, mesh_ry, ny, H0; dbgmsg=false),
  )
end

function plain_gp_map(mesh_ry, ny, H0)
  x -> VectorValue(
    x[1],
    map_vertical_GP_for_const_dep(x[2], mesh_ry, ny, H0; dbgmsg=false),
  )
end

function probe_points(xs)
  Point.(xs, 0.0)
end

_with_write_vtk(kwargs::NamedTuple, flag::Bool) = merge(kwargs, (write_vtk=flag,))

function _max_history_error(a_hist, b_hist)
  length(a_hist) == length(b_hist) || error("History arrays must have the same length.")
  maximum(maximum(abs.(a .- b)) for (a, b) in zip(a_hist, b_hist))
end

"""
    run_example(; kwargs...)

HydroElasticFEM-only empty-tank time-domain simulation with two surface
damping zones.

Returns probe histories for free-surface elevation and velocity potential.
"""
function run_example(; kwargs...)
  p = EmptyTankTimeDomainDampingParams(; kwargs...)
  _validate_time_domain_damping_params(p)

  ω = 2π / p.T
  k = 2π / p.λ
  g = WaveSpec.PhysicalConstants.g
  x_d2_min = p.x0 + p.LΩ - p.Ld

  # Generalized-alpha choice used to derive a consistent free-surface stabilization.
  γₜ = 0.5
  βₜ = 0.25
  αₕ = γₜ / (βₜ * p.Δt) / g * (1.0 - p.βₕ) / p.βₕ

  tank = G.TankDomain2D(
    L=p.LΩ,
    H=p.H0,
    nx=p.nx,
    ny=p.ny,
    map=shifted_gp_map(p.x0, p.mesh_ry, p.ny, p.H0),
    damping_zones=[
      G.DampingZone1D(L=p.Ld, x₀=[p.x0, 0.0], domain_symbol=:Γd_in),
      G.DampingZone1D(L=p.Ld, x₀=[p.x0 + p.LΩ - p.Ld, 0.0], domain_symbol=:Γd_out),
    ],
  )

  η_in(x, t) = p.η0 * cos(k * x[1] - ω * t)
  vz_in(x, t) = p.η0 * ω * sin(k * x[1] - ω * t)
  inlet_v(x, t) = -(p.η0*ω)*(cosh(k*(p.H0 + x[2])) / sinh(k*p.H0))*cos(k*x[1]-ω*t)

  μ₁_in(x::VectorValue) = p.μ₁_in * (1.0 - sin(π / 2 * (x[1] - p.x0) / p.Ld)) * (x[1] <= p.x0 + p.Ld)
  μ₁_out(x::VectorValue) = p.μ₁_out * (1.0 - cos(π / 2 * (x[1] - x_d2_min) / p.Ld)) * (x[1] >= x_d2_min)
  μ₂_in(x) = μ₁_in(x) * k
  μ₂_out(x) = μ₁_out(x) * k

  p_flow = P.PotentialFlow(
    ρw=1025.0,
    g=g,
    boundary_conditions=[
      P.PrescribedInletPotentialBC(domain=:dΓin, forcing=(t -> (x -> inlet_v(x, t))), quantity=:traction),
      P.DampingZoneBC(
        domain=:dΓd_1,
        μ₁=μ₁_in,
        μ₂=μ₂_in,
        η_in=(t -> (x -> η_in(x, t))),
        vz_in=(t -> (x -> vz_in(x, t))),
      ),
      P.DampingZoneBC(
        domain=:dΓd_2,
        μ₁=μ₁_out,
        μ₂=μ₂_out,
        η_in=(x -> 0.0),
        vz_in=(x -> 0.0),
      ),
    ],
    fe=PH.FESpaceConfig(order=p.order),
    space_domain_symbol=:Ω,
  )

  free_surface = P.FreeSurface(
    ρw=1025.0,
    g=g,
    βₕ=p.βₕ,
    fe=PH.FESpaceConfig(order=p.order),
    space_domain_symbol=:Γκ,
  )

  config = S.TimeDomainConfig(t₀=0.0, tf=p.tf)
  tconfig = S.TimeConfig(
    Δt=p.Δt,
    t₀=0.0,
    tf=p.tf,
    ρ∞=p.ρ∞,
    αₕ=αₕ,
    u0=[0.0, 0.0],
    u0t=[0.0, 0.0],
    u0tt=[0.0, 0.0],
  )

  problem = S.build_problem(tank, P.PhysicsParameters[p_flow, free_surface], config; tconfig=tconfig)
  result = S.simulate(problem, tconfig)

  probes = probe_points(p.probe_x)
  t_hist = Float64[]
  κ_hist = Vector{Any}()
  ϕ_hist = Vector{Any}()

  pvd_Ω = nothing
  pvd_Γκ = nothing
  Ω = nothing
  Γκ = nothing
  folder = joinpath(PKG_ROOT, "data", "VTK", "examples", "EmptyTankExample")

  if p.write_vtk
    Ω = S.get_triangulations(problem)[:Ω]
    Γκ = S.get_triangulations(problem)[:Γκ]
    isdir(folder) || mkpath(folder)
    pvd_Ω = createpvd(joinpath(folder, "empty_tank_time_damping_fluid"))
    pvd_Γκ = createpvd(joinpath(folder, "empty_tank_time_damping_surface"))
  end

  for (t, uh) in result.solution
    ϕₕ, κₕ = uh
    push!(t_hist, t)
    push!(κ_hist, κₕ(probes))
    push!(ϕ_hist, ϕₕ(probes))

    if p.write_vtk
      tval = round(Int, t * 1000.0)
      pvd_Ω[t] = createvtk(
        Ω,
        joinpath(folder, "empty_tank_time_damping_fluid_$(tval).vtu"),
        cellfields=["re_ϕ" => real(ϕₕ), "im_ϕ" => imag(ϕₕ), "μ₁_in" => μ₁_in, "μ₂_in" => μ₂_in, "μ₁_out" => μ₁_out, "μ₂_out" => μ₂_out],
      )
      pvd_Γκ[t] = createvtk(
        Γκ,
        joinpath(folder, "empty_tank_time_damping_surface_$(tval).vtu"),
        cellfields=["re_κ" => real(κₕ), "im_κ" => imag(κₕ)],
      )
    end
  end

  if p.write_vtk
    savepvd(pvd_Ω)
    savepvd(pvd_Γκ)
  end

  return (
    times=t_hist,
    probe_surface=κ_hist,
    probe_potential=ϕ_hist,
    αₕ=αₕ,
    problem=problem,
    result=result,
  )
end

"""
    run_plain_example(; kwargs...)

Plain Gridap implementation of the same empty-tank time-domain damping-zone
problem solved by `run_example`.
"""
function run_plain_example(; kwargs...)
  p = EmptyTankTimeDomainDampingParams(; kwargs...)
  _validate_time_domain_damping_params(p)

  ω = 2π / p.T
  k = 2π / p.λ
  g = WaveSpec.PhysicalConstants.g

  γₜ = 0.5
  βₜ = 0.25
  αₕ = γₜ / (βₜ * p.Δt) / g * (1.0 - p.βₕ) / p.βₕ

  domain = (p.x0, p.x0 + p.LΩ, -p.H0, 0.0)
  partition = (p.nx, p.ny)
  model = CartesianDiscreteModel(domain, partition, map=plain_gp_map(p.mesh_ry, p.ny, p.H0))

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "surface", [3, 4, 6])
  add_tag_from_tags!(labels, "bottom", [1, 2, 5])
  add_tag_from_tags!(labels, "inlet", [7])
  add_tag_from_tags!(labels, "outlet", [8])
  add_tag_from_tags!(labels, "water", [9])

  Ω = Interior(model)
  Γ = Boundary(model, tags="surface")
  Γin = Boundary(model, tags="inlet")

  xΓ = get_cell_coordinates(Γ)
  x_d1_min = p.x0
  x_d1_max = p.x0 + p.Ld
  x_d2_min = p.x0 + p.LΩ - p.Ld

  is_damping1(xs) = begin
    c = (1 / length(xs)) * sum(xs)
    (x_d1_min <= c[1] <= x_d1_max) && (c[2] ≈ 0.0)
  end
  is_damping2(xs) = begin
    c = (1 / length(xs)) * sum(xs)
    (x_d2_min <= c[1]) && (c[2] ≈ 0.0)
  end

  Γd1_to_Γ_mask = lazy_map(is_damping1, xΓ)
  Γd2_to_Γ_mask = lazy_map(is_damping2, xΓ)
  Γd1 = Triangulation(Γ, findall(Γd1_to_Γ_mask))
  Γd2 = Triangulation(Γ, findall(Γd2_to_Γ_mask))
  Γfs = Triangulation(Γ, findall(!, Γd1_to_Γ_mask .| Γd2_to_Γ_mask))
  Γκ = Γ

  degree = 2 * p.order
  dΩ = Measure(Ω, degree)
  dΓfs = Measure(Γfs, degree)
  dΓd1 = Measure(Γd1, degree)
  dΓd2 = Measure(Γd2, degree)
  dΓin = Measure(Γin, degree)
  nΓd1 = get_normal_vector(Γd1)
  nΓd2 = get_normal_vector(Γd2)

  reffe = ReferenceFE(lagrangian, Float64, p.order)
  VΩ = TestFESpace(Ω, reffe; conformity=:H1)
  VΓκ = TestFESpace(Γκ, reffe; conformity=:H1)
  UΩ = TransientTrialFESpace(VΩ)
  UΓκ = TransientTrialFESpace(VΓκ)
  X = TransientMultiFieldFESpace([UΩ, UΓκ])
  Y = MultiFieldFESpace([VΩ, VΓκ])

  η_in(x, t) = p.η0 * cos(k * x[1] - ω * t)
  vz_in(x, t) = p.η0 * ω * sin(k * x[1] - ω * t)
  inlet_v(x, t) = -(p.η0*ω)*(cosh(k*(p.H0 + x[2])) / sinh(k*p.H0))*cos(k*x[1]-ω*t)

  # Follow mem_time_damp_GridapUpdate.jl damping-zone profiles.
  μ₀_in = p.μ₁_in
  μ₀_out = p.μ₁_out
  μ₁_in(x::VectorValue) = μ₀_in * (1.0 - sin(π / 2 * (x[1] - p.x0) / p.Ld)) * (x[1] <= p.x0 + p.Ld)
  μ₁_out(x::VectorValue) = μ₀_out * (1.0 - cos(π / 2 * (x[1] - x_d2_min) / p.Ld)) * (x[1] >= x_d2_min)
  μ₂_in(x) = μ₁_in(x) * k
  μ₂_out(x) = μ₁_out(x) * k
  ηd(t) = x -> μ₂_in(x) * η_in(x, t)
  ∇ₙϕd(t) = x -> μ₁_in(x) * vz_in(x, t)

  ∇ₙd1(ϕ) = ∇(ϕ) ⋅ nΓd1
  ∇ₙd2(ϕ) = ∇(ϕ) ⋅ nΓd2

  m(t, (ϕₜₜ, κₜₜ), (w, u)) = ∫(0.0 * ϕₜₜ * w)dΩ
  m(t, u, ∂ₜₜu, v) = m(t, ∂ₜₜu, v)

  c(t, (ϕₜ, κₜ), (w, u)) =
    ∫(p.βₕ * (u + αₕ * w) * ϕₜ - w * κₜ)dΓfs +
    ∫(p.βₕ * (u + αₕ * w) * ϕₜ - w * κₜ)dΓd1 +
    ∫(p.βₕ * (u + αₕ * w) * ϕₜ - w * κₜ)dΓd2

  a(t, (ϕ, κ), (w, u)) =
    ∫(∇(w) ⋅ ∇(ϕ))dΩ +
    ∫(p.βₕ * (u + αₕ * w) * g * κ)dΓfs +
    ∫(p.βₕ * (u + αₕ * w) * g * κ - μ₂_in * κ * w + μ₁_in * ∇ₙd1(ϕ) * (u + αₕ * w))dΓd1 +
    ∫(p.βₕ * (u + αₕ * w) * g * κ - μ₂_out * κ * w + μ₁_out * ∇ₙd2(ϕ) * (u + αₕ * w))dΓd2

  l(t, (w, u)) =
    ∫(w * (x -> inlet_v(x, t)))dΓin -
    ∫(ηd(t) * w - ∇ₙϕd(t) * (u + αₕ * w))dΓd1

  op = TransientLinearFEOperator((a, c, m), l, X, Y; constant_forms=(true, true, true))
  ls = LUSolver()
  ode_solver = GeneralizedAlpha2(ls, p.Δt, p.ρ∞)

  u0 = interpolate_everywhere([0.0, 0.0], X(0.0))
  u0t = interpolate_everywhere([0.0, 0.0], X(0.0))
  u0tt = interpolate_everywhere([0.0, 0.0], X(0.0))
  solution = solve(ode_solver, op, 0.0, p.tf, (u0, u0t, u0tt))

  probes = probe_points(p.probe_x)
  t_hist = Float64[]
  κ_hist = Vector{Any}()
  ϕ_hist = Vector{Any}()

  pvd_Ω = nothing
  pvd_Γκ = nothing
  folder = joinpath(PKG_ROOT, "data", "VTK", "examples", "EmptyTankExample")

  if p.write_vtk
    isdir(folder) || mkpath(folder)
    pvd_Ω = createpvd(joinpath(folder, "empty_tank_time_damping_plain_fluid"))
    pvd_Γκ = createpvd(joinpath(folder, "empty_tank_time_damping_plain_surface"))
  end

  for (t, uh) in solution
    ϕₕ, κₕ = uh
    push!(t_hist, t)
    push!(κ_hist, κₕ(probes))
    push!(ϕ_hist, ϕₕ(probes))

    if p.write_vtk
      tval = round(Int, t * 1000.0)
      pvd_Ω[t] = createvtk(
        Ω,
        joinpath(folder, "empty_tank_time_damping_plain_fluid_$(tval).vtu"),
        cellfields=["re_ϕ" => real(ϕₕ), "im_ϕ" => imag(ϕₕ), "μ₁_in" => μ₁_in, "μ₂_in" => μ₂_in, "μ₁_out" => μ₁_out, "μ₂_out" => μ₂_out],
      )
      pvd_Γκ[t] = createvtk(
        Γκ,
        joinpath(folder, "empty_tank_time_damping_plain_surface_$(tval).vtu"),
        cellfields=["re_κ" => real(κₕ), "im_κ" => imag(κₕ)],
      )
    end
  end

  if p.write_vtk
    savepvd(pvd_Ω)
    savepvd(pvd_Γκ)
  end

  return (
    times=t_hist,
    probe_surface=κ_hist,
    probe_potential=ϕ_hist,
    αₕ=αₕ,
    model=model,
    solution=solution,
  )
end

"""
  compare_implementations(; atol=1e-4, rtol=1e-2, kwargs...)

Run both implementations (plain Gridap and HydroElasticFEM modules) for the
same setup and return probe-history mismatch metrics.
"""
function compare_implementations(; atol=1e-6, rtol=1e-8, kwargs...)
  kw = _with_write_vtk((; kwargs...), false)

  plain = run_plain_example(; kw...)
  structured = run_example(; kw...)

  length(plain.times) == length(structured.times) || error("Time grids differ between implementations.")
  time_err = maximum(abs.(plain.times .- structured.times))
  κ_err = _max_history_error(plain.probe_surface, structured.probe_surface)
  ϕ_err = _max_history_error(plain.probe_potential, structured.probe_potential)

  κ_scale = maximum(maximum(abs.(v)) for v in structured.probe_surface)
  ϕ_scale = maximum(maximum(abs.(v)) for v in structured.probe_potential)
  κ_ok = κ_err <= (atol + rtol * κ_scale)
  ϕ_ok = ϕ_err <= (atol + rtol * ϕ_scale)

  println("Time error: $(time_err)")
  println("Surface probe error: $(κ_err) (scale: $(κ_scale), within tolerance: $(κ_ok))")
  println("Potential probe error: $(ϕ_err) (scale: $(ϕ_scale), within tolerance: $(ϕ_ok))")

  return (
    max_time_error=time_err,
    max_surface_probe_error=κ_err,
    max_potential_probe_error=ϕ_err
  )
end

end
