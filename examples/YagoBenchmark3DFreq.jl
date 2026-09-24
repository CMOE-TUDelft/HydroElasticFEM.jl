# ============================================================
# Yago 3D VLFS Benchmark — Frequency Domain
# Reproduces Section 5.5 of Colomes et al. (2022)
# DOI: 10.1002/nme.7140
#
# Physical setup:
#   3D potential flow domain (x=propagation, y=transverse, z=depth)
#   Domain: x in [0,10L], y in [-9B,9B], z in [0,H]
#   Floating plate: x in [4.5L,5.5L], y in [-B/2,B/2], z=H
#   Regular waves from inlet (x=0), sponge layers at inlet/outlet
#
# Execution order:
#   1. Warm-up (nx=2, ny=2, nz=1, order=2)
#   2. Case 1  (nx=32, ny=4, nz=3, order=4, lambda/L=0.4)
#
# Output:
#   VTK files in examples/output/
#   Plot: |eta|/eta0 versus xi in [-1,1]
# ============================================================

module YagoBenchmark3DFreq

using Gridap
using Parameters
using Printf
using Plots

import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P

export run_yago_3d_freq, run_yago_case1, yago_tank, solve_yago

@with_kw struct YagoCaseParams
  nx::Int = 32
  ny::Int = 4
  nz::Int = 3
  order::Int = 4
  λfactor::Float64 = 0.4
  dfactor::Float64 = 4.0
  vtk_output::Bool = true
  verbose::Bool = true
end

function _constants()
  L = 300.0
  B = 60.0
  H = 58.5
  hb = 2.0
  nLΩ = 10
  nBΩ = 18
  LΩ = nLΩ * L
  BΩ = nBΩ * B
  xb₀ = 4.5 * L
  xb₁ = xb₀ + L
  yb₀ = -B / 2
  yb₁ = yb₀ + B

  g = 9.81
  ρ = 1025.0
  ρb = 256.25
  E = 11.9e9
  ν = 0.13
  D = E * (hb^3 / 12) / (1 - ν^2)
  d₀ = ρb * hb / ρ
  C = P.build_KL_tensor(E, ν, hb, ρ)

  (; L, B, H, hb, nLΩ, nBΩ, LΩ, BΩ, xb₀, xb₁, yb₀, yb₁,
   g, ρ, ρb, E, ν, D, d₀, C)
end

function _wave(c, λfactor)
  λ = c.L * λfactor
  k = 2π / λ
  ω = sqrt(c.g * k * tanh(k * c.H))

  η₀ = 0.01
  ηᵢₙ(x) = η₀ * exp(im * k * x[1])
  ϕᵢₙ(x) = -im * (η₀ * ω / k) * (cosh(k * x[3]) / sinh(k * c.H)) *
            exp(im * k * x[1])
  vᵢₙ(x) = (η₀ * ω) * (cosh(k * x[3]) / sinh(k * c.H)) *
           exp(im * k * x[1])
  vzᵢₙ(x) = -im * ω * η₀ * exp(im * k * x[1])

  (; λ, k, ω, η₀, ηᵢₙ, ϕᵢₙ, vᵢₙ, vzᵢₙ)
end

function _damping(c, wave, dfactor)
  μ₀ = 2.5
  Ld = dfactor * c.L
  Ld₀ = dfactor * c.L
  xd = c.LΩ - Ld
  xd₀ = Ld₀

  μ₁(x) =
    μ₀ * (1 - cos(π / 2 * (x[1] - xd) / Ld)) * (x[1] > xd) +
    μ₀ * (1 - cos(π / 2 * (Ld₀ - x[1]) / Ld₀)) * (x[1] < xd₀)
  μ₂(x) = μ₁(x) * wave.k
  ηd(x) = μ₂(x) * wave.ηᵢₙ(x) * (x[1] < xd₀)
  ∇ₙϕd(x) = μ₁(x) * wave.vzᵢₙ(x) * (x[1] < xd₀)

  (; μ₁, μ₂, ηd, ∇ₙϕd)
end

"""
    yago_tank(c, nx, ny, nz, dfactor) -> TankDomain{3}

Geometry of the Yago benchmark: graded 3D tank, floating plate on the free
surface (`z = H`) and full-width inlet/outlet damping zones of length
`dfactor * L`.  All triangulations and measures derive from this object.
"""
function yago_tank(c, nx, ny, nz, dfactor)
  Ld = dfactor * c.L
  cartesian = G.CartesianDomain(
    mins = (0.0, -c.BΩ / 2, 0.0),
    maxs = (c.LΩ, c.BΩ / 2, c.H),
    parts = (c.nLΩ * nx, c.nBΩ * ny, nz),
    map = x -> G.map_fn(x, c.H, nz; grading_base = 2.5),
  )
  G.TankDomain(cartesian;
    structure_domains = [
      G.StructureDomain(L = c.L, W = c.B, x₀ = [c.xb₀, c.yb₀, c.H], domain_symbol = :Γ_plate),
    ],
    damping_zones = [
      G.DampingZone(L = Ld, x₀ = [0.0, -c.BΩ / 2, c.H], domain_symbol = :Γ_d_in),
      G.DampingZone(L = Ld, x₀ = [c.LΩ - Ld, -c.BΩ / 2, c.H], domain_symbol = :Γ_d_out),
    ],
  )
end

function run_yago_3d_freq(; nx=32, ny=4, nz=3, order=4, λfactor=0.4,
                          dfactor=4.0, vtk_output=true, verbose=true)
  c = _constants()
  wave = _wave(c, λfactor)
  damping = _damping(c, wave, dfactor)

  # Geometry: every domain and measure comes from the TankDomain
  #   Ω   fluid volume           Γκ  free surface incl. damping zones (κ field)
  #   Γin inlet                  Γη  plate (η field)
  #   Λη  plate interior skeleton (C/DG terms)
  tank = yago_tank(c, nx, ny, nz, dfactor)
  model = G.build_model(tank)
  trians = G.build_triangulations(tank, model)
  dom = G.get_integration_domains(trians; degree = 2 * order)
  geo = (Ω = trians[:Ω], Γf = trians[:Γκ], Γb = trians[:Γη],
         dΩ = dom[:dΩ], dΓᵢₙ = dom[:dΓin], dΓf = dom[:dΓκ], dΓb = dom[:dΓη],
         dΛb = dom[:dΛη], nΛb = dom[:n_Λ_η])

  return solve_yago(c, wave, damping, geo; nx, order, λfactor, vtk_output, verbose)
end

"""
    solve_yago(c, wave, damping, geo; nx, order, λfactor, vtk_output, verbose)

Assemble and solve the Yago frequency-domain problem on the geometry `geo`
(triangulations `Ω, Γf, Γb` and measures `dΩ, dΓᵢₙ, dΓf, dΓb, dΛb`, skeleton
normal `nΛb`).  Returns `(ξ, |η|/η₀)` along the plate centreline.
"""
function solve_yago(c, wave, damping, geo; nx, order, λfactor, vtk_output, verbose)
  (; Ω, Γf, Γb, dΩ, dΓᵢₙ, dΓf, dΓb, dΛb, nΛb) = geo

  βₕ = 0.5
  αₕ = -im * wave.ω / c.g * (1 - βₕ) / βₕ

  # Penalty length: Δx, as in Colomes et al. (2022).  Note dom[:h_η] is the
  # isotropic sqrt(Δx·Δy), which differs for the non-square cells used here.
  h = c.LΩ / (c.nLΩ * nx)
  # Yago 3D plate uses order*(order+1)/h (Liu 2D uses order*(order-1)/h).
  γ = 1.0 * order * (order + 1) / h

  reffeη = ReferenceFE(lagrangian, Float64, order)
  reffeκ = ReferenceFE(lagrangian, Float64, order)
  reffeᵩ = ReferenceFE(lagrangian, Float64, 2)

  Vϕ = TestFESpace(Ω, reffeᵩ; conformity=:H1, vector_type=Vector{ComplexF64})
  Uϕ = TrialFESpace(Vϕ)

  Vκ = TestFESpace(Γf, reffeκ; conformity=:H1, vector_type=Vector{ComplexF64})
  Uκ = TrialFESpace(Vκ)

  Vη = TestFESpace(Γb, reffeη; conformity=:H1, vector_type=Vector{ComplexF64})
  Uη = TrialFESpace(Vη)

  Y = MultiFieldFESpace([Vϕ, Vκ, Vη])
  X = MultiFieldFESpace([Uϕ, Uκ, Uη])

  ∇ₙ(ϕ) = ∇(ϕ) ⋅ VectorValue(0.0, 0.0, 1.0)

  a((ϕ, κ, η), (w, u, v)) =
    ∫(∇(ϕ) ⋅ ∇(w))dΩ +
    ∫(
      βₕ * (c.g * κ - im * wave.ω * ϕ) * (u + αₕ * w) +
      im * wave.ω * κ * w - damping.μ₂ * κ * w +
      damping.μ₁ * ∇ₙ(ϕ) * (u + αₕ * w)
    )dΓf +
    ∫(
      ((c.g - wave.ω^2 * c.d₀) * η - im * wave.ω * ϕ) * v +
      im * wave.ω * η * w +
      (∇∇(v) ⊙ (c.C ⊙ ∇∇(η)))
    )dΓb +
    ∫(
      -jump(∇(v)) ⊙ (mean(c.C ⊙ ∇∇(η)) ⋅ nΛb.⁺) -
      (mean(c.C ⊙ ∇∇(v)) ⋅ nΛb.⁺) ⊙ jump(∇(η)) +
      c.D / c.ρ * γ * jump(∇(v)) ⊙ jump(∇(η))
    )dΛb

  b((w, u, v)) =
    ∫(wave.vᵢₙ * w)dΓᵢₙ -
    ∫(damping.ηd * w - damping.∇ₙϕd * (u + αₕ * w))dΓf

  op = AffineFEOperator(a, b, X, Y)
  sol = solve(LUSolver(), op)
  ϕₕ, κₕ, ηₕ = sol

  if vtk_output
    outdir = joinpath(@__DIR__, "output")
    mkpath(outdir)
    filename = joinpath(outdir, "Yago3D_lambda04")
    writevtk(Ω, filename * "_phi",
             cellfields=["phi_re" => real(ϕₕ), "phi_im" => imag(ϕₕ)],
             nsubcells=10)
    writevtk(Γf, filename * "_kappa",
             cellfields=["kappa_re" => real(κₕ), "kappa_im" => imag(κₕ)],
             nsubcells=10)
    writevtk(Γb, filename * "_eta",
             cellfields=["eta_re" => real(ηₕ), "eta_im" => imag(ηₕ)],
             nsubcells=10)
  end

  xs_abs = collect(range(c.xb₀, c.xb₁, length=250))
  probes = [Point(x, 0.0, c.H) for x in xs_abs]
  ηvals = ηₕ.(probes)

  ξ = -2 .* (xs_abs .- (c.xb₀ + c.L / 2)) ./ c.L
  η_rel = abs.(ηvals) ./ wave.η₀

  if verbose
    @printf("\n=== Yago 3D case lambda/L = %.2f ===\n", λfactor)
    @printf("omega = %.6f rad/s, k = %.6f 1/m\n", wave.ω, wave.k)
    @printf("peak |eta|/eta0 = %.4f at xi = %.4f\n",
            maximum(η_rel), ξ[argmax(η_rel)])

    ξ_targets = [-1.0, -0.5, 0.0, 0.5, 1.0]
    for ξt in ξ_targets
      i = argmin(abs.(ξ .- ξt))
      @printf("  xi=%5.2f -> |eta|/eta0 = %.4f\n", ξ[i], η_rel[i])
    end
  end

  return ξ, η_rel
end

function run_yago_case1(; vtk_output=true)
  return run_yago_3d_freq(
    nx = 32,
    ny = 4,
    nz = 3,
    order = 4,
    λfactor = 0.4,
    dfactor = 4.0,
    vtk_output = vtk_output,
    verbose = true,
  )
end

function _plot_case1(xs, η_rel)
  outdir = joinpath(@__DIR__, "output")
  mkpath(outdir)
  plt = plot(
    xs,
    η_rel;
    xlabel = "xi = -2(x - x_c)/L",
    ylabel = "|eta|/eta0",
    label = "Yago 3D lambda/L=0.4",
    linewidth = 2,
  )
  hline!(plt, [1.0]; linestyle=:dash, color=:black, label="incident")
  savefig(plt, joinpath(outdir, "Yago3D_lambda04.png"))
end

if abspath(PROGRAM_FILE) == @__FILE__
  println("=== WARM-UP (mandatory) ===")
  @time run_yago_3d_freq(
    nx = 2,
    ny = 2,
    nz = 1,
    order = 2,
    λfactor = 0.4,
    dfactor = 2.0,
    vtk_output = false,
    verbose = false,
  )

  println("=== PRODUCTION CASE 1 (lambda/L = 0.4) ===")
  xs, η_rel = @time run_yago_case1(vtk_output=true)
  _plot_case1(xs, η_rel)
end

end # module
