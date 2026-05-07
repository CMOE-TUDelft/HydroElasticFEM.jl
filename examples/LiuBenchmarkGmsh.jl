# ============================================================
# Liu VLFS Benchmark — Gmsh mesh version
# Reproduces Section 5.4 (test 5-3-1) of Colomés et al. (2022)
# DOI: 10.1002/nme.7140
#
# Reference implementation: oriolcg/MonolithicFEMVLFS.jl
#   src/Liu.jl, scripts/5-3-1-Liu.jl
#
# Physical setup:
#   2D hydroelastic problem — flexible beam floating on free surface
#   Domain: x ∈ [0, 13Lb], z ∈ [-H₀, 0]
#   Beam:   x ∈ [5Lb, 6Lb], z = 0
#   Incoming regular waves from the left
#   Sponge layers at inlet and outlet to absorb reflections
#
# Two frequency cases:
#   Case 1: ω = 0.4 rad/s
#   Case 2: ω = 0.8 rad/s
#
# Output: |η(x)|/η₀ along the beam, to be compared with
#   Liu et al. (1992) experimental/numerical data.
#
# Damping sponge layers are implemented via DampingZoneBC on the
# PotentialFlow entity.
# ============================================================

module LiuBenchmarkGmsh

using Gridap
using Parameters
using Printf
using Plots

import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
using HydroElasticFEM: PKG_ROOT
import GridapGmsh: gmsh as _gmsh_api

export LiuCaseParams, run_liu_case, run_liu_two_cases, generate_liu_mesh

# ──────────────────────────────────────────────────────────────────────────
# Gmsh mesh generator for the Liu domain
# ──────────────────────────────────────────────────────────────────────────

# Library path captured at module load time (required to be const for ccall).
# GridapGmsh v4.9.3 wrapper lacks the name parameter added in libgmsh >= 4.13,
# causing a segfault. Call the C function directly with the v4.15 signature.
const _GMSH_LIB = _gmsh_api.lib

function _add_physical_group(dim::Int32, tags::Vector{Int32}, name::String)
  tag  = Int32(-1)
  ierr = Ref{Cint}()
  result = ccall(
    (:gmshModelAddPhysicalGroup, _GMSH_LIB), Cint,
    (Cint, Ptr{Cint}, Csize_t, Cint, Cstring, Ptr{Cint}),
    dim, tags, length(tags), tag, name, ierr,
  )
  ierr[] != 0 && error("gmsh addPhysicalGroup failed")
  return result
end

"""
    generate_liu_mesh(; Lb, H0, filename, lc_beam, lc_fluid)

Generate a 2-D Gmsh mesh for the Liu VLFS benchmark domain.

Domain: x ∈ [0, 13Lb], z ∈ [-H0, 0].  Physical groups assigned:
`fluid`, `seabed`, `inlet`, `outlet`, `damping_in`, `damping_out`,
`free_surface`, `structure`.
"""
function generate_liu_mesh(;
  Lb::Float64       = 300.0,
  H0::Float64       = 60.0,
  filename::String  = "liu_mesh.msh",
  lc_beam::Float64  = -1.0,
  lc_fluid::Float64 = -1.0,
)
  Ld   = 4.0 * Lb
  lc_b = lc_beam  < 0.0 ? Lb / 50.0 : lc_beam
  lc_f = lc_fluid < 0.0 ? H0 / 5.0  : lc_fluid
  _generate_liu_mesh_impl(Lb, H0, Ld, filename, lc_b, lc_f)
end

function _generate_liu_mesh_impl(
  Lb::Float64, H0::Float64, Ld::Float64,
  filename::String, lc_beam::Float64, lc_fluid::Float64,
)
  gmsh = _gmsh_api
  gmsh.initialize()
  gmsh.option.setNumber("General.Terminal", 0)
  try
    x0    = 0.0
    x_d1  = Ld
    x_b0  = Ld + Lb
    x_b1  = Ld + 2 * Lb
    x_d2  = Ld + 5 * Lb
    x_end = Ld + 9 * Lb
    lc_s  = lc_beam
    lc_m  = lc_beam * 3
    lc_c  = lc_fluid
    p1 = gmsh.model.geo.addPoint(x0,    -H0, 0.0, lc_c)
    p2 = gmsh.model.geo.addPoint(x_end, -H0, 0.0, lc_c)
    p3 = gmsh.model.geo.addPoint(x0,    0.0, 0.0, lc_m)
    p4 = gmsh.model.geo.addPoint(x_d1,  0.0, 0.0, lc_m)
    p5 = gmsh.model.geo.addPoint(x_b0,  0.0, 0.0, lc_s)
    p6 = gmsh.model.geo.addPoint(x_b1,  0.0, 0.0, lc_s)
    p7 = gmsh.model.geo.addPoint(x_d2,  0.0, 0.0, lc_m)
    p8 = gmsh.model.geo.addPoint(x_end, 0.0, 0.0, lc_c)
    l_seabed = gmsh.model.geo.addLine(p1, p2)
    l_outlet = gmsh.model.geo.addLine(p2, p8)
    l_dout   = gmsh.model.geo.addLine(p8, p7)
    l_fs2    = gmsh.model.geo.addLine(p7, p6)
    l_struct = gmsh.model.geo.addLine(p6, p5)
    l_fs1    = gmsh.model.geo.addLine(p5, p4)
    l_din    = gmsh.model.geo.addLine(p4, p3)
    l_inlet  = gmsh.model.geo.addLine(p3, p1)
    cl = gmsh.model.geo.addCurveLoop([
      l_seabed, l_outlet, l_dout, l_fs2, l_struct, l_fs1, l_din, l_inlet,
    ])
    surf = gmsh.model.geo.addPlaneSurface([cl])
    gmsh.model.geo.synchronize()
    _add_physical_group(Int32(2), Int32[surf],          "fluid")
    _add_physical_group(Int32(1), Int32[l_seabed],      "seabed")
    _add_physical_group(Int32(1), Int32[l_inlet],       "inlet")
    _add_physical_group(Int32(1), Int32[l_outlet],      "outlet")
    _add_physical_group(Int32(1), Int32[l_din],         "damping_in")
    _add_physical_group(Int32(1), Int32[l_dout],        "damping_out")
    _add_physical_group(Int32(1), Int32[l_fs1, l_fs2],  "free_surface")
    _add_physical_group(Int32(1), Int32[l_struct],      "structure")
    n_beam = max(2, Int(round(abs(x_b1 - x_b0) / lc_beam)) + 1)
    gmsh.model.mesh.setTransfiniteCurve(l_struct, n_beam)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)
    gmsh.write(filename)
  finally
    gmsh.finalize()
  end
  return filename
end

# ──────────────────────────────────────────────────────────────────────────
# Parameters struct
# ──────────────────────────────────────────────────────────────────────────

@with_kw struct LiuCaseParams
  name::String                      = "LiuBenchmark"
  ω::Float64                        = 0.4       # angular frequency [rad/s]
  lc_beam::Float64                  = -1.0      # negative → use Lb/50
  vtk_output::Bool                  = false
  make_plot::Bool                   = true
  output_dir::Union{String,Nothing} = nothing   # nothing → default PKG_ROOT path
end

# ──────────────────────────────────────────────────────────────────────────
# Physical constants (Section 5.4 of Colomés et al. 2022)
# ──────────────────────────────────────────────────────────────────────────

function _constants()
  Lb  = 300.0   # [m]      beam length
  H₀  = 60.0    # [m]      water depth
  m   = 500.0   # [kg/m²]  beam mass per unit area
  EI  = 1.0e10  # [N·m]    beam bending stiffness
  g   = 9.81    # [m/s²]
  ρ   = 1025.0  # [kg/m³]  water density
  d₀  = m / ρ   # [m]      equivalent draft
  a₁  = EI / ρ  # [m⁵/s²]  structural parameter (EI/ρ)
  Ld  = 4 * Lb  # [m]      damping zone length
  η₀  = 0.01    # [m]      incident wave amplitude
  (; Lb, H₀, m, EI, g, ρ, d₀, a₁, Ld, η₀)
end

# ──────────────────────────────────────────────────────────────────────────
# Wave parameters (depend on ω)
# ──────────────────────────────────────────────────────────────────────────

function _wave_params(c, ω)
  # Monotone dispersion relation in k > 0: g*k*tanh(k*H0) - ω^2 = 0.
  f(k) = c.g * k * tanh(k * c.H₀) - ω^2

  klo = 1.0e-8
  khi = 5.0
  while f(khi) < 0.0
    khi *= 2.0
    khi > 1.0e4 && error("Failed to bracket dispersion root")
  end

  for _ in 1:80
    kmid = 0.5 * (klo + khi)
    if f(kmid) > 0.0
      khi = kmid
    else
      klo = kmid
    end
  end

  k = 0.5 * (klo + khi)
  λ = 2π / k

  η₀ = c.η₀

  ηin(x)  = η₀ * exp(im * k * x[1])
  ϕin(x)  = -im * (η₀ * ω / k) *
              (cosh(k * (x[2] + 0.075 * c.Lb)) / sinh(k * c.H₀)) *
              exp(im * k * x[1])
  vin(x)  = (η₀ * ω) *
              (cosh(k * (x[2] + 0.075 * c.Lb)) / sinh(k * c.H₀)) *
              exp(im * k * x[1])
  vzin(x) = -im * ω * η₀ * exp(im * k * x[1])

  (; k, λ, ηin, ϕin, vin, vzin)
end

# ──────────────────────────────────────────────────────────────────────────
# Numerical parameters
# ──────────────────────────────────────────────────────────────────────────

function _numerics(c, ω, lc_beam_arg)
  order  = 4
  h      = lc_beam_arg < 0.0 ? c.Lb / 50.0 : lc_beam_arg
  γ      = 1.0 * order * (order - 1) / h  # C/DG penalty parameter
  βₕ     = 0.5
  αₕ     = -im * ω / c.g * (1 - βₕ) / βₕ
  (; order, h, γ, βₕ, αₕ)
end

# ──────────────────────────────────────────────────────────────────────────
# Sponge layer (damping) functions
# ──────────────────────────────────────────────────────────────────────────

function _damping(c, wave)
  μ₀    = 6.0
  Ld    = c.Ld
  xdout = 9 * c.Lb   # x-coordinate where outlet damping begins

  μ₁in(x)  = μ₀ * (1.0 - sin(π / 2 * x[1] / Ld))
  μ₁out(x) = μ₀ * (1.0 - cos(π / 2 * (x[1] - xdout) / Ld))
  μ₂in(x)  = μ₁in(x) * wave.k
  μ₂out(x) = μ₁out(x) * wave.k

  (; μ₁in, μ₁out, μ₂in, μ₂out)
end

# ──────────────────────────────────────────────────────────────────────────
# Main solver
# ──────────────────────────────────────────────────────────────────────────

"""
    run_liu_case(params::LiuCaseParams)

Run one frequency-domain hydroelastic case for the Liu VLFS benchmark.

Returns `(xs, eta_rel, meta)` where:
- `xs`      are normalised beam coordinates `(x - 5Lb) / Lb` ∈ [0, 1]
- `eta_rel` is `|η(x)| / η₀` sampled along the beam
- `meta`    is a named tuple with derived physical constants

Reference: Section 5.4 of Colomés, Verdugo, Akkerman (2022), DOI:10.1002/nme.7140
"""
function run_liu_case(params::LiuCaseParams)
  c    = _constants()
  wave = _wave_params(c, params.ω)
  num  = _numerics(c, params.ω, params.lc_beam)
  damp = _damping(c, wave)

  # ── 1. Generate mesh ─────────────────────────────────────────────────────
  msh_dir = if params.output_dir !== nothing
    joinpath(params.output_dir, params.name)
  else
    joinpath(PKG_ROOT, "data", "VTK", "examples", "LiuBenchmarkGmsh", params.name)
  end
  isdir(msh_dir) || mkpath(msh_dir)
  msh_file = joinpath(msh_dir, "liu_mesh.msh")

  generate_liu_mesh(
    Lb       = c.Lb,
    H0       = c.H₀,
    filename = msh_file,
    lc_beam  = num.h,
    lc_fluid = c.H₀ / 5.0,
  )

  # ── 2. Load mesh ──────────────────────────────────────────────────────────
  domain = G.GmshDomain(msh_file; dim = 2)

  # ── 3. Physics entities ───────────────────────────────────────────────────
  # Incident-wave velocity at the inlet boundary (∫(w * vin)dΓin)
  vin_fn  = wave.vin

  # Zero incident conditions for the outlet damping zone
  zero_fn = x -> 0.0 + 0.0im

  potential = P.PotentialFlow(
    ρw = c.ρ,
    g  = c.g,
    boundary_conditions = P.AbstractPotentialFlowBC[
      # Inlet: prescribed normal velocity (Neumann)
      P.PrescribedInletPotentialBC(
        domain   = :dΓin,
        forcing  = vin_fn,
        quantity = :normal_gradient,
      ),
      # Inlet sponge layer (damping_in → :dΓd_1)
      P.DampingZoneBC(
        domain  = :dΓd_1,
        μ₁      = damp.μ₁in,
        μ₂      = damp.μ₂in,
        η_in    = wave.ηin,
        vz_in   = wave.vzin,
      ),
      # Outlet sponge layer (damping_out → :dΓd_2) — no incident wave
      P.DampingZoneBC(
        domain  = :dΓd_2,
        μ₁      = damp.μ₁out,
        μ₂      = damp.μ₂out,
        η_in    = zero_fn,
        vz_in   = zero_fn,
      ),
    ],
    fe             = PH.FESpaceConfig(
      order        = num.order,
      vector_type  = Vector{ComplexF64},
    ),
    space_domain_symbol = :Ω,
  )

  free_surface = P.FreeSurface(
    ρw   = c.ρ,
    g    = c.g,
    βₕ   = num.βₕ,
    fe   = PH.FESpaceConfig(
      order       = num.order,
      vector_type = Vector{ComplexF64},
    ),
    space_domain_symbol = :Γκ,
  )

  beam = P.EulerBernoulliBeam(
    L    = c.Lb,
    mᵨ   = c.d₀,
    EIᵨ  = c.a₁,
    g    = c.g,
    fe   = PH.FESpaceConfig(
      order       = num.order,
      vector_type = Vector{ComplexF64},
      γ           = num.γ,
    ),
    space_domain_symbol = :Γη,
  )

  physics = P.PhysicsParameters[potential, free_surface, beam]
  config  = S.FreqDomainConfig(ω = params.ω)

  # ── 4. Build and solve ────────────────────────────────────────────────────
  problem  = S.build_problem(domain, physics, config)
  result   = S.simulate(problem)

  ϕh, κh, ηh = result.solution

  # ── 5. Optional VTK output ────────────────────────────────────────────────
  if params.vtk_output
    trians = S.get_triangulations(problem)
    writevtk(trians[:Ω],  joinpath(msh_dir, "phi"),
             cellfields = ["phi_re" => real(ϕh), "phi_im" => imag(ϕh)])
    writevtk(trians[:Γη], joinpath(msh_dir, "eta"),
             cellfields = ["eta_re" => real(ηh), "eta_im" => imag(ηh)])
    writevtk(trians[:Γκ], joinpath(msh_dir, "kappa"),
             cellfields = ["kappa_re" => real(κh), "kappa_im" => imag(κh)])
  end

  # ── 6. Post-processing: |η| / η₀ along beam ──────────────────────────────
  x_beam_start = 5 * c.Lb
  x_beam_end   = 6 * c.Lb
  n_probes     = 400
  xs_abs  = collect(range(x_beam_start, x_beam_end, length = n_probes))
  probes  = [Point(x, 0.0) for x in xs_abs]

  ηvals   = ηh.(probes)
  xs_norm = (xs_abs .- x_beam_start) ./ c.Lb  # ∈ [0, 1]
  eta_rel = abs.(ηvals) ./ c.η₀

  meta = (;
    ω    = params.ω,
    k    = wave.k,
    λ    = wave.λ,
    η₀   = c.η₀,
    Lb   = c.Lb,
    h    = num.h,
    γ    = num.γ,
    x_beam_start,
    x_beam_end,
  )

  return xs_norm, eta_rel, meta
end

# ──────────────────────────────────────────────────────────────────────────
# Run two canonical cases and produce figures
# ──────────────────────────────────────────────────────────────────────────

"""
    run_liu_two_cases(; kwargs...)

Run both canonical Liu benchmark cases (ω = 0.4 rad/s and ω = 0.8 rad/s)
and generate comparison plots.

Plots are saved to `examples/output/Liu_omega04.png` and
`examples/output/Liu_omega08.png`.

# Keyword arguments forwarded to `LiuCaseParams`:
- `lc_beam`    : beam mesh size [m] (default Lb/50 ≈ 6 m, coarse for speed)
- `vtk_output` : write VTK files (default `false`)
"""
function run_liu_two_cases(;
  lc_beam    = -1.0,
  vtk_output = false,
)
  outdir = joinpath(PKG_ROOT, "examples", "output")
  isdir(outdir) || mkpath(outdir)

  cases = [
    (ω = 0.4, name = "omega04", label = "ω = 0.4 rad/s", file = "Liu_omega04.png"),
    (ω = 0.8, name = "omega08", label = "ω = 0.8 rad/s", file = "Liu_omega08.png"),
  ]

  results = map(cases) do c
    @printf("Running Liu benchmark case: %s\n", c.label)
    params = LiuCaseParams(
      name       = "Liu_" * c.name,
      ω          = c.ω,
      lc_beam    = lc_beam,
      vtk_output = vtk_output,
      make_plot  = true,
    )
    xs, eta_rel, meta = run_liu_case(params)
    @printf("  λ/Lb = %.3f,  peak |η|/η₀ = %.3f,  min |η|/η₀ = %.3f\n",
            meta.λ / meta.Lb, maximum(eta_rel), minimum(eta_rel))
    (xs = xs, eta_rel = eta_rel, meta = meta, c...)
  end

  # ── Individual plots ───────────────────────────────────────────────────────
  for r in results
    plt = plot(
      r.xs, r.eta_rel;
      xlabel       = "x / Lb",
      ylabel       = "|η| / η₀",
      title        = "Liu VLFS Benchmark — $(r.label)",
      label        = r.label,
      lw           = 2,
      legend       = :topright,
      framestyle   = :box,
    )
    hline!([1.0]; ls = :dash, lc = :gray, label = "η = η₀")
    savefig(plt, joinpath(outdir, r.file))
    @printf("  Saved: %s\n", joinpath(outdir, r.file))
  end

  # ── Summary table ──────────────────────────────────────────────────────────
  println()
  println("─────────────────────────────────────────────────")
  println("  Liu VLFS Benchmark — Summary")
  println("─────────────────────────────────────────────────")
  @printf("  %-16s  %-10s  %-10s\n", "Case", "peak |η|/η₀", "min |η|/η₀")
  println("  " * "─"^43)
  for r in results
    @printf("  %-16s  %-10.3f  %-10.3f\n",
            r.label, maximum(r.eta_rel), minimum(r.eta_rel))
  end
  println("─────────────────────────────────────────────────")

  return results
end

end # module LiuBenchmarkGmsh

# ── Entry point when run as a script ──────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
  LiuBenchmarkGmsh.run_liu_two_cases()
end
