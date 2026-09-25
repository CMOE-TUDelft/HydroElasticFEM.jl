module FloatingTensionedBeamExample

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
Example: floating tensioned Euler-Bernoulli beam
==================================================

Compares three structural models — `Membrane`, `EulerBernoulliBeam`, and the
new `TensionedEulerBernoulliBeam` — under identical wave conditions, on the
same tank geometry and structure span. `TensionedEulerBernoulliBeam`
interpolates between the two: setting `EIᵨ = 0` recovers the membrane
stiffness operator exactly, and setting `Tᵨ = 0` recovers the plain
Euler-Bernoulli beam stiffness operator exactly (verified at the matrix
level in `test/examples/TensionedEulerBernoulliBeamWeakFormTests.jl`; see
`docs/src/guide/theory.md` for the governing equation and asymptotic
regimes).

Follows the structured `HydroElasticFEM` API demonstrated in
`FloatingMembraneExample.run_structured_implementation`.

# Usage
```julia
include("examples/FloatingTensionedBeamExample.jl")
using .FloatingTensionedBeamExample
results = compare_structures()
```
"""

@with_kw struct FloatingTensionedBeamParams
  H0::Float64 = 10.0
  LΩ::Float64 = 60.0
  x0::Float64 = 0.0
  ω::Float64 = 1.2
  η0::Float64 = 0.25
  α::Float64 = 0.0
  βₕ::Float64 = 0.5
  # Order 2 is required for the C/DG Euler-Bernoulli bending operator
  # (order 1 gives Δη ≡ 0 everywhere); used uniformly across all three
  # structural models for a like-for-like comparison.
  order::Int = 2
  nx::Int = 240
  ny::Int = 12
  mesh_ry::Float64 = 1.08
  Lb::Float64 = 20.0
  xb0::Float64 = 20.0
  mᵨ::Float64 = 922.5 / 1025.0
  EIᵨ::Float64 = 1.0e5 / 1025.0
  Tᵨ::Float64 = 10.0 * 9.81
  τ::Float64 = 0.0
  probe_x::Vector{Float64} = collect(range(5.0, 55.0, length=9))
end

"""
    build_regular_wave_state(; H, T, h, θ=0.0)

Build a monochromatic regular-wave state using `WaveSpec` (identical helper
to `FloatingMembraneExample.build_regular_wave_state`; kept local so this
example script remains self-contained).
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
    shifted_gp_map(x0, H0, mesh_ry, ny)

Vertical mesh-grading map for the structured tank, including the horizontal
shift `x0` and the `[0, H0]` vertical convention used by `TankDomain`.
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
probe_points(xs) = Point.(xs, 0.0)

"""
    beam_indicator(xs, xb0, xb1)

Return a boolean mask identifying which probe coordinates lie on the beam
span `[xb0, xb1]`.
"""
beam_indicator(xs, xb0, xb1) = (xb0 .<= xs) .& (xs .<= xb1)

"""
    build_structure(kind::Symbol, p::FloatingTensionedBeamParams)

Construct the requested structural model (`:membrane`, `:euler_beam`, or
`:tensioned_beam`) with parameters drawn from `p`, sharing `L`, `τ`, `g`, and
`fe` across all three so the comparison in [`compare_structures`](@ref) is
like-for-like.
"""
function build_structure(kind::Symbol, p::FloatingTensionedBeamParams)
  fe = PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64})
  g  = WaveSpec.PhysicalConstants.g

  if kind == :membrane
    return P.Membrane(L=p.Lb, mᵨ=p.mᵨ, Tᵨ=p.Tᵨ, τ=p.τ, g=g, fe=fe)
  elseif kind == :euler_beam
    return P.EulerBernoulliBeam(L=p.Lb, mᵨ=p.mᵨ, EIᵨ=p.EIᵨ, τ=p.τ, g=g, fe=fe)
  elseif kind == :tensioned_beam
    return P.TensionedEulerBernoulliBeam(L=p.Lb, mᵨ=p.mᵨ, EIᵨ=p.EIᵨ, Tᵨ=p.Tᵨ, τ=p.τ, g=g, fe=fe)
  else
    error("Unknown structure kind: $kind. Use :membrane, :euler_beam, or :tensioned_beam.")
  end
end

"""
    run_case(kind::Symbol; kwargs...)

Solve the floating-structure frequency-domain problem for one structural
model (`:membrane`, `:euler_beam`, or `:tensioned_beam`) using the structured
`HydroElasticFEM` API. Keyword arguments override
[`FloatingTensionedBeamParams`](@ref).

Returns a named tuple with the combined free-surface/structure probe values
and the fluid velocity-potential probe values.
"""
function run_case(kind::Symbol; kwargs...)
  p = FloatingTensionedBeamParams(; kwargs...)

  probes = probe_points(p.probe_x)
  xb1 = p.xb0 + p.Lb
  probe_on_beam = beam_indicator(p.probe_x, p.xb0, xb1)

  tank = G.TankDomain(
    L=p.LΩ,
    H=p.H0,
    nx=p.nx,
    ny=p.ny,
    map=shifted_gp_map(p.x0, p.H0, p.mesh_ry, p.ny),
    structure_domains=[
      G.StructureDomain(L=p.Lb, x₀=[p.xb0, 0.0], domain_symbol=:Γb),
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

  structure = build_structure(kind, p)

  physics = P.PhysicsParameters[potential, free_surface, structure]
  config = S.FreqDomainConfig(ω=p.ω)
  problem = S.build_problem(tank, physics, config)

  result = S.simulate(problem)
  ϕₕ, κₕ, ηₕ = result.solution

  probe_surface = similar(ϕₕ(probes))
  probe_surface[.!probe_on_beam] = κₕ(probes[.!probe_on_beam])
  probe_surface[probe_on_beam]   = ηₕ(probes[probe_on_beam])

  return (; kind, params=p, probe_surface, probe_potential=ϕₕ(probes))
end

"""
    compare_structures(; kwargs...)

Run [`run_case`](@ref) for `:membrane`, `:euler_beam`, and `:tensioned_beam`
under identical wave conditions and geometry, print the maximum combined
surface/structure probe response for each, and return a `Dict{Symbol,Any}`
of the full results keyed by structure kind.

Keyword arguments override [`FloatingTensionedBeamParams`](@ref) and are
forwarded to every case.
"""
function compare_structures(; kwargs...)
  results = Dict{Symbol,Any}()
  @printf("%-16s %s\n", "structure", "max |combined probe response|")
  for kind in (:membrane, :euler_beam, :tensioned_beam)
    results[kind] = run_case(kind; kwargs...)
    max_resp = maximum(abs.(results[kind].probe_surface))
    @printf("%-16s %.6f\n", String(kind), max_resp)
  end
  return results
end

end # module
