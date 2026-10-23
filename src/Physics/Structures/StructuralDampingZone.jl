# =============================================================================
# StructuralDampingZone.jl
#
# Position-dependent damping on a hydroelastic structure, used as a numerical
# wave absorber ("sponge") at the ends of a structure that models an infinite
# one. A zone is a coefficient profile σ(x) over the structure's own measure;
# no separate triangulation is needed.
#
# Two mechanisms, both entering the `damping` form (−iω· in the frequency
# domain, ∂ₜ in the time domain):
#
#   dashpot      ∫ cᵨ σ v ηₜ dΓ          absorbs waves of every length
#   Kelvin-Voigt τ · k_op^σ(ηₜ, v)        stiffness-proportional; vanishes
#                                        for long waves (k → 0)
#
# where k_op^σ is the structure's elastic operator with its material
# coefficients multiplied by σ (see `stiffness_operator` with `weight`).
# =============================================================================

"""
    StructuralDampingZone

Position-dependent damping zone ("sponge") on an
[`AbstractHydroelasticStructure`](@ref), used to absorb waves at the ends of
a structure that represents an infinite one. It adds to the structure's
[`damping`](@ref) form

```math
\\int_{\\Gamma_\\eta} c_\\varrho\\,\\sigma(x)\\, v\\, \\partial_t\\eta
\\,\\mathrm{d}\\Gamma
\\;+\\; \\tau\\, k_{\\mathrm{op}}^{\\sigma}(\\partial_t\\eta, v),
```

where ``\\sigma(x) \\in [0, 1]`` is the zone profile (zero outside the zone)
and ``k_{\\mathrm{op}}^{\\sigma}`` is the structure's elastic operator with
its material coefficients multiplied by ``\\sigma``.

The viscous dashpot `cᵨ` damps waves of every length: to first order it
shifts the wavenumber by ``\\delta k \\approx i\\,\\omega c_\\varrho \\gamma_r``,
with ``\\gamma_r = 1/\\mathcal{D}_0'(k_0)`` the radiative width of the host.
The Kelvin-Voigt part `τ` is stiffness-proportional, so its effect vanishes
for long waves (``\\delta k \\propto k_0^4`` for bending).

Several zones (e.g. one at each end) are summed.

# Fields
- `profile::Function` — Zone profile `σ(x)` ∈ [0, 1], zero outside the zone;
  see [`ramp_profile`](@ref)
- `cᵨ::Float64` — Viscous dashpot coefficient / ρw [m/s]; default 0.0
- `τ::Float64` — Kelvin-Voigt (stiffness-proportional) coefficient [s];
  default 0.0
- `enabled::Bool` — Toggle for activating/deactivating this zone; default
  `true`

# Example
```julia
L, Ls = 200.0, 30.0
zones = [
  StructuralDampingZone(profile=ramp_profile(x_start=Ls, x_end=0.0), cᵨ=0.5),
  StructuralDampingZone(profile=ramp_profile(x_start=L - Ls, x_end=L), cᵨ=0.5),
]
beam = EulerBernoulliBeam(L=L, mᵨ=0.01, EIᵨ=0.49, damping_zones=zones)
```

See also: [`ramp_profile`](@ref), [`DampingZoneBC`](@ref) (the free-surface
counterpart)
"""
@with_kw struct StructuralDampingZone
  profile::Function
  cᵨ::Float64 = 0.0
  τ::Float64 = 0.0
  enabled::Bool = true
end

"""
    ramp_profile(; x_start, x_end, power=2, direction=nothing)

Ramp profile for a [`StructuralDampingZone`](@ref):

```math
\\sigma(x) = \\operatorname{clamp}\\!\\left(
\\frac{s(x) - s_{\\mathrm{start}}}{s_{\\mathrm{end}} - s_{\\mathrm{start}}},
0, 1\\right)^{p},
```

with ``s(x) = x_1`` or, when `direction` is given, the coordinate along that
direction. `σ = 0` on the inner side of `x_start` and `σ = 1` beyond `x_end`,
so `x_end < x_start` gives a zone at the lower end. A smooth ramp
(`power ≥ 2`, zone length of a few wavelengths) keeps the reflection from the
ramp itself small.

# Arguments
- `x_start::Real` — Coordinate where the ramp starts (`σ = 0`) [m]
- `x_end::Real` — Coordinate where the ramp reaches `σ = 1` [m]
- `power::Real` — Ramp exponent `p`; default 2
- `direction` — Ramp direction (any vector, normalised internally);
  default `nothing`, i.e. the first coordinate
"""
function ramp_profile(; x_start::Real, x_end::Real, power::Real=2,
                      direction=nothing)
  x_end == x_start && throw(ArgumentError("x_end must differ from x_start"))
  s0 = Float64(x_start)
  Δ = Float64(x_end) - s0
  p = Float64(power)
  if isnothing(direction)
    return x -> clamp((x[1] - s0) / Δ, 0.0, 1.0)^p
  end
  d = Tuple(Float64.(direction)) ./ sqrt(sum(abs2, direction))
  coordinate(x) = sum(d[i] * x[i] for i in eachindex(d))
  return x -> clamp((coordinate(x) - s0) / Δ, 0.0, 1.0)^p
end

_enabled_zones(zones) = filter(z -> z.enabled, zones)

# Sum of `coef(z) * z.profile(x)` over the zones with a nonzero coefficient,
# or `nothing` when there is none.
function _zone_coefficient(zones, coef::Function)
  active = [z for z in zones if !iszero(coef(z))]
  isempty(active) && return nothing
  return x -> sum(coef(z) * z.profile(x) for z in active)
end
