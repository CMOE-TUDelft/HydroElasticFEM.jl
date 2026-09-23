# ─────────────────────────────────────────────────────────────
# Physical reflection / transmission / absorption coefficients, built on
# top of the generic fitting utilities in WaveDecomposition.jl.
# ─────────────────────────────────────────────────────────────

"""
    ReflectionTransmissionResult

Bundled result of a single-frequency reflection/transmission/absorption
computation, returned by
[`reflection_transmission_coefficients`](@ref).

# Fields
- `R::ComplexF64` — complex reflection coefficient (reflected / fitted
  incident amplitude).
- `T::ComplexF64` — complex transmission coefficient (transmitted / fitted
  incident amplitude).
- `A::Float64` — energy-balance absorption coefficient,
  `1 - abs2(R) - abs2(T)` (see [`absorption_coefficient`](@ref)).
- `k::Float64` — wavenumber used.
- `incident_fit::WaveProbeFit` — fit on the up-wave probe line. Its `.A` is
  the fitted incident amplitude (compare against a nominal `η0` as a
  sanity check); its `.B` is the reflected amplitude used for `R`.
- `transmitted_fit::WaveProbeFit` — fit on the down-wave probe line. Its
  `.A` is the transmitted amplitude used for `T`; its `.B` should be small
  and is a diagnostic for spurious reflections re-entering from the far
  boundary.
"""
struct ReflectionTransmissionResult
    R::ComplexF64
    T::ComplexF64
    A::Float64
    k::Float64
    incident_fit::WaveProbeFit
    transmitted_fit::WaveProbeFit
end

function Base.show(io::IO, r::ReflectionTransmissionResult)
    @printf(io, "ReflectionTransmissionResult(|R|=%.4f, |T|=%.4f, |R|²+|T|²=%.4f, A=%.4f)",
        abs(r.R), abs(r.T), abs2(r.R) + abs2(r.T), r.A)
end

"""
    absorption_coefficient(R, T) -> Float64

Energy-balance absorption coefficient `A = 1 - |R|² - |T|²`.

Valid whenever the fluid domain and its (possibly damped) structure
conserve energy except for what leaves as reflected/transmitted propagating
waves — i.e. this measures *all* energy not accounted for by `R` and `T`,
whatever its source (structural damping, a `ResonatorSingle` PTO, or
numerical/discretisation loss masquerading as physical dissipation). For a
purely elastic, undamped structure this should come out `≈ 0`, which
doubles as a solver/post-processing verification check — see the energy
check printed by [`reflection_transmission_coefficients`](@ref).
"""
absorption_coefficient(R::Number, T::Number) = 1 - abs2(R) - abs2(T)

"""
    reflection_transmission_coefficients(xs_in, vals_in, xs_out, vals_out, k;
                                          η0=nothing, cond_tol=50.0) -> ReflectionTransmissionResult
    reflection_transmission_coefficients(field, xs_in, xs_out, k;
                                          y=0.0, η0=nothing, cond_tol=50.0) -> ReflectionTransmissionResult

Compute the reflection coefficient `R`, transmission coefficient `T`, and an
energy-balance absorption coefficient `A` from a known, single-frequency
free-surface elevation field, given probe points on the up-wave ("in") and
down-wave ("out") sides of a structure.

Two methods are provided:

1. **Raw-vector method** — pass already-sampled complex values directly
   (`vals_in`, `vals_out`). Use this if you obtained the free-surface
   samples some other way (e.g. a time-domain FFT amplitude at the forcing
   frequency), or already sampled a Gridap field yourself.
2. **Field method** — pass a Gridap `CellField`/`FEFunction` (or any
   callable accepting a `Gridap.Point`) plus x-coordinates; this internally
   calls [`sample_probe_line`](@ref) at vertical coordinate `y`.

# Arguments
- `xs_in`, `xs_out::AbstractVector{<:Real}`: probe x-coordinates on the
  up-wave and down-wave sides, `N ≥ 2` each, in the constant-depth region
  beyond evanescent-mode decay (see [`probe_positions`](@ref)).
- `vals_in`, `vals_out::AbstractVector{<:Number}`: sampled complex
  free-surface elevation at `xs_in`/`xs_out` (raw-vector method only).
- `field`: a Gridap `CellField`/`FEFunction` giving the complex free-surface
  elevation — typically `κh` (open water) or `ηh` (over a structure) from
  `result.solution` — or any callable `Gridap.Point -> Complex`
  (field method only).
- `k::Real`: wavenumber of the incident wave. Use the same value the
  simulation itself used for the radiation boundary condition (e.g.
  `sea_state.k[1]`), not an independently recomputed one, so post-processing
  can never drift out of sync with the solve.

# Keyword arguments
- `y::Real=0.0`: vertical coordinate to sample `field` at (field method
  only). **This must match the free-surface elevation for your mesh/domain
  convention.** In particular, `TankDomain`'s default (unshifted) Cartesian
  mesh has the free surface at `z = H0`, seabed at `z = 0` — *not* the more
  common SWL-at-zero convention some hand-written incident-wave formulas
  assume. Check where your `Γκ`/`Γη` triangulation actually sits before
  trusting the default.
- `η0::Union{Real,Nothing}=nothing`: nominal incident amplitude. If given, a
  warning is issued when it differs from the *fitted* incident amplitude
  `abs(incident_fit.A)` by more than 5% — a useful mesh/probe-placement
  sanity check. Does not affect `R`/`T`, which are always normalised by the
  fitted amplitude for self-consistency (so `R`/`T` remain correct even if
  you don't know `η0` precisely).
- `cond_tol::Real=50.0`: throws an informative error if either probe line's
  fit condition number exceeds this. Almost always caused by a probe
  spacing landing near a multiple of half the local wavelength — see
  [`suggest_probe_offsets`](@ref).

# Returns
A [`ReflectionTransmissionResult`](@ref).

# Method
Each side is decomposed independently with [`fit_wave_components`](@ref):
the up-wave line as `A_in*exp(ikx) + B_in*exp(-ikx)` (incident + reflected),
the down-wave line as `A_out*exp(ikx) + B_out*exp(-ikx)` (transmitted +
anything spuriously reflected back from downstream). Then

    R = B_in  / A_in
    T = A_out / A_in
    A = 1 - |R|² - |T|²

`|B_out|` is not used in the coefficients but is available via
`result.transmitted_fit.B` as a diagnostic: it should be small relative to
`|A_out|`, and a large value signals spurious reflection re-entering the
transmission-side probe line (e.g. from an imperfect outlet condition or a
downstream obstruction).

# Example
```julia
ϕh, κh, ηh = result.solution
Lwave  = 2π / k
xs_in  = probe_positions(xb0, :upwave,   Lwave, H0)
xs_out = probe_positions(xb1, :downwave, Lwave, H0)

rt = reflection_transmission_coefficients(κh, xs_in, xs_out, k; y=H0, η0=0.1)
println(rt)                       # |R|=..., |T|=..., |R|²+|T|²=..., A=...
abs2(rt.R) + abs2(rt.T)            # → 1 for an undamped structure
```

# References
- Goda, Y., & Suzuki, Y. (1976). Estimation of incident and reflected waves
  in random wave experiments. *Proc. 15th Int. Conf. Coastal Eng.*, 828-845.
- Zelt, J. A., & Skjelbreia, J. E. (1992). Estimating incident and reflected
  wave fields using an arbitrary number of wave gauges.
  *Proc. 23rd Int. Conf. Coastal Eng.*, 777-789.
"""
function reflection_transmission_coefficients(
    xs_in::AbstractVector{<:Real}, vals_in::AbstractVector{<:Number},
    xs_out::AbstractVector{<:Real}, vals_out::AbstractVector{<:Number},
    k::Real;
    η0::Union{Real,Nothing}=nothing,
    cond_tol::Real=50.0,
)
    fit_in = fit_wave_components(xs_in, vals_in, k)
    fit_out = fit_wave_components(xs_out, vals_out, k)

    fit_in.cond > cond_tol && error(
        "Up-wave probe fit is ill-conditioned (cond=$(round(fit_in.cond, digits=1)) > " *
        "$cond_tol). Probe spacings likely land near a multiple of half the local " *
        "wavelength (λ = $(round(2π / k, digits=3))); see `suggest_probe_offsets`.",
    )
    fit_out.cond > cond_tol && error(
        "Down-wave probe fit is ill-conditioned (cond=$(round(fit_out.cond, digits=1)) > " *
        "$cond_tol). See the up-wave error message for the likely cause.",
    )

    R = fit_in.B / fit_in.A
    T = fit_out.A / fit_in.A
    Acoef = absorption_coefficient(R, T)

    if η0 !== nothing
        rel_err = abs(abs(fit_in.A) - η0) / η0
        if rel_err > 0.05
            @warn "Fitted incident amplitude differs from the nominal η0 by " *
                  "$(round(100 * rel_err, digits=1))%. Check probe placement, mesh " *
                  "resolution, or that `k` matches the simulation's sea_state." fitted =
                abs(fit_in.A) nominal = η0
        end
    end

    return ReflectionTransmissionResult(R, T, Acoef, Float64(k), fit_in, fit_out)
end

function reflection_transmission_coefficients(
    field, xs_in::AbstractVector{<:Real}, xs_out::AbstractVector{<:Real}, k::Real;
    y::Real=0.0, kwargs...,
)
    vals_in = sample_probe_line(field, xs_in, y)
    vals_out = sample_probe_line(field, xs_out, y)
    return reflection_transmission_coefficients(xs_in, vals_in, xs_out, vals_out, k; kwargs...)
end

# ─────────────────────────────────────────────────────────────
# Energy flux / dissipation-based absorption
#
# An alternative (and cross-check) to the energy-balance `absorption_coefficient`
# above: compute A directly from the power dissipated by a damped structural
# element (a `ResonatorSingle` PTO, or beam/membrane structural damping),
# divided by the incident wave power. The two routes should roughly agree
# for a converged solve; a persistent mismatch usually means the far-field
# probes are contaminated by evanescent modes or too close to the structure.
# ─────────────────────────────────────────────────────────────

"""
    group_velocity(k, H; g=9.81) -> Float64

Linear finite-depth group velocity,

    cg = (ω / 2k) * (1 + 2kH / sinh(2kH)),   ω = sqrt(g k tanh(kH))

reducing to `cg → ω/(2k)` (half the phase speed) in deep water (`kH ≫ 1`)
and `cg → sqrt(gH)` in shallow water (`kH ≪ 1`).
"""
function group_velocity(k::Real, H::Real; g::Real=9.81)
    ω = sqrt(g * k * tanh(k * H))
    return (ω / (2k)) * (1 + 2k * H / sinh(2k * H))
end

"""
    energy_flux(η0, k, H; ρw=1025.0, g=9.81) -> Float64

Time-averaged incident linear-wave energy flux per unit crest length,

    P = 1/2 * ρw * g * η0^2 * cg

with `cg` from [`group_velocity`](@ref).
"""
function energy_flux(η0::Real, k::Real, H::Real; ρw::Real=1025.0, g::Real=9.81)
    return 0.5 * ρw * g * η0^2 * group_velocity(k, H; g=g)
end

"""
    absorption_from_dissipated_power(P_diss, η0, k, H; ρw=1025.0, g=9.81) -> Float64

`A = P_diss / P_incident`, an alternative route to the energy-balance
estimate `1 - |R|² - |T|²` in [`absorption_coefficient`](@ref), for use when
a specific dissipative element's power is known (e.g. from
[`resonator_dissipated_power`](@ref)). Cross-checking the two is a good
verification of both the FE solve and the far-field probe placement.
"""
function absorption_from_dissipated_power(P_diss::Real, η0::Real, k::Real, H::Real;
                                           ρw::Real=1025.0, g::Real=9.81)
    return P_diss / energy_flux(η0, k, H; ρw=ρw, g=g)
end

"""
    resonator_dissipated_power(qh, ω, C) -> Float64

Time-averaged power dissipated by a linearly-damped resonator or
power-take-off element with damping coefficient `C` (the `C` field of
`ResonatorSingle`) and complex frequency-domain response amplitude
`qh` (from `result.solution`), under the `x(t) = Re(x̂ * exp(iωt))`
convention used throughout HydroElasticFEM.jl:

    P = 1/2 * C * ω^2 * |qh|^2
"""
resonator_dissipated_power(qh::Number, ω::Real, C::Real) = 0.5 * C * ω^2 * abs2(qh)
