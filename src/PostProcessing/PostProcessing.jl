"""
    module PostProcessing

Frequency-domain post-processing utilities for wave-structure interaction
results: reflection, transmission, and absorption coefficients.

This module is deliberately independent of `Physics`/`Simulation` — it never
builds or solves anything. It operates purely on values you already have: a
known wavenumber `k` and a free-surface elevation field (or a set of
already-sampled complex values) at a handful of probe points. That means it
works the same whether the field came from `result.solution` (`κh`/`ηh`),
from re-loaded VTK output, from a different solver entirely, or from a
time-domain simulation post-processed with your own FFT.

## Two layers

**`WaveDecomposition.jl` building blocks** (numerical, physics-agnostic):
- [`fit_wave_components`](@ref) — least-squares fit of `N ≥ 2` probe samples
  to two counter-propagating plane waves. Generalizes the classic 3-gauge
  Goda & Suzuki (1976) method to an arbitrary number of gauges.
- [`suggest_probe_offsets`](@ref) / [`probe_positions`](@ref) — probe-layout
  helpers that avoid the spacings that make the fit singular.
- [`sample_probe_line`](@ref) — evaluate a Gridap field at a line of probes.

**`Coefficients.jl` physical layer** (built on the above):
- [`reflection_transmission_coefficients`](@ref) — the main entry point:
  `field`/vectors + probe points + `k` → `R`, `T`, energy-balance `A`.
- [`absorption_coefficient`](@ref) — the `1 - |R|² - |T|²` energy relation
  in isolation.
- [`group_velocity`](@ref), [`energy_flux`](@ref),
  [`absorption_from_dissipated_power`](@ref),
  [`resonator_dissipated_power`](@ref) — an independent, power-based route
  to `A` for structures with an explicit dissipative element (e.g. a
  `ResonatorSingle` PTO), useful as a cross-check on the energy-balance `A`.

## Quick start

```julia
import HydroElasticFEM.PostProcessing as PP

# ϕh, κh, ηh = result.solution  (from a HydroElasticFEM.Simulation solve)
# k = sea_state.k[1]            (the wavenumber the solve itself used)

Lwave  = 2π / k
xs_in  = PP.probe_positions(xb0, :upwave,   Lwave, H0)   # upstream of the structure
xs_out = PP.probe_positions(xb1, :downwave, Lwave, H0)   # downstream of the structure

rt = PP.reflection_transmission_coefficients(κh, xs_in, xs_out, k; y=H0, η0=η0)
println(rt)                            # |R|=..., |T|=..., |R|²+|T|²=..., A=...
```

See the "Reflection, Transmission & Absorption" guide page for a full worked
example on a `TankDomain`, including the vertical-coordinate convention
gotcha (`y` above) and how to cross-check `A` against a `ResonatorSingle`'s
dissipated power.
"""
module PostProcessing

using LinearAlgebra
using Printf
using Gridap: Point

include("WaveDecomposition.jl")
include("Coefficients.jl")

export WaveProbeFit, fit_wave_components
export suggest_probe_offsets, probe_positions, sample_probe_line
export ReflectionTransmissionResult, reflection_transmission_coefficients
export absorption_coefficient
export group_velocity, energy_flux
export absorption_from_dissipated_power, resonator_dissipated_power

end # module PostProcessing
