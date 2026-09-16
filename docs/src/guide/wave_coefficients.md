# Reflection, Transmission & Absorption

This walkthrough computes the reflection coefficient `R`, transmission
coefficient `T`, and absorption coefficient `A` for a floating elastic beam
on a flat-bed `TankDomain`, using the `PostProcessing` module. It assumes
you can already build and solve a frequency-domain problem — see
[Run Your First Hydroelastic Simulation](@ref) if not.

`PostProcessing` never touches the solve itself. Once you have a
frequency-domain free-surface field (`κh` or `ηh` from `result.solution`),
everything below is pure post-processing: sample the field at a handful of
probe points and fit it to counter-propagating plane waves.

## 1. A coordinate-convention gotcha, up front

`TankDomain`'s default (unshifted) mesh has the **seabed at `z = 0` and the
free surface at `z = H`**, not the SWL-at-zero convention some hand-written
Airy-wave formulas assume. This matters here because every probe point you
evaluate `κh` at must actually sit on the free surface. Concretely:

```julia
tank = TankDomain(L=120.0, H=10.0, nx=360, ny=12,
    structure_domains=[StructureDomain(L=10.0, x₀=[50.0, 10.0])])   # x₀[2] = H, not 0
...
κh(Point(35.0, 10.0))   # correct: z = H = 10.0
κh(Point(35.0, 0.0))    # wrong domain for this field on this mesh
```

If you instead used a custom `map=` that shifts the mesh to the `z ∈
[-H,0]` convention, use `y=0.0` (or whatever your shift produces) instead.
Either way, `PostProcessing` takes `y` as an explicit keyword precisely
because it can't know your convention — get it from wherever you built the
`StructureDomain`/incident-wave forcing, don't guess.

## 2. Set up and solve

```julia
using HydroElasticFEM
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
import HydroElasticFEM.PostProcessing as PP
using WaveSpec, Gridap, Printf

H0, ω, η0 = 10.0, 1.2, 0.1
xb0, Lb   = 50.0, 10.0
xb1       = xb0 + Lb

# Sea state (carries k, consistent with what RadiationBC will use internally)
sea_state = let
    spec   = WaveSpec.ContinuousSpectrums.RegularWave(2η0, 2π / ω)
    ds     = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
    spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
    WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, [ω], [WaveSpec.AiryWaves.solve_wavenumber(ω, H0)], [0.0], H0, 1)
end
k = sea_state.k[1]

# Incident-wave kinematics, native TankDomain z (0 at bed, H0 at surface)
ϕin(x) = -im * (η0 * ω / k) * (cosh(k * x[2]) / sinh(k * H0)) * exp(im * k * x[1])
vin(x) = VectorValue(
    (η0 * ω) * (cosh(k * x[2]) / sinh(k * H0)) * exp(im * k * x[1]),
    -im * (η0 * ω) * (sinh(k * x[2]) / sinh(k * H0)) * exp(im * k * x[1]),
)
f_in(x) = (vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * k * ϕin(x)

tank = TankDomain(L=120.0, H=H0, nx=360, ny=12,
    structure_domains=[StructureDomain(L=Lb, x₀=[xb0, H0], domain_symbol=:Γbeam)])

fluid = PotentialFlow(g=9.81, sea_state=sea_state,
    boundary_conditions=[
        RadiationBC(domain=:dΓin), RadiationBC(domain=:dΓout),
        PrescribedInletPotentialBC(domain=:dΓin, forcing=f_in, quantity=:traction),
    ],
    fe=FESpaceConfig(order=2, vector_type=Vector{ComplexF64}), space_domain_symbol=:Ω)

fsurf = FreeSurface(g=9.81, βₕ=0.5,
    fe=FESpaceConfig(order=2, vector_type=Vector{ComplexF64}), space_domain_symbol=:Γκ)

h, order = 120.0 / 360, 2
beam = EulerBernoulliBeam(L=Lb, mᵨ=0.5, EIᵨ=1.0e4,
    fe=FESpaceConfig(order=order, vector_type=Vector{ComplexF64}, γ=order * (order - 1) / h),
    space_domain_symbol=:Γη)

config  = S.FreqDomainConfig(ω=ω)
problem = build_problem(tank, PhysicsParameters[fluid, fsurf, beam], config)
result  = simulate(problem)
ϕh, κh, ηh = result.solution
```

No sponge/damping zones are needed here: the bed is flat and the problem is
single-frequency, so `RadiationBC` is exact for the scattered field, and the
inlet traction supplies exactly the correction needed to keep the incident
wave present in the total solution. (`FloatingMembraneExample.jl` in
`examples/` uses the same trick.)

## 3. Lay out probes and compute R, T, A

```julia
Lwave = 2π / k

xs_in  = PP.probe_positions(xb0, :upwave,   Lwave, H0)   # defaults: margin=1.5H, n=4
xs_out = PP.probe_positions(xb1, :downwave, Lwave, H0)

rt = PP.reflection_transmission_coefficients(κh, xs_in, xs_out, k; y=H0, η0=η0)
println(rt)
@printf("|R| = %.4f   |T| = %.4f   |R|²+|T|² = %.4f   A = %.4f\n",
    abs(rt.R), abs(rt.T), abs2(rt.R) + abs2(rt.T), rt.A)
```

`reflection_transmission_coefficients` does three things: fits the up-wave
probe line to `incident + reflected` and the down-wave line to
`transmitted (+ any spurious back-reflection)`, using
[`fit_wave_components`](@ref HydroElasticFEM.PostProcessing.fit_wave_components)
(a generalized Goda & Suzuki, 1976 decomposition); normalises `R` and `T` by
the *fitted* incident amplitude rather than the nominal `η0` you passed in
(so they're correct even if your probes or mesh introduce a small amplitude
error); and warns if that fitted amplitude drifts more than 5% from `η0`,
which usually means the mesh needs refining or the probes are too close to
the beam.

## 4. Sanity-check with the energy balance

The beam here has no structural damping (`τ = 0`, the `EulerBernoulliBeam`
default), so energy in equals energy out:

```julia
@assert isapprox(abs2(rt.R) + abs2(rt.T), 1.0; atol=1e-2) "check probe placement / mesh resolution"
```

If this fails by more than a percent or two, first check
`rt.incident_fit.cond` and `rt.transmitted_fit.cond`
(both should be well under 50 — `reflection_transmission_coefficients`
raises an error above `cond_tol` for exactly this reason) and
`rt.transmitted_fit.B` (should be small — a large value means the down-wave
probe line is picking up a reflection from somewhere downstream, not just
the transmitted wave).

## 5. With a dissipative structure: cross-check `A` two ways

Swap the beam for a `ResonatorSingle` with `C > 0` (or give the beam
`τ > 0`) and the energy balance no longer closes to 1 — `A` measures real
absorption. You can cross-check the energy-balance `A` against a direct,
independent power calculation:

```julia
qh = result.solution[3]   # resonator response amplitude, if PhysicsParameters = [fluid, fsurf, resonator]
P_diss = PP.resonator_dissipated_power(qh(Point(x_resonator, H0)), ω, resonator.C)
A_from_power = PP.absorption_from_dissipated_power(P_diss, η0, k, H0)

@printf("A (energy balance) = %.4f   A (dissipated power) = %.4f\n", rt.A, A_from_power)
```

The two routes are independent (one from far-field wave amplitudes, one
from the near-field resonator response), so persistent disagreement between
them is a strong signal to refine the mesh or revisit probe placement
before trusting either number.

## Common mistakes

**Wrong `y`:** the single most common error — see §1. Symptoms: `κh(...)`
throws a point-location error, or returns values with no resemblance to a
plane wave (`fit.residual` will be large, or `fit.cond` may look fine while
the physical answer is nonsense — the residual is the more reliable tell).

**Probe spacing near a wavelength multiple:** `reflection_transmission_coefficients`
raises an error when `cond > cond_tol` rather than silently returning a
noisy `R`/`T`. Don't raise `cond_tol` to make the error go away — move the
probes (`PP.suggest_probe_offsets` avoids the worst spacings by
construction, but a coincidence between your `margin` and `k` can still land
badly; try a different `n` or `margin`).

**Reusing an independently-recomputed `k`:** always take `k` from the same
`sea_state` (or dispersion solve) the simulation used for its radiation
boundary condition. A `k` that's even slightly off makes
`fit_wave_components` fit the wrong plane waves — usually visible as a
large `residual` even though `cond` looks fine.

**Probes outside the domain:** `probe_positions` doesn't check tank bounds
for you (it doesn't know them). A long wavelength or a short tank can push
probes past the inlet/outlet — reduce `n` or `margin`, or extend the tank.
