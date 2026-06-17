# Run Your First Hydroelastic Simulation

This walkthrough gets you from a fresh Julia install to a working frequency-domain
hydroelastic simulation in about 15 minutes.
No prior Gridap experience is required, but basic Julia familiarity is assumed.

## 1. Installation

Start Julia 1.10 or newer and activate a project environment:

```julia
using Pkg
Pkg.activate(".")                          # or Pkg.activate("path/to/myproject")
Pkg.add("HydroElasticFEM")
```

HydroElasticFEM pulls in [Gridap.jl](https://github.com/gridap/Gridap.jl) and
[WaveSpec.jl](https://github.com/CMOE/WaveSpec.jl) automatically.
After the install, test that the package loads:

```julia
using HydroElasticFEM
```

If you see `WaveSpec initialized: Ready for spectral sea state synthesis.`, you are good.

## 2. Empty Tank: Potential Flow with a Free Surface

The simplest meaningful simulation is a 2D rectangular tank driven by a monochromatic
incident wave at the inlet.
There is no structure — just fluid and a free surface.

Copy the following script into a file or run it in the REPL:

```julia
using HydroElasticFEM
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S
using WaveSpec
using Gridap
using Printf
using LinearAlgebra

# ── 1. Geometry ──────────────────────────────────────────────────────────
# TankDomain builds a 2D Cartesian mesh.
# L = tank length [m], H = still-water depth [m],
# nx, ny = number of cells in x and y.
tank = TankDomain(L = 60.0, H = 10.0, nx = 120, ny = 10)

# ── 2. Incident wave ─────────────────────────────────────────────────────
ω  = 2.0           # angular frequency [rad/s]
H0 = 10.0          # depth [m]  — must match TankDomain H
η0 = 0.25          # wave amplitude [m]
k  = WaveSpec.AiryWaves.solve_wavenumber(ω, H0)   # wave number [1/m]

# Airy-wave functions needed for the inlet traction BC.
ϕin(x) = -im * (η0 * ω / k) * (cosh(k * x[2]) / sinh(k * H0)) * exp(im * k * x[1])
vin(x) = VectorValue(
    (η0 * ω) * (cosh(k * x[2]) / sinh(k * H0)) * exp(im * k * x[1]),
    -im * ω * η0 * exp(im * k * x[1]),
)
# Inlet traction: normal velocity - ik·ϕ  (Robin condition)
f_in(x) = -(vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * k * ϕin(x)

# Build a sea-state object (carries k and ω for downstream physics)
sea_state = begin
    spec   = WaveSpec.ContinuousSpectrums.RegularWave(2η0, 2π/ω)
    ds     = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
    spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
    WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, [ω], [k], [0.0], H0, 1)
end

# ── 3. Physics ───────────────────────────────────────────────────────────
# PotentialFlow owns the bulk fluid (Laplace + BCs on inlet/outlet).
# FESpaceConfig sets the polynomial order and the field type.
# vector_type = Vector{ComplexF64} is required for frequency-domain problems.
fluid = PotentialFlow(
    ρw = 1025.0,               # water density [kg/m³]
    g  = 9.81,                 # gravity [m/s²]
    sea_state = sea_state,
    boundary_conditions = [
        RadiationBC(domain = :dΓin),
        RadiationBC(domain = :dΓout),
        PrescribedInletPotentialBC(domain = :dΓin, forcing = f_in, quantity = :traction),
    ],
    fe = FESpaceConfig(order = 1, vector_type = Vector{ComplexF64}),
    space_domain_symbol = :Ω,
)

# FreeSurface owns the free-surface elevation variable κ.
# βₕ ∈ (0,1] is the time-splitting parameter; 0.5 is a safe default.
fsurf = FreeSurface(
    ρw = 1025.0,
    g  = 9.81,
    βₕ = 0.5,
    fe = FESpaceConfig(order = 1, vector_type = Vector{ComplexF64}),
    space_domain_symbol = :Γκ,
)

# ── 4. Simulation ─────────────────────────────────────────────────────────
# FreqDomainConfig sets the driving frequency.
# build_problem assembles the monolithic FE system.
# simulate calls Gridap's linear solver.
config  = S.FreqDomainConfig(ω = ω)
problem = build_problem(tank, PhysicsParameters[fluid, fsurf], config)
result  = simulate(problem)

ϕh, κh = result.solution    # velocity potential and free-surface elevation

# ── 5. Inspect the result ─────────────────────────────────────────────────
# Probe the free-surface amplitude at several x-positions.
probes = [Point(x, 0.0) for x in range(5.0, 55.0, length = 9)]
κ_vals = κh(probes)

println("x [m]   |κ|/η₀")
for (pt, κ) in zip(probes, κ_vals)
    @printf("  %5.1f  %6.4f\n", pt[1], abs(κ) / η0)
end
```

**Expected output** (approximately, for the defaults above):

```
x [m]   |κ|/η₀
    5.0  1.0032
   11.2  0.9981
   17.5  1.0018
   23.8  0.9994
   30.0  1.0008
   36.2  0.9979
   42.5  1.0021
   48.8  0.9986
   55.0  1.0010
```

`|κ|/η₀ ≈ 1` everywhere means that the wave propagates through the tank without
significant reflection — exactly what a well-conditioned empty-tank solution should show.
Small deviations from 1 are normal numerical discretisation error.

### Common mistakes

**Wrong name for `TankDomain`:** The old names `TankDomain2D` and `TankDomain3D`
were removed. Using them gives:

```
ERROR: UndefVarError: `TankDomain2D` not defined
```

Use `TankDomain(L=..., H=..., nx=..., ny=...)` for 2D
or `TankDomain(L=..., W=..., H=..., nx=..., ny=..., nz=...)` for 3D.

**Forgetting to set `ω`:** If you omit `FreqDomainConfig(ω = ω)` and pass only
`S.FreqDomainConfig()`, you will get the default `ω = 0`, which makes the
free-surface equation degenerate. Always match `ω` with the wave you defined above.

**Real-valued FE spaces:** Passing `vector_type = Vector{Float64}` to a
frequency-domain problem causes a `DimensionMismatch` when complex-valued boundary
data is assembled. Always use `Vector{ComplexF64}` in frequency-domain runs.

**Units:** All parameters are in SI. `ρw = 1025.0` kg/m³, `g = 9.81` m/s², lengths
in meters, `ω` in rad/s. Mixing unit systems is the most common source of silently
wrong results.

## 3. Add a Floating Membrane

Now add a thin elastic membrane that covers the middle third of the free surface.
The membrane is a tension-only structural entity (no bending stiffness), coupled to
the fluid through a kinematic condition on its wetted face.

Replace the geometry and physics blocks above with:

```julia
# Membrane spans x ∈ [20, 40] m inside a 60 m tank.
xm0 = 20.0
Lm  = 20.0

tank = TankDomain(
    L = 60.0, H = 10.0, nx = 120, ny = 10,
    structure_domains = [
        StructureDomain(L = Lm, x₀ = [xm0, 0.0], domain_symbol = :Γm),
    ],
)

# PotentialFlow and FreeSurface are unchanged (omit structure_domains from free surface).
# Add the membrane entity:
membrane = Membrane(
    L  = Lm,
    mᵨ = 0.9,       # surface mass density / ρ_w [dimensionless]
    Tᵨ = 98.1,      # surface tension / ρ_w [m³/s²]
    g  = 9.81,
    fe = FESpaceConfig(order = 1, vector_type = Vector{ComplexF64}),
    space_domain_symbol = :Γη,
)

physics = PhysicsParameters[fluid, fsurf, membrane]
problem = build_problem(tank, physics, config)
result  = simulate(problem)

ϕh, κh, ηh = result.solution
```

The solution now has three fields: velocity potential `ϕh`, free-surface elevation
`κh` on the open water, and membrane displacement `ηh`.

**What does `|η|/η₀` mean physically?** It is the ratio of the membrane's displacement
amplitude to the incident wave amplitude. A value greater than 1 indicates amplification
(resonance), less than 1 indicates the membrane is attenuating the wave energy, and
exactly 1 would mean the membrane moves exactly with the undisturbed wave.

### Save VTK output for Paraview

```julia
trians = get_triangulations(problem)
writevtk(trians[:Ω],  "fluid",      cellfields = ["phi_re" => real(ϕh), "phi_im" => imag(ϕh)])
writevtk(trians[:Γκ], "free_surf",  cellfields = ["kap_re" => real(κh), "kap_im" => imag(κh)])
writevtk(trians[:Γη], "membrane",   cellfields = ["eta_re" => real(ηh), "eta_im" => imag(ηh)])
```

Open the `.vtu` files in [ParaView](https://www.paraview.org). Plot `phi_re` in the
fluid domain to see the standing-wave pattern, and `eta_re` on the membrane to see its
deflection shape.

## 4. What to Read Next

- [Architecture Guide](@ref "Architecture Guide") — explains the three-layer design,
  the trait-based entity system, and the four-step simulation pipeline in detail.
- [Adding a New Structural Entity](@ref "How to Add a New Structural Entity") — step-by-step
  recipe for implementing your own physics entity.
- [`TankDomain`](@ref), [`PotentialFlow`](@ref), [`FreeSurface`](@ref),
  [`Membrane`](@ref) — full API docstrings.
