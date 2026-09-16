# ─────────────────────────────────────────────────────────────
# Generic complex-exponential probe-line decomposition.
#
# This file has no knowledge of reflection/transmission physics — it just
# fits `N ≥ 2` samples of a complex, single-frequency field along a straight,
# constant-depth line to a sum of two counter-propagating plane waves. This
# is the numerical core shared by every coefficient in Coefficients.jl.
# ─────────────────────────────────────────────────────────────

"""
    WaveProbeFit

Result of decomposing `N ≥ 2` free-surface samples, taken along a line at
constant depth, into two counter-propagating plane-wave components.

# Fields
- `A::ComplexF64` — complex amplitude of the `+x`-travelling component.
- `B::ComplexF64` — complex amplitude of the `-x`-travelling component.
- `cond::Float64` — condition number of the fit's design matrix. Large
  values (rule of thumb: `> 50`) indicate the probe spacings are close to a
  `k`-dependent singularity (a spacing at or near a multiple of half the
  local wavelength) and the fit should not be trusted.
- `residual::Float64` — relative least-squares residual,
  `norm(M*[A;B] - vals) / norm(vals)`. Non-zero even for exact FE data
  (interpolation + evanescent-mode leakage); values `≳ 1e-2` usually mean
  the probes are too close to the structure or the mesh is under-resolved.
- `xs::Vector{Float64}` — probe x-coordinates used for the fit.
- `vals::Vector{ComplexF64}` — sampled complex values used for the fit.
- `k::Float64` — wavenumber used for the fit.

See also [`fit_wave_components`](@ref).
"""
struct WaveProbeFit
    A::ComplexF64
    B::ComplexF64
    cond::Float64
    residual::Float64
    xs::Vector{Float64}
    vals::Vector{ComplexF64}
    k::Float64
end

function Base.show(io::IO, f::WaveProbeFit)
    @printf(io, "WaveProbeFit(A=%.4g∠%.1f°, B=%.4g∠%.1f°, cond=%.1f, residual=%.2e, N=%d)",
        abs(f.A), rad2deg(angle(f.A)), abs(f.B), rad2deg(angle(f.B)), f.cond, f.residual, length(f.xs))
end

"""
    fit_wave_components(xs, vals, k) -> WaveProbeFit

Least-squares fit of

    vals[i] ≈ A * exp(im*k*xs[i]) + B * exp(-im*k*xs[i])

over `N ≥ 2` probes at constant depth, for a single, known wavenumber `k`.

This is the generalized Goda & Suzuki (1976) decomposition: for `N = 3` it
reduces exactly to their classic three-gauge method; for `N > 3` it is the
natural least-squares extension (cf. Zelt & Skjelbreia, 1992), which damps
out evanescent-mode leakage and finite-element discretisation noise instead
of relying on a single gauge triplet. `N = 2` is the minimal (exactly
determined) case and has no averaging benefit.

`xs` and `vals` need not be sorted or equally spaced, but see
[`suggest_probe_offsets`](@ref) for spacing guidance: some spacings make
the fit singular (or numerically close to it) regardless of `N`.

# Arguments
- `xs::AbstractVector{<:Real}`: probe x-coordinates.
- `vals::AbstractVector{<:Number}`: sampled complex field values at `xs`
  (real-valued input is accepted and promoted to `ComplexF64`).
- `k::Real`: wavenumber of the plane waves being fitted.

# Returns
A [`WaveProbeFit`](@ref).

# Example
```julia
xs   = [10.0, 13.7, 18.2, 25.9]
vals = κh.(Point.(xs, H0))          # sample a Gridap CellField
fit  = fit_wave_components(xs, vals, k)
fit.cond > 50 && @warn "ill-conditioned probe layout" fit.cond
```

# References
- Goda, Y., & Suzuki, Y. (1976). Estimation of incident and reflected waves
  in random wave experiments. *Proc. 15th Int. Conf. Coastal Eng.*, 828-845.
- Zelt, J. A., & Skjelbreia, J. E. (1992). Estimating incident and reflected
  wave fields using an arbitrary number of wave gauges.
  *Proc. 23rd Int. Conf. Coastal Eng.*, 777-789.
"""
function fit_wave_components(xs::AbstractVector{<:Real},
                              vals::AbstractVector{<:Number},
                              k::Real)
    N = length(xs)
    N == length(vals) || error("xs and vals must have the same length (got $N and $(length(vals))).")
    N >= 2 || error("Need at least 2 probes to fit A and B; got $N.")

    M = Matrix{ComplexF64}(undef, N, 2)
    for i in 1:N
        M[i, 1] = exp(im * k * xs[i])
        M[i, 2] = exp(-im * k * xs[i])
    end
    valsC = collect(ComplexF64, vals)

    sol = M \ valsC   # exact solve for N=2; least squares via QR for N>2
    A, B = sol[1], sol[2]

    resid_abs = norm(M * sol - valsC)
    resid_rel = resid_abs / max(norm(valsC), eps())

    return WaveProbeFit(A, B, cond(M), resid_rel, collect(Float64, xs), valsC, Float64(k))
end

"""
    suggest_probe_offsets(Lwave; n=4) -> Vector{Float64}

Return `n` increasing probe offsets, in metres, measured from a reference
point, spread over less than one wavelength `Lwave = 2π/k`.

The fractions of `Lwave` are deliberately irregular (`0.10, 0.20, 0.32,
0.46, 0.63, 0.85`) so that no pairwise spacing among them — nor any
individual offset — lands near a multiple of `Lwave/2`, which is exactly
where [`fit_wave_components`](@ref)'s design matrix becomes singular
(`sin(kΔx) → 0`). Goda & Suzuki (1976) recommend keeping at least one gauge
spacing between `0.05` and `0.45` of the wavelength; these fractions satisfy
that for any `n`.

`n` must be between 2 and 6 (the length of the built-in fraction set). For
more probes, call this twice with different `n` and concatenate, or build a
custom offset vector directly — [`fit_wave_components`](@ref) accepts any
probe layout, ordered or not.

# Example
```julia
Lwave = 2π / k
offs  = suggest_probe_offsets(Lwave; n=4)   # 4 offsets, increasing
xs_up = (x_structure_edge - standoff) .- offs
```
"""
function suggest_probe_offsets(Lwave::Real; n::Int=4)
    n in 2:6 || error("n must be between 2 and 6 (got $n); see the docstring for how to extend this.")
    fracs = (0.10, 0.20, 0.32, 0.46, 0.63, 0.85)
    return collect(fracs[1:n]) .* Lwave
end

"""
    probe_positions(edge_x, direction, Lwave, H; margin=1.5, n=4) -> Vector{Float64}

Build a full set of far-field probe x-coordinates on one side of a
structure, combining an evanescent-mode standoff with
[`suggest_probe_offsets`](@ref).

# Arguments
- `edge_x::Real`: x-coordinate of the structure edge facing this side (the
  upstream edge for the incident/up-wave side, the downstream edge for the
  transmitted/down-wave side).
- `direction::Symbol`: `:upwave` places probes at `x < edge_x`, walking away
  from the structure; `:downwave` places them at `x > edge_x`.
- `Lwave::Real`: wavelength, `2π/k`.
- `H::Real`: still-water depth, used to set the evanescent-mode standoff
  (evanescent modes decay over a length scale of order `H`).

# Keyword arguments
- `margin::Real=1.5`: standoff from `edge_x`, in units of `H`.
- `n::Int=4`: number of probes (`2`–`6`, see [`suggest_probe_offsets`](@ref)).

# Returns
A `Vector{Float64}` of probe x-coordinates, ordered by increasing distance
from `edge_x`. Callers are responsible for checking the result stays inside
the fluid domain (a short domain and a long wavelength can push probes past
the tank boundary).

# Example
```julia
Lwave  = 2π / k
xs_in  = probe_positions(xb0, :upwave,   Lwave, H0)
xs_out = probe_positions(xb1, :downwave, Lwave, H0)
```
"""
function probe_positions(edge_x::Real, direction::Symbol, Lwave::Real, H::Real;
                          margin::Real=1.5, n::Int=4)
    direction in (:upwave, :downwave) ||
        error("direction must be :upwave or :downwave (got :$direction).")
    standoff = margin * H
    offs = suggest_probe_offsets(Lwave; n=n)
    return direction == :upwave ? (edge_x - standoff) .- offs : (edge_x + standoff) .+ offs
end

"""
    sample_probe_line(field, xs, y) -> Vector{ComplexF64}

Evaluate a Gridap `CellField`/`FEFunction` (or any callable accepting a
`Gridap.Point`) at points `(xs[i], y)`, returning a plain `ComplexF64`
vector suitable for [`fit_wave_components`](@ref).

`y` must match the actual elevation of the surface `field` lives on for your
mesh/domain convention. In particular, `TankDomain`'s default (unshifted)
Cartesian mesh has the free surface at `z = H0`, *not* `z = 0` — see the
coordinate-convention note in [`reflection_transmission_coefficients`](@ref).
"""
function sample_probe_line(field, xs::AbstractVector{<:Real}, y::Real)
    return ComplexF64[ComplexF64(field(Point(x, y))) for x in xs]
end
