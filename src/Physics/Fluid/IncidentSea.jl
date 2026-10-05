# ==========================================================================
# IncidentSea — linear (Airy) incident wave field for time-domain forcing
# ==========================================================================

"""
    IncidentSea(ω, kx, ky, amplitude, phase, k, h; ramp = t -> 1.0, z_surface = 0.0, dim = 3, g)
    IncidentSea(realization; kwargs...)
    IncidentSea(state::WaveSpec.AiryWaves.AiryState; kwargs...)

Linear incident sea made of Airy components with angular frequency `ω`,
wavenumber vector `(kx, ky)`, amplitude and phase, in water of depth `h`.
The `AiryState` constructor follows the component order and phases used by
`WaveSpec.realize`.
"""
struct IncidentSea
    ω::Vector{Float64}
    kx::Vector{Float64}
    ky::Vector{Float64}
    A::Vector{Float64}
    φ::Vector{Float64}
    k::Vector{Float64}
    h::Float64
    e2kh::Vector{Float64}
    cϕ::Vector{Float64}
    cu::Vector{Float64}
    cv::Vector{Float64}
    cw::Vector{Float64}
    ramp::Function
    z_surface::Float64
    dim::Int
    ω_groups::Vector{Float64}
    groups::Vector{Vector{Int}}
    mode::Base.RefValue{Tuple{Symbol,Int}}
end

function IncidentSea(ω, kx, ky, amplitude, phase, k, h;
                     ramp::Function = t -> 1.0, z_surface::Real = 0.0, dim::Int = 3,
                     g::Real = WaveSpec.PhysicalConstants.g)
    n = length(ω)
    all(length(v) == n for v in (kx, ky, amplitude, phase, k)) ||
        throw(DimensionMismatch("IncidentSea: all component vectors must have the same length"))
    dim in (2, 3) || throw(ArgumentError("IncidentSea: dim must be 2 or 3"))
    ω, kx, ky, A, φ, k = (collect(Float64, v) for v in (ω, kx, ky, amplitude, phase, k))
    e2kh = @. exp(-2 * k * h)
    cϕ = @. g / (ω * (1 + e2kh))
    cw = @. ω / (1 - e2kh)
    cu = @. cw * kx / k
    cv = @. cw * ky / k
    ω_groups = unique(ω)
    groups = [findall(==(w), ω) for w in ω_groups]
    IncidentSea(ω, kx, ky, A, φ, k, Float64(h), e2kh, cϕ, cu, cv, cw, ramp,
                Float64(z_surface), dim, ω_groups, groups, Ref((:time, 0)))
end

function IncidentSea(r; kwargs...)
    hasproperty(r, :components) && hasproperty(r, :k) && hasproperty(r, :h) ||
        throw(ArgumentError("IncidentSea: expected a WaveSpec AiryRealization or an AiryState"))
    c = r.components
    IncidentSea(c.ω, c.kx, c.ky, c.amplitude, c.phase, r.k, r.h; kwargs...)
end

function IncidentSea(state::WaveSpec.AiryWaves.AiryState; kwargs...)
    nω, nθ = state.nω, state.nθ
    A = vec(permutedims(WaveSpec.AiryWaves.get_amplitudes(state)))
    φ = vec(permutedims(WaveSpec.AiryWaves.get_random_phases(state)))
    ω = repeat(state.ω, inner = nθ)
    k = repeat(state.k, inner = nθ)
    kx = k .* repeat(cos.(state.θ), outer = nω)
    ky = k .* repeat(sin.(state.θ), outer = nω)
    IncidentSea(ω, kx, ky, A, φ, k, state.h; kwargs...)
end

num_frequency_groups(sea::IncidentSea) = length(sea.groups)

@inline _horizontal(sea::IncidentSea, x) = sea.dim == 3 ? (x[1], x[2]) : (x[1], 0.0)
@inline _depth_coord(sea::IncidentSea, x) = x[sea.dim] - sea.z_surface

@inline function _profiles(sea::IncidentSea, i, z)
    e⁺ = exp(sea.k[i] * z)
    e⁻ = ifelse(iszero(e⁺), zero(e⁺), sea.e2kh[i] / e⁺)
    return e⁺ + e⁻, e⁺ - e⁻
end

@inline function _sea_sum(coef::F, sea::IncidentSea, kind::Symbol, x, t) where {F}
    mode, grp = sea.mode[]
    mode === :off && return 0.0
    xh, yh = _horizontal(sea, x)
    s = 0.0
    if mode === :time
        @inbounds for i in eachindex(sea.ω)
            θ = sea.kx[i] * xh + sea.ky[i] * yh - sea.ω[i] * t + sea.φ[i]
            s += coef(i) * (kind === :cos ? cos(θ) : sin(θ))
        end
        return sea.ramp(t) * s
    end
    @inbounds for i in sea.groups[grp]
        ψ = sea.kx[i] * xh + sea.ky[i] * yh + sea.φ[i]
        sψ, cψ = sincos(ψ)
        trig = kind === :cos ? (mode === :cos ? cψ : sψ) : (mode === :cos ? sψ : -cψ)
        s += coef(i) * trig
    end
    return s
end

"""
    IncidentSeaField{K} <: Function

An `(x, t)` callable returning one field of an `IncidentSea`.
"""
struct IncidentSeaField{K} <: Function
    sea::IncidentSea
end

incident_sea(f::IncidentSeaField) = f.sea

function (f::IncidentSeaField{:η})(x, t)
    sea = f.sea
    _sea_sum(i -> sea.A[i], sea, :cos, x, t)
end

function (f::IncidentSeaField{:ϕ})(x, t)
    sea = f.sea
    z = _depth_coord(sea, x)
    _sea_sum(i -> sea.A[i] * sea.cϕ[i] * _profiles(sea, i, z)[1], sea, :sin, x, t)
end

function (f::IncidentSeaField{:w})(x, t)
    sea = f.sea
    z = _depth_coord(sea, x)
    _sea_sum(i -> sea.A[i] * sea.cw[i] * _profiles(sea, i, z)[2], sea, :sin, x, t)
end

function (f::IncidentSeaField{:velocity})(x, t)
    sea = f.sea
    z = _depth_coord(sea, x)
    u = _sea_sum(i -> sea.A[i] * sea.cu[i] * _profiles(sea, i, z)[1], sea, :cos, x, t)
    w = _sea_sum(i -> sea.A[i] * sea.cw[i] * _profiles(sea, i, z)[2], sea, :sin, x, t)
    if sea.dim == 3
        v = _sea_sum(i -> sea.A[i] * sea.cv[i] * _profiles(sea, i, z)[1], sea, :cos, x, t)
        return VectorValue(u, v, w)
    end
    return VectorValue(u, w)
end

"""
    incident_elevation(sea) -> (x, t) -> η

Free-surface elevation of the incident sea.
"""
incident_elevation(sea::IncidentSea) = IncidentSeaField{:η}(sea)

"""
    incident_potential(sea) -> (x, t) -> ϕ

Velocity potential of the incident sea.
"""
incident_potential(sea::IncidentSea) = IncidentSeaField{:ϕ}(sea)

"""
    incident_vertical_velocity(sea) -> (x, t) -> w

Vertical fluid velocity of the incident sea.
"""
incident_vertical_velocity(sea::IncidentSea) = IncidentSeaField{:w}(sea)

"""
    incident_velocity(sea) -> (x, t) -> VectorValue

Fluid velocity vector of the incident sea.
"""
incident_velocity(sea::IncidentSea) = IncidentSeaField{:velocity}(sea)

function with_sea_mode(f::Function, sea::IncidentSea, mode::Symbol, group::Int = 0)
    old = sea.mode[]
    sea.mode[] = (mode, group)
    try
        return f()
    finally
        sea.mode[] = old
    end
end
