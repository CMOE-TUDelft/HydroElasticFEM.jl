"""
    ResonatorSingle <: PhysicsParameters

Parameters for a single locally resonant mass-spring-damper.

# Fields
- `M::Float64` — Mass
- `K::Float64` — Stiffness
- `C::Float64` — Damping
- `ρw::Float64` — Density of water
- `η_sym::Symbol` — Symbol of the structural field this resonator couples to
- `XZ::VectorValue{2,Float64}` — Position
- `ωn1::Float64` — Natural frequency (derived: `√(K/M)`)
"""
@with_kw struct ResonatorSingle <: PhysicsParameters
    M::Float64
    K::Float64
    C::Float64     = 0.0
    ρw::Float64    = 1025.0
    η_sym::Symbol  = :η_m
    XZ::VectorValue{2,Float64} = VectorValue(0.0, 0.0)
    ωn1::Float64   = sqrt(K / M)
end

function print_parameters(resn::ResonatorSingle)
    @printf("\n[MSG] Resonator Properties:\n")
    @printf("[VAL] M = %.4f kg\n", resn.M)
    @printf("[VAL] K = %.4f N/m\n", resn.K)
    @printf("[VAL] C = %.4f Ns/m\n", resn.C)
    @printf("[VAL] XZ = (%.4f, %.4f) m\n", resn.XZ[1], resn.XZ[2])
    @printf("[VAL] ωn1 = %.4f rad/s\n", resn.ωn1)
    println()
end

function print_parameters(resn::Vector{ResonatorSingle})
    print_parameters.(resn)
end

"""
    resonator_array(N, M::Real, K::Real, C::Real, XZ; ρw=1025.0, η_sym=:η_m) -> Vector{ResonatorSingle}

Create `N` identical resonators at positions `XZ`.
"""
function resonator_array(N::Int, M::Real, K::Real, C::Real,
                         XZ::Vector{VectorValue{2,Float64}};
                         ρw::Real=1025.0, η_sym::Symbol=:η_m)
    length(XZ) == N || throw(ArgumentError("XZ must be of length N"))
    [ResonatorSingle(M=M, K=K, C=C, ρw=ρw, η_sym=η_sym, XZ=xz) for xz in XZ]
end

"""
    resonator_array(N, M::Vector, K::Vector, C::Vector, XZ; ρw=1025.0, η_sym=:η_m) -> Vector{ResonatorSingle}

Create `N` resonators with individual parameters.
"""
function resonator_array(N::Int, M::Vector{<:Real}, K::Vector{<:Real},
                         C::Vector{<:Real},
                         XZ::Vector{VectorValue{2,Float64}};
                         ρw::Real=1025.0, η_sym::Symbol=:η_m)
    (length(M) == N && length(K) == N &&
     length(C) == N && length(XZ) == N) ||
        throw(ArgumentError("M, K, C, and XZ must be of length N"))
    [ResonatorSingle(M=m, K=k, C=c, ρw=ρw, η_sym=η_sym, XZ=xz) for (m, k, c, xz) in zip(M, K, C, XZ)]
end

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only q_i terms — no coupling to structure η

function mass(resn::Vector{ResonatorSingle}, dom::WeakFormDomains, x_tt, y)
    δ_p = dom[:δ_p]
    ξ1  = y[Symbol("q_1")]
    q1  = x_tt[Symbol("q_1")]
    val = ∫((ξ1 ⋅ q1) * 0.0)dom[:dΩ]
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qₜₜi = x_tt[Symbol("q_$i")]
        ξi   = y[Symbol("q_$i")]
        val += ri.M * δi(qₜₜi ⋅ ξi)
    end
    return val
end

function damping(resn::Vector{ResonatorSingle}, dom::WeakFormDomains, x_t, y)
    δ_p = dom[:δ_p]
    ξ1  = y[Symbol("q_1")]
    q1  = x_t[Symbol("q_1")]
    val = ∫((ξ1 ⋅ q1) * 0.0)dom[:dΩ]
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qₜi = x_t[Symbol("q_$i")]
        ξi  = y[Symbol("q_$i")]
        val += ri.C * δi(qₜi ⋅ ξi)
    end
    return val
end

function stiffness(resn::Vector{ResonatorSingle}, dom::WeakFormDomains, x, y)
    δ_p = dom[:δ_p]
    ξ1  = y[Symbol("q_1")]
    q1  = x[Symbol("q_1")]
    val = ∫((ξ1 ⋅ q1) * 0.0)dom[:dΩ]
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qi = x[Symbol("q_$i")]
        ξi = y[Symbol("q_$i")]
        val += ri.K * δi(qi ⋅ ξi)
    end
    return val
end

function rhs(resn::Vector{ResonatorSingle}, dom::WeakFormDomains, f, y)
    δ_p = dom[:δ_p]
    ξ1  = y[Symbol("q_1")]
    q1  = f[Symbol("q_1")]
    val = ∫((ξ1 ⋅ q1) * 0.0)dom[:dΩ]
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        fi = f[Symbol("q_$i")]
        ξi = y[Symbol("q_$i")]
        val += δi(fi ⋅ ξi)
    end
    return val
end
