"""
    ResonatorSingle <: PhysicsParameters

Parameters for a single locally resonant mass-spring-damper.

# Fields
- `M::Float64` — Mass
- `K::Float64` — Stiffness
- `C::Float64` — Damping
- `XZ::VectorValue{2,Float64}` — Position
- `ωn1::Float64` — Natural frequency (derived: `√(K/M)`)
"""
@with_kw struct ResonatorSingle <: PhysicsParameters
    M::Float64
    K::Float64
    C::Float64     = 0.0
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
    resonator_array(N, M::Real, K::Real, C::Real, XZ) -> Vector{ResonatorSingle}

Create `N` identical resonators at positions `XZ`.
"""
function resonator_array(N::Int, M::Real, K::Real, C::Real,
                         XZ::Vector{VectorValue{2,Float64}})
    length(XZ) == N || throw(ArgumentError("XZ must be of length N"))
    [ResonatorSingle(M=M, K=K, C=C, XZ=xz) for xz in XZ]
end

"""
    resonator_array(N, M::Vector, K::Vector, C::Vector, XZ) -> Vector{ResonatorSingle}

Create `N` resonators with individual parameters.
"""
function resonator_array(N::Int, M::Vector{<:Real}, K::Vector{<:Real},
                         C::Vector{<:Real},
                         XZ::Vector{VectorValue{2,Float64}})
    (length(M) == N && length(K) == N &&
     length(C) == N && length(XZ) == N) ||
        throw(ArgumentError("M, K, C, and XZ must be of length N"))
    [ResonatorSingle(M=m, K=k, C=c, XZ=xz) for (m, k, c, xz) in zip(M, K, C, XZ)]
end
