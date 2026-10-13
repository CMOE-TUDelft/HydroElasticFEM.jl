# ─────────────────────────────────────────────────────────────
# Precomputed spectral forcing for IncidentSea-driven problems
# ─────────────────────────────────────────────────────────────
#
# When every time-dependent input of a problem is a field of one
# `P.IncidentSea`, the right-hand side is exactly
#
#     F(t) = F₀ + ramp(t) Σ_g [cos(ω_g t) F_c,g + sin(ω_g t) F_s,g]
#
# (see IncidentSea.jl).  The vectors are assembled once, with the sea
# switched to `:off` (F₀) and to the cos/sin part of each frequency group,
# and stored as sparse (index, value) vectors: the forcing only touches
# boundary DOFs.

# Non-zero entries of an assembled load vector.
struct _CompactVector
    idx::Vector{Int}
    val::Vector{Float64}
end
function _CompactVector(v::AbstractVector)
    idx = findall(!iszero, v)
    _CompactVector(idx, Float64.(v[idx]))
end

"""
    SpectralForcing

Forcing `F(t) = F₀ + ramp(t) Σ_g [cos(ω_g t) F_c,g + sin(ω_g t) F_s,g]`
with precomputed sparse load vectors.
"""
struct SpectralForcing
    F0::_CompactVector
    ω::Vector{Float64}
    Fc::Vector{_CompactVector}
    Fs::Vector{_CompactVector}
    ramp::Function
end

# r .-= F(t)
function _subtract_forcing!(r::AbstractVector, sf::SpectralForcing, t::Real)
    _axpy_sparse!(r, -1.0, sf.F0)
    a = sf.ramp(t)
    iszero(a) && return r
    for (ω, Fc, Fs) in zip(sf.ω, sf.Fc, sf.Fs)
        s, c = sincos(ω * t)
        _axpy_sparse!(r, -a * c, Fc)
        _axpy_sparse!(r, -a * s, Fs)
    end
    r
end

function _axpy_sparse!(r::AbstractVector, α::Real, x::_CompactVector)
    nz, vals = x.idx, x.val
    @inbounds for j in eachindex(nz)
        r[nz[j]] += α * vals[j]
    end
    r
end

# True if `f` has a method taking exactly `n` positional arguments.
_takes_args(f, n::Int) = any(m -> m.nargs == n + 1 && !m.isva, methods(f))

# A generic user function is time dependent if it is `(x, t) -> …` or
# `t -> (x -> …)`.
function _is_time_dependent(f::Function)
    _takes_args(f, 2) && return true
    _takes_args(f, 1) || return false
    try
        return f(0.0) isa Function
    catch
        return false
    end
end

# Collect the IncidentSea objects referenced by the entities' inputs (fields
# of the entities and of their boundary-condition vectors).  Returns
# `nothing` if some time-dependent input is not an IncidentSea field.
function _incident_seas(entities)
    seas = P.IncidentSea[]
    function visit(v)
        if v isa P.IncidentSeaField
            s = P.incident_sea(v)
            any(x -> x === s, seas) || push!(seas, s)
            return true
        elseif v isa Function
            return !_is_time_dependent(v)
        end
        return true
    end
    for e in entities, name in fieldnames(typeof(e))
        v = getfield(e, name)
        items = v isa AbstractVector ? v : (v,)
        for item in items
            if item isa P.AbstractPotentialFlowBC
                for bname in fieldnames(typeof(item))
                    visit(getfield(item, bname)) || return nothing
                end
            else
                visit(item) || return nothing
            end
        end
    end
    return seas
end

"""
    build_spectral_forcing(entities, l, Y, assembler) -> Union{SpectralForcing, Nothing}

Precompute the forcing of the linear form `l(t, y)` when all time
dependence comes from a single `IncidentSea`; `nothing` otherwise.
"""
function build_spectral_forcing(entities, l, Y, assembler)
    seas = _incident_seas(entities)
    (seas === nothing || length(seas) != 1) && return nothing
    sea = only(seas)
    vec_of() = _CompactVector(assemble_vector(y -> l(0.0, y), assembler, Y))
    F0 = P.with_sea_mode(vec_of, sea, :off)
    ng = P.num_frequency_groups(sea)
    Fc = [P.with_sea_mode(vec_of, sea, :cos, g) for g in 1:ng]
    Fs = [P.with_sea_mode(vec_of, sea, :sin, g) for g in 1:ng]
    SpectralForcing(F0, copy(sea.ω_groups), Fc, Fs, sea.ramp)
end
