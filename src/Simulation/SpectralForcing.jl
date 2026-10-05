# ─────────────────────────────────────────────────────────────
# Precomputed spectral forcing for IncidentSea-driven problems
# ─────────────────────────────────────────────────────────────

using LinearAlgebra: mul!

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

_takes_args(f, n::Int) = any(m -> m.nargs == n + 1 && !m.isva, methods(f))

function _is_time_dependent(f::Function)
    _takes_args(f, 2) && return true
    _takes_args(f, 1) || return false
    try
        return f(0.0) isa Function
    catch
        return false
    end
end

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

Precompute the forcing of `l(t, y)` when all time dependence comes from a
single `IncidentSea`; return `nothing` otherwise.
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

const _ODEs = Gridap.ODEs

"""
    SpectralForcingTFEOperator(op, forcing)

Transient operator adapter that substitutes a precomputed forcing vector
while delegating constant-form residual optimization to Gridap.
"""
struct SpectralForcingTFEOperator{T<:_ODEs.AbstractLinearODE,O,F} <: TransientFEOperator{T}
    op::O
    forcing::F
    function SpectralForcingTFEOperator(op::TransientFEOperator{T}, forcing::SpectralForcing) where {T<:_ODEs.AbstractLinearODE}
        new{T,typeof(op),typeof(forcing)}(op, forcing)
    end
end

Gridap.FESpaces.get_test(w::SpectralForcingTFEOperator) = Gridap.FESpaces.get_test(w.op)
Gridap.FESpaces.get_trial(w::SpectralForcingTFEOperator) = Gridap.FESpaces.get_trial(w.op)
Gridap.Polynomials.get_order(w::SpectralForcingTFEOperator) = Gridap.Polynomials.get_order(w.op)
_ODEs.get_res(w::SpectralForcingTFEOperator) = _ODEs.get_res(w.op)
_ODEs.get_jacs(w::SpectralForcingTFEOperator) = _ODEs.get_jacs(w.op)
_ODEs.get_forms(w::SpectralForcingTFEOperator) = _ODEs.get_forms(w.op)
_ODEs.get_num_forms(w::SpectralForcingTFEOperator) = _ODEs.get_num_forms(w.op)
_ODEs.is_form_constant(w::SpectralForcingTFEOperator, k::Integer) = _ODEs.is_form_constant(w.op, k)
_ODEs.get_assembler(w::SpectralForcingTFEOperator) = _ODEs.get_assembler(w.op)
_ODEs.allocate_tfeopcache(w::SpectralForcingTFEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}) =
    _ODEs.allocate_tfeopcache(w.op, t, us)
_ODEs.update_tfeopcache!(c, w::SpectralForcingTFEOperator, t::Real) =
    _ODEs.update_tfeopcache!(c, w.op, t)

Gridap.FESpaces.get_algebraic_operator(w::SpectralForcingTFEOperator) =
    SpectralForcingODEOperator(_ODEs.ODEOpFromTFEOp(w), w.forcing)

struct SpectralForcingODEOperator{T<:_ODEs.AbstractLinearODE,O,F} <: _ODEs.ODEOperator{T}
    inner::O
    forcing::F
    function SpectralForcingODEOperator(inner::_ODEs.ODEOperator{T}, forcing::SpectralForcing) where {T<:_ODEs.AbstractLinearODE}
        new{T,typeof(inner),typeof(forcing)}(inner, forcing)
    end
end

Gridap.Polynomials.get_order(o::SpectralForcingODEOperator) = Gridap.Polynomials.get_order(o.inner)
_ODEs.get_num_forms(o::SpectralForcingODEOperator) = _ODEs.get_num_forms(o.inner)
_ODEs.get_forms(o::SpectralForcingODEOperator) = _ODEs.get_forms(o.inner)
_ODEs.is_form_constant(o::SpectralForcingODEOperator, k::Integer) = _ODEs.is_form_constant(o.inner, k)
_ODEs.allocate_odeopcache(o::SpectralForcingODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}) =
    _ODEs.allocate_odeopcache(o.inner, t, us)
_ODEs.update_odeopcache!(c, o::SpectralForcingODEOperator, t::Real) =
    _ODEs.update_odeopcache!(c, o.inner, t)
Gridap.Algebra.allocate_residual(o::SpectralForcingODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}, c) =
    Gridap.Algebra.allocate_residual(o.inner, t, us, c)
Gridap.Algebra.allocate_jacobian(o::SpectralForcingODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}, c) =
    Gridap.Algebra.allocate_jacobian(o.inner, t, us, c)

function _ODEs.jacobian_add!(J::AbstractMatrix, o::SpectralForcingODEOperator, t::Real,
                             us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}}, c)
    _ODEs.jacobian_add!(J, o.inner, t, us, ws, c)
end

function Gridap.Algebra.residual!(r::AbstractVector, o::SpectralForcingODEOperator, t::Real,
                                  us::Tuple{Vararg{AbstractVector}}, odeopcache; add::Bool = false)
    inner = o.inner
    tfeop = inner.tfeop
    V = Gridap.FESpaces.get_test(tfeop)
    v = get_fe_basis(V)
    uh = _ODEs._make_uh_from_us(inner, us, odeopcache.Us)
    !add && fill!(r, zero(eltype(r)))

    dc = Gridap.CellData.DomainContribution()
    forms = _ODEs.get_forms(tfeop)
    ∂tkuh = uh
    for k in 0:Gridap.Polynomials.get_order(inner)
        if _can_use_cached_form(inner, odeopcache, k)
            mul!(r, odeopcache.const_forms[k + 1], us[k + 1], true, true)
        else
            dc = dc + forms[k + 1](t, ∂tkuh, v)
        end
        k < Gridap.Polynomials.get_order(inner) && (∂tkuh = ∂t(∂tkuh))
    end

    if Gridap.CellData.num_domains(dc) > 0
        vecdata = Gridap.FESpaces.collect_cell_vector(V, dc)
        Gridap.FESpaces.assemble_vector_add!(r, _ODEs.get_assembler(tfeop), vecdata)
    end
    _subtract_forcing!(r, o.forcing, t)
end

function _can_use_cached_form(odeop, cache, k::Integer)
    _ODEs.is_form_constant(odeop, k) && _has_zero_dirichlet_values(cache.Us[k + 1])
end

function _has_zero_dirichlet_values(U::Gridap.FESpaces.SingleFieldFESpace)
    Gridap.FESpaces.num_dirichlet_dofs(U) == 0 && return true
    iszero(Gridap.FESpaces.get_dirichlet_dof_values(U))
end

function _has_zero_dirichlet_values(U::Gridap.MultiField.MultiFieldFESpace)
    all(_has_zero_dirichlet_values, U.spaces)
end

_has_zero_dirichlet_values(::Any) = false
