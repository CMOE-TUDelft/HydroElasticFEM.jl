# ─────────────────────────────────────────────────────────────
# Algebraic residual for linear transient operators
# ─────────────────────────────────────────────────────────────
#
# For a `TransientLinearFEOperator` with constant forms, Gridap caches the
# form matrices A_k (stiffness, damping, mass) in `odeopcache.const_forms`
# but its `residual!` still re-assembles every bilinear form cell by cell at
# each time step:
#
#     r = Σ_k form_k(t, ∂tᵏu, v) - res(t, v)
#
# `AlgebraicLinearTFEOperator` wraps such an operator and computes the same
# residual algebraically,
#
#     r = Σ_k A_k ∂tᵏu - F(t),
#
# so only the forcing vector F(t) is assembled per step.  This is exact when
# all forms are constant and every Dirichlet value is zero (the FE function
# then equals its free values); otherwise the stock Gridap path is used.

const _ODEs = Gridap.ODEs

"""
    AlgebraicLinearTFEOperator(op::TransientFEOperator; forcing = nothing)

Wrapper around a linear `TransientFEOperator` (e.g. a
`TransientLinearFEOperator` with `constant_forms = (true, …)`) whose
algebraic operator evaluates the residual as `Σ_k A_k ∂tᵏu - F(t)` from the
cached constant form matrices, assembling only the forcing `F(t)` at each
step.  Falls back to Gridap's full re-assembly when a form is not constant
or a Dirichlet value is non-zero.  Everything else is delegated to `op`.

`forcing` may be a precomputed [`SpectralForcing`](@ref) that replaces the
per-step assembly of `F(t)` (see [`build_spectral_forcing`](@ref)).
"""
struct AlgebraicLinearTFEOperator{T<:_ODEs.AbstractLinearODE,O,F} <: TransientFEOperator{T}
    op::O
    forcing::F
    function AlgebraicLinearTFEOperator(op::TransientFEOperator{T}; forcing = nothing) where {T<:_ODEs.AbstractLinearODE}
        new{T,typeof(op),typeof(forcing)}(op, forcing)
    end
end

Gridap.FESpaces.get_test(w::AlgebraicLinearTFEOperator) = Gridap.FESpaces.get_test(w.op)
Gridap.FESpaces.get_trial(w::AlgebraicLinearTFEOperator) = Gridap.FESpaces.get_trial(w.op)
Gridap.Polynomials.get_order(w::AlgebraicLinearTFEOperator) = Gridap.Polynomials.get_order(w.op)
_ODEs.get_res(w::AlgebraicLinearTFEOperator) = _ODEs.get_res(w.op)
_ODEs.get_jacs(w::AlgebraicLinearTFEOperator) = _ODEs.get_jacs(w.op)
_ODEs.get_forms(w::AlgebraicLinearTFEOperator) = _ODEs.get_forms(w.op)
_ODEs.get_num_forms(w::AlgebraicLinearTFEOperator) = _ODEs.get_num_forms(w.op)
_ODEs.is_form_constant(w::AlgebraicLinearTFEOperator, k::Integer) = _ODEs.is_form_constant(w.op, k)
_ODEs.get_assembler(w::AlgebraicLinearTFEOperator) = _ODEs.get_assembler(w.op)
_ODEs.allocate_tfeopcache(w::AlgebraicLinearTFEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}) =
    _ODEs.allocate_tfeopcache(w.op, t, us)
_ODEs.update_tfeopcache!(c, w::AlgebraicLinearTFEOperator, t::Real) = _ODEs.update_tfeopcache!(c, w.op, t)

Gridap.FESpaces.get_algebraic_operator(w::AlgebraicLinearTFEOperator) =
    AlgebraicLinearODEOperator(_ODEs.ODEOpFromTFEOp(w), w.forcing)

"""
    AlgebraicLinearODEOperator

`ODEOperator` returned by `get_algebraic_operator(::AlgebraicLinearTFEOperator)`:
delegates to Gridap's `ODEOpFromTFEOp`, except for `residual!`.
"""
struct AlgebraicLinearODEOperator{T<:_ODEs.AbstractLinearODE,O,F} <: _ODEs.ODEOperator{T}
    inner::O
    forcing::F
    function AlgebraicLinearODEOperator(inner::_ODEs.ODEOperator{T}, forcing = nothing) where {T<:_ODEs.AbstractLinearODE}
        new{T,typeof(inner),typeof(forcing)}(inner, forcing)
    end
end

Gridap.Polynomials.get_order(o::AlgebraicLinearODEOperator) = Gridap.Polynomials.get_order(o.inner)
_ODEs.get_num_forms(o::AlgebraicLinearODEOperator) = _ODEs.get_num_forms(o.inner)
_ODEs.get_forms(o::AlgebraicLinearODEOperator) = _ODEs.get_forms(o.inner)
_ODEs.is_form_constant(o::AlgebraicLinearODEOperator, k::Integer) = _ODEs.is_form_constant(o.inner, k)
_ODEs.allocate_odeopcache(o::AlgebraicLinearODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}) =
    _ODEs.allocate_odeopcache(o.inner, t, us)
_ODEs.update_odeopcache!(c, o::AlgebraicLinearODEOperator, t::Real) = _ODEs.update_odeopcache!(c, o.inner, t)
Gridap.Algebra.allocate_residual(o::AlgebraicLinearODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}, c) =
    Gridap.Algebra.allocate_residual(o.inner, t, us, c)
Gridap.Algebra.allocate_jacobian(o::AlgebraicLinearODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}, c) =
    Gridap.Algebra.allocate_jacobian(o.inner, t, us, c)
function _ODEs.jacobian_add!(J::AbstractMatrix, o::AlgebraicLinearODEOperator, t::Real,
                             us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}}, c)
    _ODEs.jacobian_add!(J, o.inner, t, us, ws, c)
end

function Gridap.Algebra.residual!(r::AbstractVector, o::AlgebraicLinearODEOperator, t::Real,
                                  us::Tuple{Vararg{AbstractVector}}, odeopcache; add::Bool = false)
    if !_algebraic_residual_applies(odeopcache)
        return Gridap.Algebra.residual!(r, o.inner, t, us, odeopcache; add = add)
    end
    !add && fill!(r, zero(eltype(r)))

    # Forcing: residual = Σ_k A_k ∂tᵏu - F(t)
    if o.forcing === nothing
        tfeop = o.inner.tfeop
        V = Gridap.FESpaces.get_test(tfeop)
        v = get_fe_basis(V)
        uh = _ODEs._make_uh_from_us(o.inner, us, odeopcache.Us)
        dc = (-1) * _ODEs.get_res(tfeop)(t, uh, v)
        vecdata = Gridap.FESpaces.collect_cell_vector(V, dc)
        Gridap.FESpaces.assemble_vector_add!(r, _ODEs.get_assembler(tfeop), vecdata)
    else
        _subtract_forcing!(r, o.forcing, t)
    end

    for (A, u) in zip(odeopcache.const_forms, us)
        mul!(r, A, u, true, true)
    end
    r
end

# The algebraic residual equals the assembled one when every form matrix is
# cached (constant forms) and all Dirichlet values are zero.
function _algebraic_residual_applies(odeopcache)
    all(!isnothing, odeopcache.const_forms) || return false
    all(_has_zero_dirichlet, odeopcache.Us)
end

_has_zero_dirichlet(U::MultiFieldFESpace) = all(_has_zero_dirichlet, U.spaces)
function _has_zero_dirichlet(U::SingleFieldFESpace)
    num_dirichlet_dofs(U) == 0 && return true
    iszero(get_dirichlet_dof_values(U))
end
_has_zero_dirichlet(U) = false
