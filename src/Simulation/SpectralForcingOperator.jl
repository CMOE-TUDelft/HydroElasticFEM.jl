# ─────────────────────────────────────────────────────────────────────────────
# Spectral forcing adapter for Gridap transient linear operators
# ─────────────────────────────────────────────────────────────────────────────

const _ODEs = Gridap.ODEs

struct SpectralForcingTFEOperator{T<:_ODEs.AbstractLinearODE,O,F,Z} <:
  TransientFEOperator{T}
  op::O
  forcing::F
  zero_residual::Z
end

function SpectralForcingTFEOperator(
  op::TransientFEOperator{T}, forcing, zero_residual
) where {T<:_ODEs.AbstractLinearODE}
    wrapper_type = SpectralForcingTFEOperator{
      T,typeof(op),typeof(forcing),typeof(zero_residual)}
    wrapper_type(op, forcing, zero_residual)
end

Gridap.FESpaces.get_test(w::SpectralForcingTFEOperator) =
  Gridap.FESpaces.get_test(w.op)
Gridap.FESpaces.get_trial(w::SpectralForcingTFEOperator) =
  Gridap.FESpaces.get_trial(w.op)
Gridap.Polynomials.get_order(w::SpectralForcingTFEOperator) =
  Gridap.Polynomials.get_order(w.op)
_ODEs.get_res(w::SpectralForcingTFEOperator) = w.zero_residual
_ODEs.get_jacs(w::SpectralForcingTFEOperator) = _ODEs.get_jacs(w.op)
_ODEs.get_forms(w::SpectralForcingTFEOperator) = _ODEs.get_forms(w.op)
_ODEs.get_num_forms(w::SpectralForcingTFEOperator) = _ODEs.get_num_forms(w.op)
_ODEs.is_form_constant(w::SpectralForcingTFEOperator, k::Integer) =
  _ODEs.is_form_constant(w.op, k)
_ODEs.get_assembler(w::SpectralForcingTFEOperator) = _ODEs.get_assembler(w.op)
_ODEs.allocate_tfeopcache(
  w::SpectralForcingTFEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}
) =
  _ODEs.allocate_tfeopcache(w.op, t, us)
_ODEs.update_tfeopcache!(c, w::SpectralForcingTFEOperator, t::Real) =
  _ODEs.update_tfeopcache!(c, w.op, t)

Gridap.FESpaces.get_algebraic_operator(w::SpectralForcingTFEOperator) =
  SpectralForcingODEOperator(_ODEs.ODEOpFromTFEOp(w), w.forcing)

struct SpectralForcingODEOperator{T<:_ODEs.AbstractLinearODE,O,F} <:
  _ODEs.ODEOperator{T}
  inner::O
  forcing::F
end

function SpectralForcingODEOperator(
  inner::_ODEs.ODEOperator{T}, forcing
) where {T<:_ODEs.AbstractLinearODE}
  SpectralForcingODEOperator{T,typeof(inner),typeof(forcing)}(inner, forcing)
end

Gridap.Polynomials.get_order(o::SpectralForcingODEOperator) =
  Gridap.Polynomials.get_order(o.inner)
_ODEs.get_num_forms(o::SpectralForcingODEOperator) =
  _ODEs.get_num_forms(o.inner)
_ODEs.get_forms(o::SpectralForcingODEOperator) = _ODEs.get_forms(o.inner)
_ODEs.is_form_constant(o::SpectralForcingODEOperator, k::Integer) =
  _ODEs.is_form_constant(o.inner, k)
_ODEs.allocate_odeopcache(
  o::SpectralForcingODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}
) =
  _ODEs.allocate_odeopcache(o.inner, t, us)
_ODEs.update_odeopcache!(c, o::SpectralForcingODEOperator, t::Real) =
  _ODEs.update_odeopcache!(c, o.inner, t)
Gridap.Algebra.allocate_residual(
  o::SpectralForcingODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}, c
) =
  Gridap.Algebra.allocate_residual(o.inner, t, us, c)
Gridap.Algebra.allocate_jacobian(
  o::SpectralForcingODEOperator, t::Real, us::Tuple{Vararg{AbstractVector}}, c
) =
  Gridap.Algebra.allocate_jacobian(o.inner, t, us, c)

function _ODEs.jacobian_add!(
  J::AbstractMatrix, o::SpectralForcingODEOperator, t::Real,
  us::Tuple{Vararg{AbstractVector}}, ws::Tuple{Vararg{Real}}, c
)
  _ODEs.jacobian_add!(J, o.inner, t, us, ws, c)
end

function Gridap.Algebra.residual!(
  r::AbstractVector, o::SpectralForcingODEOperator, t::Real,
  us::Tuple{Vararg{AbstractVector}}, odeopcache; add::Bool = false
)
  Gridap.Algebra.residual!(r, o.inner, t, us, odeopcache; add = add)
  _subtract_forcing!(r, o.forcing, t)
end