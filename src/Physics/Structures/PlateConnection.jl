# ==========================================================================
# PlateConnection — shear and rotational connection between two plates
# ==========================================================================
#
# Two Kirchhoff-Love plates with separate fields η_a (on Γ_a) and η_b (on
# Γ_b) meet along the interface Λ_ab declared by a
# `Geometry.StructureConnection` (plus side in a, n⁺ from a to b).  The
# connection adds interface terms to the stiffness form only; it has no field
# of its own.
#
# Notation: [[η]] = η_a⁺ - η_b⁻ (deflection jump),
#           [[∂ₙη]] = ∇η_a⁺⋅n⁺ + ∇η_b⁻⋅n⁻ (slope jump across the edge).
#
# Shear (vertical force transfer):
#   spring k_s :  ∫ k_s [[v]] [[η]] dΛ
#   :rigid     :  penalty ∫ (β D_ρ / h³) [[v]] [[η]] dΛ  (β = shear_penalty)
# Rotation (moment transfer):
#   :free      :  no term (hinge)
#   spring kᵣ  :  ∫ kᵣ [[∂ₙv]] [[∂ₙη]] dΛ
#   :rigid     :  the plate's C/DG slope-continuity terms across Λ, with the
#                 fields of a on the plus side and of b on the minus side.
#
# Following the Gridap caveat of hand-built interface skeletons, constants
# (C, unit vectors) are combined with each field before its ⁺/⁻ trace is
# taken, and slopes use the side normals.

"""
    PlateConnection(; plate_a, plate_b, interface, normal,
                    shear = :rigid, rotation = :rigid, shear_penalty = 1.0e4)

Connection between two [`KirchhoffLovePlate`](@ref)s that have separate
fields, along the interface declared by a
`Geometry.StructureConnection(a, b, domain_symbol = interface, normal_symbol = normal)`
(plus side on `plate_a`).  The entity has no field of its own and contributes
only interface stiffness terms, divided by ρ like the plates.

# Keywords
- `plate_a`, `plate_b` — the connected plates (their `symbol`s give the fields;
  their `space_domain_symbol`s give `h` and must match the connection's `a`, `b`).
- `interface`, `normal` — `IntegrationDomains` keys of the interface measure and normal.
- `shear` — `:rigid` (penalty on the deflection jump, see below) or a spring
  stiffness `k_s` per unit length, divided by ρ (`0.0` = no shear transfer).
- `rotation` — `:rigid` (C/DG slope continuity, as inside a plate), `:free`
  (hinge), or a rotational spring `kᵣ` (moment per unit length per radian, divided by ρ).
- `shear_penalty` — `β` in the rigid-shear penalty `β D_ρ / h³`
  (`D_ρ`, `h` the larger bending stiffness and the smaller element size of the
  two plates).  The penalty is not variationally consistent: the deflection
  jump scales like `1/β`.

`rotation = :rigid` with `shear = :rigid` reproduces a continuous plate (to
the penalty tolerance); `rotation = :free` gives a hinge between separate plates.
"""
@with_kw struct PlateConnection <: PhysicsParameters
    plate_a::KirchhoffLovePlate
    plate_b::KirchhoffLovePlate
    interface::Symbol
    normal::Symbol
    shear::Union{Symbol,Float64} = :rigid
    rotation::Union{Symbol,Float64} = :rigid
    shear_penalty::Float64 = 1.0e4
    fe::Nothing = nothing
    function PlateConnection(plate_a, plate_b, interface, normal, shear, rotation, shear_penalty, fe)
        shear isa Real || shear === :rigid ||
            throw(ArgumentError("PlateConnection: shear must be :rigid or a stiffness, got $shear"))
        rotation isa Real || rotation in (:rigid, :free) ||
            throw(ArgumentError("PlateConnection: rotation must be :rigid, :free or a stiffness, got $rotation"))
        variable_symbol(plate_a) == variable_symbol(plate_b) &&
            throw(ArgumentError("PlateConnection: the two plates must have different field symbols"))
        new(plate_a, plate_b, interface, normal,
            shear isa Real ? Float64(shear) : shear,
            rotation isa Real ? Float64(rotation) : rotation, shear_penalty, fe)
    end
end

function print_parameters(c::PlateConnection)
    @printf("\n[MSG] Plate connection :%s ↔ :%s on :%s\n",
            variable_symbol(c.plate_a), variable_symbol(c.plate_b), c.interface)
    println("[VAL] shear = ", c.shear, ", rotation = ", c.rotation)
end

# Field-less coupling entity: no unknowns, stiffness terms only.
variable_symbol(::PlateConnection) = nothing
variable_symbols(::PlateConnection) = ()
field_fe_configs(::PlateConnection) = ()
has_mass_form(::PlateConnection) = false
has_damping_form(::PlateConnection) = false
has_stiffness_form(::PlateConnection) = true
has_rhs_form(::PlateConnection) = false

function stiffness(c::PlateConnection, dom::IntegrationDomains, x, y)
    pa, pb = c.plate_a, c.plate_b
    ηa, ηb = x[variable_symbol(pa)], x[variable_symbol(pb)]
    va, vb = y[variable_symbol(pa)], y[variable_symbol(pb)]
    dΛ, n = dom[c.interface], dom[c.normal]
    h = min(dom[skeleton_keys(pa.space_domain_symbol)[3]], dom[skeleton_keys(pb.space_domain_symbol)[3]])
    D_ρ = max(pa.C[1, 1, 1, 1], pb.C[1, 1, 1, 1])

    # Constant factors k are applied to the fields (and after ∇) before the
    # ⁺/⁻ traces are taken.
    jump_η(ua, ub, k = 1.0) = (k * ua).⁺ - (k * ub).⁻
    slope_jump(ua, ub, k = 1.0) = (k * ∇(ua)).⁺ ⋅ n.⁺ + (k * ∇(ub)).⁻ ⋅ n.⁻
    grad_jump(ua, ub, k = 1.0) = (k * ∇(ua)).⁺ - (k * ∇(ub)).⁻

    val = nothing
    # Shear
    k_s = c.shear === :rigid ? c.shear_penalty * D_ρ / h^3 : c.shear
    iszero(k_s) || (val = _add_contribution(val,
        ∫(jump_η(va, vb) * jump_η(ηa, ηb, k_s))dΛ))

    # Rotation
    if c.rotation isa Float64
        kᵣ = c.rotation
        iszero(kᵣ) || (val = _add_contribution(val,
            ∫(slope_jump(va, vb) * slope_jump(ηa, ηb, kᵣ))dΛ))
    elseif c.rotation === :rigid
        γ = max(pa.fe.γ, pb.fe.γ)
        τ = D_ρ * γ / h
        mean_M(ua, ub) = (0.5 * (pa.C ⊙ ∇∇(ua))).⁺ + (0.5 * (pb.C ⊙ ∇∇(ub))).⁻
        # (vectors: ⋅ rather than ⊙, which Gridap does not support between
        # blocks of different fields)
        val = _add_contribution(val, ∫(
            -grad_jump(va, vb) ⋅ (mean_M(ηa, ηb) ⋅ n.⁺)
            - (mean_M(va, vb) ⋅ n.⁺) ⋅ grad_jump(ηa, ηb)
            + grad_jump(va, vb) ⋅ grad_jump(ηa, ηb, τ))dΛ)
    end
    val === nothing && return ∫(jump_η(va, vb) * jump_η(ηa, ηb, 0.0))dΛ
    return val
end
