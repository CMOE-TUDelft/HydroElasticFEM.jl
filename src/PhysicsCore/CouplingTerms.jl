# ─────────────────────────────────────────────────────────────
# Coupling weak forms between pairs of physics entities
#
# Each coupling function takes (entityA, entityB, dom, x/x_t/x_tt, y)
# and returns cross-terms that involve fields from BOTH entities.
# ─────────────────────────────────────────────────────────────

# =========================================================================
# Fluid ↔ Structure coupling (kinematic / dynamic BC on Γ_s)
#
# Time-domain coupling lives only in damping:
#     ∫( y[η] * x_t[:ϕ] − y[:ϕ] * x_t[η] ) dΓ_s
# =========================================================================

function mass(f::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, x_tt, y)
    w = y[:ϕ]
    ∫(0.0 * w)dom[:dΓ_s]
end

function damping(f::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, x_t, y)
    ηs = η_symbol(s)
    ϕₜ = x_t[:ϕ];  ηₜ = x_t[ηs]
    w  = y[:ϕ];     v  = y[ηs]
    ∫(v * ϕₜ - w * ηₜ)dom[:dΓ_s]
end

function stiffness(f::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, x, y)
    w = y[:ϕ]
    ∫(0.0 * w)dom[:dΓ_s]
end

function rhs(f::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, x, y)
    w = y[:ϕ]
    ∫(0.0 * w)dom[:dΓ_s]
end

# =========================================================================
# Resonator ↔ Structure coupling
#
# Cross-terms between q_i DOFs and the structural displacement η.
# Uses ri.η_sym to identify which structure field to couple to.
# =========================================================================

function mass(resn::Vector{ResonatorSingle}, s::AbstractStructure,
              dom::WeakFormDomains, x_tt, y)
    ξ1 = y[Symbol("q_1")]
    q1 = x_tt[Symbol("q_1")]
    ∫((ξ1 ⋅ q1) * 0.0)dom[:dΩ]
end

function damping(resn::Vector{ResonatorSingle}, s::AbstractStructure,
                 dom::WeakFormDomains, x_t, y)
    δ_p   = dom[:δ_p]
    η_sym = first(resn).η_sym
    ηₜ    = x_t[η_sym]
    v     = y[η_sym]
    î1    = VectorValue(1.0)
    ξ1    = y[Symbol("q_1")]
    q1    = x_t[Symbol("q_1")]
    val   = ∫((ξ1 ⋅ q1) * 0.0)dom[:dΩ]
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qₜi = x_t[Symbol("q_$i")]
        ξi  = y[Symbol("q_$i")]
        # force on structure from resonator velocity
        val += (ri.C / ri.ρw) * δi(v * ((qₜi ⋅ î1) - ηₜ))
        # force on resonator from structure velocity
        val += -ri.C * δi((ξi ⋅ î1) * ηₜ)
    end
    return val
end

function stiffness(resn::Vector{ResonatorSingle}, s::AbstractStructure,
                   dom::WeakFormDomains, x, y)
    δ_p   = dom[:δ_p]
    η_sym = first(resn).η_sym
    η     = x[η_sym]
    v     = y[η_sym]
    î1    = VectorValue(1.0)
    ξ1    = y[Symbol("q_1")]
    q1    = x[Symbol("q_1")]
    val   = ∫((ξ1 ⋅ q1) * 0.0)dom[:dΩ]
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qi = x[Symbol("q_$i")]
        ξi = y[Symbol("q_$i")]
        # force on structure from resonator displacement
        val += (-ri.K / ri.ρw) * δi(v * ((qi ⋅ î1) - η))
        # force on resonator from structure displacement
        val += -ri.K * δi((ξi ⋅ î1) * η)
    end
    return val
end

function rhs(resn::Vector{ResonatorSingle}, s::AbstractStructure,
             dom::WeakFormDomains, x, y)
    ξ1 = y[Symbol("q_1")]
    ∫((ξ1 ⋅ ξ1) * 0.0)dom[:dΩ]
end
