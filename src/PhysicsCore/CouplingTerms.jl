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

function mass(pf::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, x_tt, y)
    w = y[variable_symbol(pf)]
    ∫(0.0 * w)dom[:dΓ_s]
end

function damping(pf::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, x_t, y)
    ϕ_sym = variable_symbol(pf)
    η_sym = variable_symbol(s)
    ϕₜ = x_t[ϕ_sym];  ηₜ = x_t[η_sym]
    w  = y[ϕ_sym];     v  = y[η_sym]
    ∫(v * ϕₜ - w * ηₜ)dom[:dΓ_s]
end

function stiffness(pf::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, x, y)
    w = y[variable_symbol(pf)]
    ∫(0.0 * w)dom[:dΓ_s]
end

function rhs(pf::PotentialFlow, s::AbstractStructure, dom::WeakFormDomains, f, y)
    w = y[variable_symbol(pf)]
    ∫(0.0 * w)dom[:dΓ_s]
end

# =========================================================================
# Fluid ↔ FreeSurface coupling (linearised free-surface BC on Γ_fs)
#
# Frequency domain:
#   ∫( βₕ*(u + αₕ*w)*(g*κ − iω*ϕ) + iω*w*κ ) dΓ_fs
#   with αₕ = −iω/g * (1−βₕ)/βₕ
#
# Decomposed into mass/damping/stiffness so that the composed
# weakform = −ω²·mass + (−iω)·damping + stiffness reproduces the above.
# =========================================================================

function mass(pf::PotentialFlow, fs::FreeSurface, dom::WeakFormDomains, x_tt, y)
    ϕ_sym = variable_symbol(pf)
    ϕₜₜ = x_tt[ϕ_sym]
    w    = y[ϕ_sym]
    βₕ   = fs.βₕ
    g    = fs.g
    ∫((1 - βₕ) / g * w * ϕₜₜ)dom[:dΓ_fs]
end

function damping(pf::PotentialFlow, fs::FreeSurface, dom::WeakFormDomains, x_t, y)
    ϕ_sym = variable_symbol(pf)
    κ_sym = variable_symbol(fs)
    ϕₜ = x_t[ϕ_sym]
    κₜ = x_t[κ_sym]
    w  = y[ϕ_sym]
    u  = y[κ_sym]
    βₕ = fs.βₕ
    ∫(βₕ * u * ϕₜ - βₕ * w * κₜ)dom[:dΓ_fs]
end

function stiffness(pf::PotentialFlow, fs::FreeSurface, dom::WeakFormDomains, x, y)
    κ_sym = variable_symbol(fs)
    κ = x[κ_sym]
    u = y[κ_sym]
    βₕ = fs.βₕ
    g  = fs.g
    ∫(βₕ * g * u * κ)dom[:dΓ_fs]
end

function rhs(pf::PotentialFlow, fs::FreeSurface, dom::WeakFormDomains, f, y)
    u = y[variable_symbol(fs)]
    ∫(0.0 * u)dom[:dΓ_fs]
end

# =========================================================================
# Resonator ↔ Structure coupling
#
# Cross-terms between q_i DOFs and the structural displacement η.
# Uses variable_symbol(s) to identify which structure field to couple to.
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
    η_sym = variable_symbol(s)
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
    η_sym = variable_symbol(s)
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
             dom::WeakFormDomains, f, y)
    ξ1 = y[Symbol("q_1")]
    ∫((ξ1 ⋅ ξ1) * 0.0)dom[:dΩ]
end
