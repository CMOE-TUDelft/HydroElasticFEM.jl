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

# Fluid-structure coupling only contributes through damping terms.
has_damping_form(::PotentialFlow, ::Structure) = true

function damping(pf::PotentialFlow, s::Structure, dom::IntegrationDomains, x_t, y)
    ϕ_sym = variable_symbol(pf)
    η_sym = variable_symbol(s)
    ϕₜ = x_t[ϕ_sym];  ηₜ = x_t[η_sym]
    w  = y[ϕ_sym];     v  = y[η_sym]
    ∫(v * ϕₜ - w * ηₜ)dom[:dΓη]
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

# Free-surface coupling has mass and damping terms only.
has_mass_form(::PotentialFlow, ::FreeSurface) = true
has_damping_form(::PotentialFlow, ::FreeSurface) = true
has_stiffness_form(::PotentialFlow, ::FreeSurface) = true
has_rhs_form(::PotentialFlow, ::FreeSurface) = true

function mass(pf::PotentialFlow, fs::FreeSurface, dom::IntegrationDomains, x_tt, y)
    if _simulation_domain(dom) == :time && haskey(dom, :αₕ)
        return ∫(0.0 * y[variable_symbol(pf)] * x_tt[variable_symbol(pf)])dom[:dΓκ]
    end
    ϕ_sym = variable_symbol(pf)
    ϕₜₜ = x_tt[ϕ_sym]
    w    = y[ϕ_sym]
    βₕ   = fs.βₕ
    g    = fs.g
    ∫((1 - βₕ) / g * w * ϕₜₜ)dom[:dΓκ]
end

function damping(pf::PotentialFlow, fs::FreeSurface, dom::IntegrationDomains, x_t, y)
    ϕ_sym = variable_symbol(pf)
    κ_sym = variable_symbol(fs)
    ϕₜ = x_t[ϕ_sym]
    κₜ = x_t[κ_sym]
    w  = y[ϕ_sym]
    u  = y[κ_sym]
    βₕ = fs.βₕ
    if _simulation_domain(dom) == :time && haskey(dom, :αₕ)
        αₕ = stabilization_parameter(fs, dom)
        return ∫(βₕ * (u + αₕ * w) * ϕₜ - w * κₜ)dom[:dΓκ]
    end
    ∫(βₕ * u * ϕₜ - βₕ * w * κₜ)dom[:dΓκ]
end

function stiffness(pf::PotentialFlow, fs::FreeSurface, dom::IntegrationDomains, x, y)
    ϕ_sym = variable_symbol(pf)
    κ_sym = variable_symbol(fs)
    ϕ = x[ϕ_sym]
    κ = x[κ_sym]
    w = y[ϕ_sym]
    u = y[κ_sym]

    val = nothing
    if _simulation_domain(dom) == :time && haskey(dom, :αₕ)
        αₕ = stabilization_parameter(fs, dom)
        val = ∫(fs.βₕ * fs.g * αₕ * w * κ)dom[:dΓκ]
    end

    for bc in _active_damping_zone_bcs(pf)
        dΓ = dom[bc.domain]
        nΓ = _boundary_normal(dom, bc.domain)
        μ₁ = _as_space_function(bc.μ₁)
        μ₂ = _as_space_function(bc.μ₂)
        ∇ₙϕ = ∇(ϕ) ⋅ nΓ
        zone_val = ∫(μ₁ * ∇ₙϕ * u - (μ₂ * κ * w))dΓ
        val = _add_contribution(val, zone_val)
    end

    isnothing(val) && return ∫(0.0 * w * κ)dom[:dΓκ]
    return val
end

function rhs(pf::PotentialFlow, fs::FreeSurface, dom::IntegrationDomains, f, y)
    κ_sym = variable_symbol(fs)
    u = y[κ_sym]

    val = nothing
    for bc in _active_damping_zone_bcs(pf)
        dΓ = dom[bc.domain]
        ∇ₙϕd = _normal_damped(bc, dom)
        zone_val = ∫(∇ₙϕd * u)dΓ
        val = _add_contribution(val, zone_val)
    end

    isnothing(val) && return ∫(0.0 * u)dom[:dΓκ]
    return val
end


# =========================================================================
# Resonator ↔ Structure coupling
#
# Cross-terms between q_i DOFs and the structural displacement η.
# Uses variable_symbol(s) to identify which structure field to couple to.
# =========================================================================

# Resonator-structure coupling has damping and stiffness terms only.
has_damping_form(::Vector{ResonatorSingle}, ::Structure) = true
has_stiffness_form(::Vector{ResonatorSingle}, ::Structure) = true

function damping(resn::Vector{ResonatorSingle}, s::Structure,
                 dom::IntegrationDomains, x_t, y)
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

function stiffness(resn::Vector{ResonatorSingle}, s::Structure,
                   dom::IntegrationDomains, x, y)
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
