# ─────────────────────────────────────────────────────────────
# IntegrationDomains: Dict-based container for Gridap measures,
# normals, and metadata used by weak form assembly.
# ─────────────────────────────────────────────────────────────

"""
    IntegrationDomains(; key=value, ...)
    IntegrationDomains(dict::Dict{Symbol,Any})

Dict-based container for Gridap measures, normals, DiracDeltas,
and any other domain data needed by weak form methods.

Each `weakform` dispatch accesses only the keys it needs via
`dom[:key]`.  No fixed schema — new keys can be added without
changing this type.

# Standard key conventions (not enforced)
- `:dΩ`     — fluid interior measure
- `:dΓ_fs`  — free-surface measure (outside structure)
- `:dΓ_s`   — structure surface measure
- `:dΛ_s`, `:n_Λ_s`, `:h_s` — beam skeleton measures/normals + mesh size
- `:dΛ_sb`, `:n_Λ_sb`  — structure boundary (fixed BC Neumann)
- `:dΓ_in`, `:dΓ_out`  — inlet / outlet radiation boundaries
- `:dΓ_d_1`, `:dΓ_d_2` — damping zone measures
- `:δ_p`    — vector of DiracDelta functionals (resonator points)
"""
struct IntegrationDomains
    data::Dict{Symbol, Any}
end

IntegrationDomains(; kwargs...) =
    IntegrationDomains(Dict{Symbol, Any}(k => v for (k, v) in pairs(kwargs)))

Base.getindex(d::IntegrationDomains, k::Symbol)            = d.data[k]
Base.haskey(d::IntegrationDomains, k::Symbol)               = haskey(d.data, k)
Base.get(d::IntegrationDomains, k::Symbol, default)         = get(d.data, k, default)
Base.setindex!(d::IntegrationDomains, val, k::Symbol)       = (d.data[k] = val)
Base.keys(d::IntegrationDomains)                            = keys(d.data)
