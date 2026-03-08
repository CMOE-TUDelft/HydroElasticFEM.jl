# ─────────────────────────────────────────────────────────────
# WeakFormDomains: Dict-based container for Gridap measures,
# normals, and metadata used by weak form assembly.
# ─────────────────────────────────────────────────────────────

"""
    WeakFormDomains(; key=value, ...)

Dict-based container for Gridap measures, normals, DiracDeltas,
and any other domain data needed by weak form methods.

Each `weakform` dispatch accesses only the keys it needs via
`dom[:key]`.  No fixed schema — new keys can be added without
changing this type.

# Standard key conventions (not enforced)
- `:dΩ`     — fluid interior measure
- `:dΓ_fs`  — free-surface measure (outside structure)
- `:dΓ_s`   — structure surface measure
- `:dΛ_s`, `:n_Λ_s`, `:h_s`, `:γ_s` — beam skeleton + Nitsche params
- `:dΛ_sb`, `:n_Λ_sb`  — structure boundary (fixed BC Neumann)
- `:dΓ_in`, `:dΓ_ot`   — inlet / outlet radiation boundaries
- `:dΓ_d1`, `:dΓ_d2`   — damping zone measures
- `:δ_p`    — vector of DiracDelta functionals (resonator points)
"""
struct WeakFormDomains
    data::Dict{Symbol, Any}
end

WeakFormDomains(; kwargs...) =
    WeakFormDomains(Dict{Symbol, Any}(k => v for (k, v) in pairs(kwargs)))

Base.getindex(d::WeakFormDomains, k::Symbol)            = d.data[k]
Base.haskey(d::WeakFormDomains, k::Symbol)               = haskey(d.data, k)
Base.get(d::WeakFormDomains, k::Symbol, default)         = get(d.data, k, default)
Base.setindex!(d::WeakFormDomains, val, k::Symbol)       = (d.data[k] = val)
Base.keys(d::WeakFormDomains)                            = keys(d.data)
