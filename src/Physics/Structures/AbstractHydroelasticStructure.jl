# =============================================================================
# AbstractHydroelasticStructure.jl
#
# `Membrane` and `EulerBernoulliBeam` (and, by construction, the new
# `TensionedEulerBernoulliBeam`) share an identical pattern across their
# `mass`, `damping`, and `stiffness` weak forms:
#
#   mass(s, dom, x_tt, y)      = ∫ mᵨ v ηₜₜ dΩ                       (identical)
#   stiffness(s, dom, x, y)    = ∫ g v η dΩ  +  <elastic operator>    (g-term identical)
#   damping(s, dom, x_t, y)    = τ · <elastic operator>(x_t, y)       (stiffness-proportional
#                                                                       Rayleigh damping)
#   rhs(s, dom, f, y)          = ∫ v f[sym] dΩ                        (identical)
#
# where `<elastic operator>` is the one piece that actually differs between
# structures — membrane tension (`Tρ·∇v⋅∇η`) vs. Euler-Bernoulli bending
# (`EIρ·Δv·Δη` plus its C/DG skeleton terms). `AbstractHydroelasticStructure`
# factors the four identical forms above into shared default methods, so a
# new structure only has to implement `stiffness_operator` (and, optionally,
# `extra_stiffness_form` for contributions — such as joint springs — that are
# *not* subject to Rayleigh damping).
#
# `KirchhoffLovePlate`, `TimoshenkoBeam`, and `ResonatorArray` do not follow
# this pattern closely enough to benefit (no damping form, multi-field, or a
# fundamentally different mass operator) and are left as direct `Structure`
# subtypes.
# =============================================================================

"""
    AbstractHydroelasticStructure <: Structure

Shared base type for single-field hydroelastic structures whose `mass`,
`damping`, and `stiffness` weak forms follow the common pattern

```math
\\begin{aligned}
m(\\ddot\\eta, v) &= \\int_{\\Gamma_\\eta} m_\\varrho \\, v \\, \\ddot\\eta \\, \\mathrm{d}\\Gamma \\\\
c(\\dot\\eta, v)   &= \\tau \\cdot k_{\\mathrm{op}}(\\dot\\eta, v) \\\\
k(\\eta, v)        &= \\int_{\\Gamma_\\eta} g \\, v \\, \\eta \\, \\mathrm{d}\\Gamma + k_{\\mathrm{op}}(\\eta, v) + k_{\\mathrm{extra}}(\\eta, v) \\\\
l(v)               &= \\int_{\\Gamma_\\eta} v \\, f \\, \\mathrm{d}\\Gamma
\\end{aligned}
```

where ``k_{\\mathrm{op}}`` is the structure-specific elastic stiffness
operator (membrane tension, Euler-Bernoulli bending, or a combination) and
``k_{\\mathrm{extra}}`` collects any additional stiffness contribution that
should *not* be scaled by the stiffness-proportional Rayleigh damping
coefficient `τ` (e.g. rotational-spring joints).

Concrete subtypes get `mass`, `damping`, `stiffness`, and `rhs` for free by
implementing only:

- [`mass_density(s)`](@ref) — returns `s`'s `mᵨ` field.
- [`damping_parameter(s)`](@ref) — returns `s`'s `τ` field.
- [`stiffness_operator(s, dom, x, y)`](@ref) — the elastic operator, reused
  (scaled by `τ`) for the damping form.

and, only if needed:

- [`gravitational_acceleration(s)`](@ref) — defaults to `s.g`.
- [`extra_stiffness_form(s, dom, x, y)`](@ref) — defaults to `nothing`
  (no contribution).

# See also
[`Membrane`](@ref), [`EulerBernoulliBeam`](@ref),
[`TensionedEulerBernoulliBeam`](@ref)
"""
abstract type AbstractHydroelasticStructure <: Structure end

"""
    mass_density(s::AbstractHydroelasticStructure) -> Float64

Mass per unit manifold measure, normalised by fluid density `ρw` (the
structure's `mᵨ` field). Must be implemented by every concrete subtype.
"""
function mass_density(s::AbstractHydroelasticStructure)
    error("mass_density not implemented for $(typeof(s))")
end

"""
    damping_parameter(s::AbstractHydroelasticStructure) -> Float64

Stiffness-proportional (Kelvin-Voigt/Rayleigh) structural damping
coefficient `τ`. Must be implemented by every concrete subtype; return `0.0`
to disable structural damping (the resulting `damping` contribution is then
identically zero, matching `has_damping_form(s) = true` with no effect).
"""
function damping_parameter(s::AbstractHydroelasticStructure)
    error("damping_parameter not implemented for $(typeof(s))")
end

"""
    gravitational_acceleration(s::AbstractHydroelasticStructure) -> Float64

Gravitational acceleration used in the hydrostatic restoring term
`∫ g v η dΓ`. Defaults to `s.g`; override if a subtype stores it under a
different field name.
"""
gravitational_acceleration(s::AbstractHydroelasticStructure) = s.g

"""
    stiffness_operator(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x, y)

Elastic stiffness contribution specific to structure `s` (membrane tension,
Euler-Bernoulli bending, or a combination), *excluding* the hydrostatic
restoring term and any [`extra_stiffness_form`](@ref) contribution.

This is the single mandatory specialisation point for a new
`AbstractHydroelasticStructure`: [`stiffness`](@ref) adds the shared
hydrostatic term and `extra_stiffness_form` around it, and [`damping`](@ref)
reuses it — scaled by [`damping_parameter`](@ref) — to assemble
stiffness-proportional Rayleigh damping without any additional code.

Must be implemented by every concrete subtype.

# Returns
`Gridap.FESpaces.DomainContribution`
"""
function stiffness_operator(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x, y)
    error("stiffness_operator not implemented for $(typeof(s))")
end

"""
    extra_stiffness_form(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x, y)

Optional additional stiffness contribution that is *not* subject to
stiffness-proportional Rayleigh damping (e.g. rotational-spring joints on a
beam, which are purely elastic connections). Defaults to `nothing`, i.e. no
extra contribution; [`stiffness`](@ref) combines it with the damped part via
`_add_contribution`, which treats `nothing` as an additive identity.
"""
extra_stiffness_form(::AbstractHydroelasticStructure, dom::IntegrationDomains, x, y) = nothing

"""
    mass(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x_tt, y)

Shared inertia (mass) bilinear form for every `AbstractHydroelasticStructure`:

```math
\\int_{\\Gamma_\\eta} m_\\varrho \\, v \\, \\partial_{tt}\\eta \\, \\mathrm{d}\\Gamma_\\eta
```

Uses [`mass_density(s)`](@ref); concrete subtypes do not need to implement
`mass` themselves.
"""
function mass(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x_tt, y)
    sym = variable_symbol(s)
    ηₜₜ = x_tt[sym]
    v   = y[sym]
    dΩ  = _space_measure(dom, s)
    ∫(mass_density(s) * v * ηₜₜ)dΩ
end

"""
    damping(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x_t, y)

Shared stiffness-proportional (Kelvin-Voigt/Rayleigh) structural damping
bilinear form, obtained by evaluating [`stiffness_operator`](@ref) at the
velocity trial field `x_t` and scaling it by [`damping_parameter(s)`](@ref):

```math
\\tau \\cdot k_{\\mathrm{op}}(\\partial_t \\eta, v)
```

This is exact — not approximate — because `stiffness_operator` is linear in
its material coefficients (a scalar factor commutes with `mean`, `jump`, and
`∫`), so scaling the whole operator by `τ` is algebraically identical to
scaling each material coefficient inside it by `τ` before assembling.
Concrete subtypes do not need to implement `damping` themselves.
"""
function damping(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x_t, y)
    damping_parameter(s) * stiffness_operator(s, dom, x_t, y)
end

"""
    stiffness(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x, y)

Shared stiffness bilinear form: hydrostatic restoring term plus the
structure-specific [`stiffness_operator`](@ref) plus any
[`extra_stiffness_form`](@ref):

```math
\\int_{\\Gamma_\\eta} g \\, v \\, \\eta \\, \\mathrm{d}\\Gamma_\\eta
+ k_{\\mathrm{op}}(\\eta, v) + k_{\\mathrm{extra}}(\\eta, v)
```

Concrete subtypes do not need to implement `stiffness` themselves.
"""
function stiffness(s::AbstractHydroelasticStructure, dom::IntegrationDomains, x, y)
    sym = variable_symbol(s)
    η = x[sym]
    v = y[sym]
    dΩ = _space_measure(dom, s)
    gravity_term = ∫(gravitational_acceleration(s) * v * η)dΩ
    total = gravity_term + stiffness_operator(s, dom, x, y)
    return _add_contribution(total, extra_stiffness_form(s, dom, x, y))
end

"""
    rhs(s::AbstractHydroelasticStructure, dom::IntegrationDomains, f, y)

Shared right-hand side (applied load) linear form:

```math
\\int_{\\Gamma_\\eta} v \\, f_\\eta \\, \\mathrm{d}\\Gamma_\\eta
```

Concrete subtypes do not need to implement `rhs` themselves.
"""
function rhs(s::AbstractHydroelasticStructure, dom::IntegrationDomains, f, y)
    sym = variable_symbol(s)
    v = y[sym]
    dΩ = _space_measure(dom, s)
    ∫(v * f[sym])dΩ
end
