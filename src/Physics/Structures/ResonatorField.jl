"""
    ResonatorField <: PhysicsParameters

Homogenised (continuous) field of mass-spring-damper resonators attached to a
structure, the continuum limit of a dense [`ResonatorArray`](@ref).

The resonator displacement `q` is a field on the same structural surface as
the host displacement `η` (the structure whose `space_domain_symbol` equals
`space_domain_symbol` of the field), so it works for line structures in 2D
(beams, membranes) and surface structures in 3D (plates, membranes) alike.
`n(x)` is the number of resonators per unit structure measure: per unit
length on a 1D structure, per unit area on a 2D one.

Weak forms, normalised by the fluid density `ρw` like the structural forms
(test functions `ξ` of `q` and `v` of `η`):

```math
\\begin{aligned}
m(\\ddot q, \\xi) &= \\int_{\\Gamma} \\frac{n M}{\\rho_w} \\, \\ddot q \\, \\xi \\, \\mathrm{d}\\Gamma \\\\
c(\\dot q, \\dot\\eta; \\xi, v) &= \\int_{\\Gamma} \\frac{n C}{\\rho_w} (\\dot q - \\dot\\eta)(\\xi - v) \\, \\mathrm{d}\\Gamma \\\\
k(q, \\eta; \\xi, v) &= \\int_{\\Gamma} \\frac{n K}{\\rho_w} (q - \\eta)(\\xi - v) \\, \\mathrm{d}\\Gamma
\\end{aligned}
```

The `q`-only parts are the single-entity `mass`/`damping`/`stiffness`; the
cross and `η`-`η` parts are the `(ResonatorField, Structure)` coupling forms.
The assembled operator is symmetric. In the frequency domain (time factor
``e^{-i\\omega t}``) eliminating `q` gives the host an effective mass density

```math
m_{\\mathrm{eff},\\varrho} = m_\\varrho + \\frac{n M}{\\rho_w}
\\frac{K - i\\omega C}{K - i\\omega C - \\omega^2 M},
```

with ``\\omega_r^2 = K/M`` and ``\\zeta = C/(2\\sqrt{KM})``, as for
[`ResonatorSingle`](@ref). A `ResonatorArray` with spacing `a` (area `a²`
in 3D) corresponds to `n = 1/a` (`1/a²`). The field enters the same
`physics` vector as the host structure; it is time- and frequency-domain
ready (no extra forcing form).

# Fields
- `M`, `K`, `C` — mass [kg], stiffness [N/m] and viscous damping [N·s/m] of
  one resonator; each a number or a function `x -> value` (graded
  populations). `C` defaults to 0.
- `n` — number density [1/m or 1/m²]: a number or a function `x -> n(x)`;
  default 1.
- `ρw::Float64` — fluid density [kg/m³]; default 1025.
- `symbol::Symbol` — field unknown; default `:q`.
- `space_domain_symbol::Symbol` — triangulation of the field, and the
  host-structure selector; default `:Γη`.
- `fe::FESpaceConfig` — FE space of `q`. Only `L²` products of `q` appear, so
  a continuous Lagrangian space of the host's order is a good choice (no
  Dirichlet tags needed).

# Example
```julia
beam = EulerBernoulliBeam(L=L, mᵨ=m, EIᵨ=EI, fe=fe)
res  = ResonatorField(M=10.0, K=10.0*ωr^2, C=2ζ*10.0*ωr, n=0.5, fe=fe)
physics = PhysicsParameters[PotentialFlow(...), beam, res]
```

See also: [`ResonatorArray`](@ref), [`effective_mass_density`](@ref)
"""
@with_kw struct ResonatorField <: PhysicsParameters
    M::Union{Real,Function}
    K::Union{Real,Function}
    C::Union{Real,Function} = 0.0
    n::Union{Real,Function} = 1.0
    ρw::Float64 = 1025.0
    symbol::Symbol = :q
    space_domain_symbol::Symbol = :Γη
    fe::FESpaceConfig = FESpaceConfig()
end

variable_symbol(rf::ResonatorField) = rf.symbol

function print_parameters(rf::ResonatorField)
    _str(p) = p isa Function ? "x -> …" : @sprintf("%.4f", p)
    @printf("\n[MSG] Resonator Field Properties:\n")
    @printf("[VAL] M = %s kg\n", _str(rf.M))
    @printf("[VAL] K = %s N/m\n", _str(rf.K))
    @printf("[VAL] C = %s Ns/m\n", _str(rf.C))
    @printf("[VAL] n = %s 1/m^d\n", _str(rf.n))
    @printf("[VAL] space_domain_symbol = %s\n", rf.space_domain_symbol)
    println()
end

"""
    effective_mass_density(rf::ResonatorField, ω, x=nothing)

Added dynamic mass per unit structure measure, divided by `ρw`, that the
field contributes to its host at angular frequency `ω` (time factor
``e^{-i\\omega t}``): ``n M (K - i\\omega C) / (\\rho_w (K - i\\omega C - \\omega^2 M))``.
`x` is needed only when a parameter is a function of position.
"""
function effective_mass_density(rf::ResonatorField, ω, x=nothing)
    M, K, C, n = (_rf_eval(p, x) for p in (rf.M, rf.K, rf.C, rf.n))
    Kc = K - im * ω * C
    return n * M / rf.ρw * Kc / (Kc - ω^2 * M)
end

_rf_eval(p::Real, x) = p
_rf_eval(p::Function, x) = p(x)

_rf_iszero(p) = p isa Real && iszero(p)

# n·p/ρw as a number, or as a CellField on `trian` if any factor depends on x.
function _rf_coefficient(rf::ResonatorField, p, trian)
    if rf.n isa Real && p isa Real
        return Float64(rf.n * p / rf.ρw)
    end
    return CellField(x -> _rf_eval(rf.n, x) * _rf_eval(p, x) / rf.ρw, trian)
end

function _rf_self_form(rf::ResonatorField, p, dom::IntegrationDomains, x, y)
    sym = variable_symbol(rf)
    q, ξ = x[sym], y[sym]
    c = _rf_coefficient(rf, p, get_triangulation(ξ))
    ∫(c * q * ξ)_space_measure(dom, rf)
end

"""
    mass(rf::ResonatorField, dom, x_tt, y)

Resonator inertia `∫ n M/ρw q̈ ξ dΓ`.
"""
mass(rf::ResonatorField, dom::IntegrationDomains, x_tt, y) =
    _rf_self_form(rf, rf.M, dom, x_tt, y)

"""
    damping(rf::ResonatorField, dom, x_t, y)

Dashpot, resonator-only part `∫ n C/ρw q̇ ξ dΓ`.
"""
damping(rf::ResonatorField, dom::IntegrationDomains, x_t, y) =
    _rf_self_form(rf, rf.C, dom, x_t, y)

"""
    stiffness(rf::ResonatorField, dom, x, y)

Spring, resonator-only part `∫ n K/ρw q ξ dΓ`.
"""
stiffness(rf::ResonatorField, dom::IntegrationDomains, x, y) =
    _rf_self_form(rf, rf.K, dom, x, y)

has_mass_form(::ResonatorField) = true
has_damping_form(rf::ResonatorField) = !_rf_iszero(rf.C)
has_stiffness_form(::ResonatorField) = true
has_rhs_form(::ResonatorField) = false
