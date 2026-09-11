"""
    ResonatorSingle <: PhysicsParameters

Parameters for a single lumped-parameter mass-spring-damper resonator,
used to model locally resonant meta-structures.

The resonator couples to the surrounding fluid via a delta-function
(`δ_p`) point interaction at `XZ`.

# Fields
- `M::Float64`   — Resonator mass [kg]
- `K::Float64`   — Spring stiffness [N/m]
- `C::Float64`   — Viscous damping coefficient [N·s/m]; default 0.0
- `ρw::Float64`  — Fluid density [kg/m³]; default 1025.0
- `XZ::VectorValue{2,Float64}` — Resonator position `(x, z)` [m]
- `symbol::Symbol` — Field unknown symbol; default `:q`
- `space_domain_symbol::Symbol` — Triangulation key; default `:Ω`
- `fe::FESpaceConfig` — FE space parameters
- `ωn1::Float64` — Undamped natural frequency [rad/s], derived as `√(K/M)`

See also: [`resonator_array`](@ref)
"""
@with_kw struct ResonatorSingle <: PhysicsParameters
    M::Float64
    K::Float64
    C::Float64     = 0.0
    ρw::Float64    = 1025.0
    XZ::VectorValue{2,Float64} = VectorValue(0.0, 0.0)
    ωn1::Float64   = sqrt(K / M)
    symbol::Symbol = :q
    space_domain_symbol::Symbol = :Ω
    fe::FESpaceConfig = FESpaceConfig()
end

variable_symbol(s::ResonatorSingle) = s.symbol

"""
    ResonatorArray <: PhysicsParameters

Collects a set of [`ResonatorSingle`](@ref) entities into a single physics
entity, so a resonator array can be passed to the assembler like any other
`PhysicsParameters` subtype (in place of a bare `Vector{ResonatorSingle}`).

Each resonator contributes an independent constant-space DOF `q_i` (indexed
by position in `resonators`). Supports `length`, `isempty`, indexing, and
iteration, delegating to the wrapped `resonators` vector.

# Fields
- `resonators::Vector{ResonatorSingle}` — the individual resonators, in DOF order

See also: [`resonator_array`](@ref)
"""
struct ResonatorArray <: PhysicsParameters
    resonators::Vector{ResonatorSingle}
end

Base.length(ra::ResonatorArray) = length(ra.resonators)
Base.isempty(ra::ResonatorArray) = isempty(ra.resonators)
Base.getindex(ra::ResonatorArray, i) = ra.resonators[i]
Base.iterate(ra::ResonatorArray, state...) = iterate(ra.resonators, state...)

function print_parameters(resn::ResonatorSingle)
    @printf("\n[MSG] Resonator Properties:\n")
    @printf("[VAL] M = %.4f kg\n", resn.M)
    @printf("[VAL] K = %.4f N/m\n", resn.K)
    @printf("[VAL] C = %.4f Ns/m\n", resn.C)
    @printf("[VAL] XZ = (%.4f, %.4f) m\n", resn.XZ[1], resn.XZ[2])
    @printf("[VAL] ωn1 = %.4f rad/s\n", resn.ωn1)
    println()
end

"""
    print_parameters(ra::ResonatorArray)

Print parameters for every resonator in the array by delegating to the
single-resonator `print_parameters` method.
"""
function print_parameters(ra::ResonatorArray)
    print_parameters.(ra.resonators)
end

"""
    resonator_array(N, M, K, C, XZ; ρw=1025.0, fe=FESpaceConfig()) -> ResonatorArray

Create `N` identical resonators positioned at the locations in `XZ`.

# Arguments
- `N::Int`  — Number of resonators (must be positive)
- `M::Real` — Mass [kg]
- `K::Real` — Spring stiffness [N/m]
- `C::Real` — Viscous damping [N·s/m]
- `XZ::Vector{VectorValue{2,Float64}}` — Position list of length `N` [m]
- `ρw::Real` — Fluid density [kg/m³]; default 1025.0
- `fe::FESpaceConfig` — FE space parameters; default `FESpaceConfig()`
"""
function resonator_array(N::Int, M::Real, K::Real, C::Real,
                         XZ::Vector{VectorValue{2,Float64}};
                         ρw::Real=1025.0,
                         fe::FESpaceConfig=FESpaceConfig())
    N > 0 || throw(ArgumentError("N must be positive (got $N)"))
    length(XZ) == N || throw(ArgumentError("XZ must be of length N"))
    ResonatorArray([ResonatorSingle(M=M, K=K, C=C, ρw=ρw, XZ=xz, fe=fe) for xz in XZ])
end

"""
    resonator_array(N, M, K, C, XZ; ρw=1025.0, fe=FESpaceConfig()) -> ResonatorArray

Create `N` resonators with individually specified parameters.

# Arguments
- `N::Int`            — Number of resonators (must be positive)
- `M::Vector{<:Real}` — Masses [kg], length `N`
- `K::Vector{<:Real}` — Spring stiffnesses [N/m], length `N`
- `C::Vector{<:Real}` — Viscous damping coefficients [N·s/m], length `N`
- `XZ::Vector{VectorValue{2,Float64}}` — Positions [m], length `N`
- `ρw::Real`          — Fluid density [kg/m³]; default 1025.0
- `fe::FESpaceConfig` — FE space parameters; default `FESpaceConfig()`
"""
function resonator_array(N::Int, M::Vector{<:Real}, K::Vector{<:Real},
                         C::Vector{<:Real},
                         XZ::Vector{VectorValue{2,Float64}};
                         ρw::Real=1025.0, 
                         fe::FESpaceConfig=FESpaceConfig())
    N > 0 || throw(ArgumentError("N must be positive (got $N)"))
    (length(M) == N && length(K) == N &&
     length(C) == N && length(XZ) == N) ||
        throw(ArgumentError("M, K, C, and XZ must be of length N"))
    ResonatorArray([ResonatorSingle(M=m, K=k, C=c, ρw=ρw, XZ=xz, fe=fe) for (m, k, c, xz) in zip(M, K, C, XZ)])
end

# ── Single-variable weak forms: mass, damping, stiffness, rhs ──
#    Only q_i terms — no coupling to structure η

"""
    mass(ra::ResonatorArray, dom::IntegrationDomains, x_tt, y)

Resonator array inertia (mass) bilinear form.

Assembles the sum of point-mass contributions over all resonators at their
attachment locations using the delta-function distributions `δ_p`:

```math
\\sum_i M_i \\, \\delta_i(\\partial_{tt} q_i \\cdot \\xi_i)
```

# Arguments
- `ra::ResonatorArray`: array of resonator entities
- `dom::IntegrationDomains`: integration measures (requires `:delta_p`, `:dΩ`)
- `x_tt`: second time-derivative trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`
"""
function mass(ra::ResonatorArray, dom::IntegrationDomains, x_tt, y)
    resn = ra.resonators
    δ_p = dom[:δ_p]
    ξ1  = y[variable_symbol(resn[1])]
    q1  = x_tt[variable_symbol(resn[1])]
    dΩ  = _space_measure(dom, resn)
    val = ∫((ξ1 ⋅ q1) * 0.0)dΩ
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qₜₜi = x_tt[variable_symbol(ri)]
        ξi   = y[variable_symbol(ri)]
        val += ri.M * δi(qₜₜi ⋅ ξi)
    end
    return val
end

"""
    damping(ra::ResonatorArray, dom::IntegrationDomains, x_t, y)

Resonator array self-damping bilinear form.

Assembles the sum of viscous damping contributions over all resonators:

```math
\\sum_i C_i \\, \\delta_i(\\partial_t q_i \\cdot \\xi_i)
```

# Arguments
- `ra::ResonatorArray`: array of resonator entities
- `dom::IntegrationDomains`: integration measures (requires `:delta_p`, `:dΩ`)
- `x_t`: first time-derivative trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`
"""
function damping(ra::ResonatorArray, dom::IntegrationDomains, x_t, y)
    resn = ra.resonators
    δ_p = dom[:δ_p]
    ξ1  = y[variable_symbol(resn[1])]
    q1  = x_t[variable_symbol(resn[1])]
    dΩ  = _space_measure(dom, resn)
    val = ∫((ξ1 ⋅ q1) * 0.0)dΩ
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qₜi = x_t[variable_symbol(ri)]
        ξi  = y[variable_symbol(ri)]
        val += ri.C * δi(qₜi ⋅ ξi)
    end
    return val
end

"""
    stiffness(ra::ResonatorArray, dom::IntegrationDomains, x, y)

Resonator array self-stiffness bilinear form.

Assembles the sum of spring contributions over all resonators:

```math
\\sum_i K_i \\, \\delta_i(q_i \\cdot \\xi_i)
```

# Arguments
- `ra::ResonatorArray`: array of resonator entities
- `dom::IntegrationDomains`: integration measures (requires `:delta_p`, `:dΩ`)
- `x`: trial `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`
"""
function stiffness(ra::ResonatorArray, dom::IntegrationDomains, x, y)
    resn = ra.resonators
    δ_p = dom[:δ_p]
    ξ1  = y[variable_symbol(resn[1])]
    q1  = x[variable_symbol(resn[1])]
    dΩ  = _space_measure(dom, resn)
    val = ∫((ξ1 ⋅ q1) * 0.0)dΩ
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        qi = x[variable_symbol(ri)]
        ξi = y[variable_symbol(ri)]
        val += ri.K * δi(qi ⋅ ξi)
    end
    return val
end

"""
    rhs(ra::ResonatorArray, dom::IntegrationDomains, f, y)

Resonator array right-hand side (external forcing) linear form.

Assembles point loads on each resonator DOF using the delta distributions:

```math
\\sum_i \\delta_i(f_i \\cdot \\xi_i)
```

Returns `false` for the form-presence trait (`has_rhs_form` returns `false`);
this method is only called when an external resonator forcing is explicitly
provided.

# Arguments
- `ra::ResonatorArray`: array of resonator entities
- `dom::IntegrationDomains`: integration measures (requires `:delta_p`, `:dΩ`)
- `f`: forcing `FieldMap`
- `y`: test `FieldMap`

# Returns
- `Gridap.FESpaces.DomainContribution`
"""
function rhs(ra::ResonatorArray, dom::IntegrationDomains, f, y)
    resn = ra.resonators
    δ_p = dom[:δ_p]
    ξ1  = y[variable_symbol(resn[1])]
    q1  = f[variable_symbol(resn[1])]
    dΩ  = _space_measure(dom, resn)
    val = ∫((ξ1 ⋅ q1) * 0.0)dΩ
    for (i, (δi, ri)) in enumerate(zip(δ_p, resn))
        fi = f[variable_symbol(ri)]
        ξi = y[variable_symbol(ri)]
        val += δi(fi ⋅ ξi)
    end
    return val
end

# Form-presence traits for ResonatorArray
# (dedicated methods since defaults for PhysicsParameters assume a single field)
has_mass_form(ra::ResonatorArray) = !isempty(ra.resonators)
has_damping_form(ra::ResonatorArray) = !isempty(ra.resonators)
has_stiffness_form(ra::ResonatorArray) = !isempty(ra.resonators)
has_rhs_form(::ResonatorArray) = false
