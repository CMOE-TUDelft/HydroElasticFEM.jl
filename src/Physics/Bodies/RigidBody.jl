# ==========================================================================
# RigidBody — floating rigid body with global (constant) DOFs
# ==========================================================================
#
# A rigid body owns one vector unknown ξ (one component per active mode),
# held in a `ConstantFESpace` on Ω like a resonator DOF.  With the reference
# point x_r, translations ξ_t and small rotations α, a body point moves by
#
#     u(x) = ξ_t + α × (x − x_r).
#
# In 2D, x = (x, z) and the only rotation is about y, so
#     u = (ξ₁ + α r_z, ξ₃ − α r_x),   r = x − x_r,
# which is the 3D motion restricted to the (x, z) plane.
#
# On the wetted surface Γ_b, with n the outward normal of the fluid (into the
# body), the normal velocity is u̇·n = ξ̇·N(x), where N is the generalized
# normal
#     3D:  N = (n, r × n)
#     2D:  N = (n_x, n_z, r_z n_x − r_x n_z)
# restricted to the active modes.
#
# Weak forms (normalized by ρw, like every other entity), test function ζ:
#     mass       ζ·(M̂ ξ̈)                         M̂ = M_rb / ρw
#     damping    ζ·(B̂ ξ̇)                         B̂ = B_ext / ρw
#     stiffness  ζ·(K̂ ξ) [+ ∫_Γb g (ζ·N)(ξ·N)]    K̂ = (C_hs + K_ext) / ρw
# The fluid coupling lives in `CouplingTerms.jl`.  Constant matrices are
# assembled as ∫_Γb ζ·(A ξ) / |Γ_b| dΓ, which needs no reference-point
# geometry (x_r may lie inside the body, outside the fluid mesh).

"""
    AbstractRigidBody <: PhysicsParameters

Base type for rigid bodies: entities whose unknowns are a few global DOFs
(rigid-body modes) rather than a field.  Not a `Structure`: the structure
couplings assume a scalar deflection field.
"""
abstract type AbstractRigidBody <: PhysicsParameters end

const RIGID_BODY_MODES_2D = (:surge, :heave, :pitch)
const RIGID_BODY_MODES_3D = (:surge, :sway, :heave, :roll, :pitch, :yaw)

"""
    rigid_body_modes(ambient_dim) -> NTuple{N,Symbol}

All rigid-body modes in `ambient_dim` dimensions, in DOF order:
`(:surge, :heave, :pitch)` in 2D, `(:surge, :sway, :heave, :roll, :pitch, :yaw)`
in 3D.
"""
rigid_body_modes(ambient_dim::Int) =
    ambient_dim == 2 ? RIGID_BODY_MODES_2D :
    ambient_dim == 3 ? RIGID_BODY_MODES_3D :
    throw(ArgumentError("rigid bodies need ambient_dim 2 or 3, got $ambient_dim"))

"""
    RigidBody(; mass, inertia, x_G, kwargs...)

Floating rigid body with the DOFs `modes` (a subset of
[`rigid_body_modes`](@ref)), in 2D `(x, z)` or 3D `(x, y, z)`.

Dimensional inputs (kg, N/m, …) are given per unit width in 2D; the weak
forms divide them by `ρw`, like the other entities.

# Keywords
- `mass::Real` — body mass `m`.
- `inertia` — moment of inertia about the centre of gravity: a scalar `I_yy`
  in 2D, a 3×3 matrix in 3D.
- `x_G` — centre of gravity; its length sets `ambient_dim`.
- `x_ref = x_G` — reference point of the rotational modes.
- `modes` — active modes (default: all for `ambient_dim`).
- `hydrostatics = :zero_draft` — `:zero_draft` (restoring `∫ g (ζ·N)(ξ·N)` on
  the wetted surface, exact for a body of zero draft and identical to the
  hydrostatic term of the elastic structures), `:none`, or a full-mode
  restoring matrix (dimensional).
- `K_ext`, `B_ext` — external stiffness (moorings) and damping (PTO), full-mode
  matrices about `x_ref` (dimensional; default zero).
- `ρw = 1025.0`, `g = 9.81`.
- `symbol = :ξ` — unknown symbol.
- `space_domain_symbol = :Γ_b` — triangulation key of the wetted surface.
- `fe = FESpaceConfig()` — only `vector_type` is used.

A mode with no restoring and no fluid coupling (surge, sway, yaw of a
zero-draft body) makes the system singular; it must be restrained through
`K_ext` or left out of `modes`.

# Example
```julia
body = RigidBody(mass = 2.0e4, inertia = 6.7e5, x_G = [20.0, 0.0],
                 modes = (:heave, :pitch), symbol = :ξ, space_domain_symbol = :Γ_b)
```
"""
struct RigidBody <: AbstractRigidBody
    mass::Float64
    inertia::Matrix{Float64}            # about x_G; 1×1 (I_yy) in 2D, 3×3 in 3D
    x_G::Vector{Float64}
    x_ref::Vector{Float64}
    ambient_dim::Int
    modes::Tuple{Vararg{Symbol}}
    hydrostatics::Union{Symbol,Matrix{Float64}}
    K_ext::Matrix{Float64}
    B_ext::Matrix{Float64}
    ρw::Float64
    g::Float64
    symbol::Symbol
    space_domain_symbol::Symbol
    fe::FESpaceConfig
end

function RigidBody(; mass::Real, inertia, x_G::AbstractVector{<:Real},
                   x_ref::AbstractVector{<:Real} = x_G,
                   modes = rigid_body_modes(length(x_G)),
                   hydrostatics::Union{Symbol,AbstractMatrix{<:Real}} = :zero_draft,
                   K_ext::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
                   B_ext::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
                   ρw::Real = 1025.0, g::Real = 9.81,
                   symbol::Symbol = :ξ, space_domain_symbol::Symbol = :Γ_b,
                   fe::FESpaceConfig = FESpaceConfig())
    D = length(x_G)
    all_modes = rigid_body_modes(D)
    nfull = length(all_modes)
    length(x_ref) == D ||
        throw(ArgumentError("RigidBody: x_ref has length $(length(x_ref)), x_G has length $D"))
    mass > 0 || throw(ArgumentError("RigidBody: mass must be positive, got $mass"))
    Iner = inertia isa Real ? fill(Float64(inertia), 1, 1) : Matrix{Float64}(inertia)
    size(Iner) == (D == 2 ? (1, 1) : (3, 3)) ||
        throw(ArgumentError("RigidBody: inertia must be a scalar I_yy in 2D or a 3×3 matrix in 3D"))
    modes = Tuple(Symbol.(modes))
    isempty(modes) && throw(ArgumentError("RigidBody: no active modes"))
    allunique(modes) || throw(ArgumentError("RigidBody: repeated modes in $modes"))
    for m in modes
        m in all_modes ||
            throw(ArgumentError("RigidBody: unknown mode :$m for ambient_dim = $D (valid: $all_modes)"))
    end
    _full(A, name) = isnothing(A) ? zeros(nfull, nfull) :
        (size(A) == (nfull, nfull) ? Matrix{Float64}(A) :
         throw(ArgumentError("RigidBody: $name must be $nfull×$nfull (all modes), got $(size(A))")))
    K = _full(K_ext, "K_ext")
    B = _full(B_ext, "B_ext")
    hs = hydrostatics isa Symbol ? hydrostatics : _full(hydrostatics, "hydrostatics")
    if hs isa Symbol
        hs in (:zero_draft, :none) ||
            throw(ArgumentError("RigidBody: hydrostatics must be :zero_draft, :none or a matrix, got :$hs"))
    end
    body = RigidBody(Float64(mass), Iner, Vector{Float64}(x_G), Vector{Float64}(x_ref), D,
                     modes, hs, K, B, Float64(ρw), Float64(g), symbol, space_domain_symbol, fe)
    _check_restrained(body)
    return body
end

# Zero-draft bodies have no normal velocity in the horizontal modes, so these
# modes need an external restraint.
function _check_restrained(b::RigidBody)
    b.hydrostatics === :zero_draft || return nothing
    free = b.ambient_dim == 2 ? (:surge,) : (:surge, :sway, :yaw)
    all_modes = rigid_body_modes(b.ambient_dim)
    for m in b.modes
        m in free || continue
        j = findfirst(==(m), all_modes)
        (iszero(b.K_ext[j, :]) && iszero(b.K_ext[:, j])) && throw(ArgumentError(
            "RigidBody: mode :$m of a zero-draft body has no restoring and no fluid " *
            "coupling; restrain it through K_ext or remove it from `modes`"))
    end
    return nothing
end

variable_symbol(b::RigidBody) = b.symbol
ambient_dimension(b::RigidBody) = b.ambient_dim

"""
    num_modes(b::AbstractRigidBody) -> Int

Number of active modes (the size of the body unknown).
"""
num_modes(b::RigidBody) = length(b.modes)

# Indices of the active modes in the full mode list.
_mode_indices(b::RigidBody) =
    [findfirst(==(m), rigid_body_modes(b.ambient_dim)) for m in b.modes]

# Restrict a full-mode matrix to the active modes.
_active_block(b::RigidBody, A::AbstractMatrix) = A[_mode_indices(b), _mode_indices(b)]

function print_parameters(b::RigidBody)
    @printf("\n[MSG] Rigid body :%s on :%s (%dD)\n", b.symbol, b.space_domain_symbol, b.ambient_dim)
    println("[VAL] modes = ", b.modes)
    @printf("[VAL] mass = %.4f kg\n", b.mass)
    println("[VAL] x_G = ", b.x_G, ", x_ref = ", b.x_ref)
    println("[VAL] hydrostatics = ", b.hydrostatics isa Symbol ? b.hydrostatics : "<matrix>")
    println()
end

# ── Kinematics ───────────────────────────────────────────────

_cross_matrix(r) = [0.0 -r[3] r[2]; r[3] 0.0 -r[1]; -r[2] r[1] 0.0]

"""
    rigid_body_mass_matrix(b::RigidBody) -> Matrix{Float64}

Full-mode rigid-body mass matrix about `b.x_ref` (dimensional):

```math
M_{rb} = \\begin{bmatrix} m I & -m [r_G]_\\times \\\\ m [r_G]_\\times & I_G + m(|r_G|^2 I - r_G r_G^T) \\end{bmatrix},
\\quad r_G = x_G - x_r .
```

In 2D this is the 3D matrix restricted to (surge, heave, pitch):
`[m 0 m r_z; 0 m −m r_x; m r_z −m r_x I_yy + m |r_G|²]`.
"""
function rigid_body_mass_matrix(b::RigidBody)
    m = b.mass
    if b.ambient_dim == 2
        rx, rz = b.x_G .- b.x_ref
        return [m 0.0 m*rz;
                0.0 m -m*rx;
                m*rz -m*rx b.inertia[1, 1] + m * (rx^2 + rz^2)]
    end
    r = b.x_G .- b.x_ref
    Rx = _cross_matrix(r)
    E = _identity3()
    Iref = b.inertia + m * (sum(abs2, r) * E - r * r')
    return [m*E -m*Rx; m*Rx Iref]
end

_identity3() = [i == j ? 1.0 : 0.0 for i in 1:3, j in 1:3]

"""
    rigid_body_displacement(b::RigidBody, ξ, x) -> Vector

Displacement `u(x) = ξ_t + α × (x − x_r)` of the body point `x` for the
active-mode amplitudes `ξ` (ordered as `b.modes`).
"""
function rigid_body_displacement(b::RigidBody, ξ::AbstractVector, x::AbstractVector)
    q = zeros(eltype(ξ), length(rigid_body_modes(b.ambient_dim)))
    q[_mode_indices(b)] .= ξ
    r = x .- b.x_ref
    if b.ambient_dim == 2
        return [q[1] + q[3] * r[2], q[2] - q[3] * r[1]]
    end
    α = q[4:6]
    return q[1:3] .+ [α[2] * r[3] - α[3] * r[2], α[3] * r[1] - α[1] * r[3], α[1] * r[2] - α[2] * r[1]]
end

# Full-mode generalized normal (tuple) at x with normal n.
_generalized_normal_full(x_r, x, n, ::Val{2}) =
    (n[1], n[2], (x[2] - x_r[2]) * n[1] - (x[1] - x_r[1]) * n[2])
function _generalized_normal_full(x_r, x, n, ::Val{3})
    r = (x[1] - x_r[1], x[2] - x_r[2], x[3] - x_r[3])
    (n[1], n[2], n[3],
     r[2] * n[3] - r[3] * n[2], r[3] * n[1] - r[1] * n[3], r[1] * n[2] - r[2] * n[1])
end

"""
    generalized_normal(b::RigidBody) -> Function

`(x, n) -> VectorValue` giving the generalized normal `N` of the active modes
at the point `x` with fluid normal `n`, so that `u·n = ξ·N`.  Use as
`Operation(generalized_normal(b))(x_Γ, n_Γ)`.
"""
function generalized_normal(b::RigidBody)
    x_r = Tuple(b.x_ref)
    idx = Tuple(_mode_indices(b))
    D = Val(b.ambient_dim)
    (x, n) -> begin
        Nf = _generalized_normal_full(x_r, x, n, D)
        VectorValue(map(i -> Nf[i], idx))
    end
end

"""
    vertical_shape(b::RigidBody) -> Function

`x -> VectorValue` with `ξ ⋅ vertical_shape(b)(x) = e_z ⋅ u(x)`, the vertical
displacement of the body point `x`.
"""
function vertical_shape(b::RigidBody)
    x_r = Tuple(b.x_ref)
    idx = Tuple(_mode_indices(b))
    if b.ambient_dim == 2
        x -> (Uf = (0.0, 1.0, -(x[1] - x_r[1])); VectorValue(map(i -> Uf[i], idx)))
    else
        # u_z = ξ₃ + α₁ r_y − α₂ r_x
        x -> (Uf = (0.0, 0.0, 1.0, x[2] - x_r[2], -(x[1] - x_r[1]), 0.0);
              VectorValue(map(i -> Uf[i], idx)))
    end
end

# Generalized normal as a CellField on the wetted surface.
function _generalized_normal_field(b::RigidBody, dom::IntegrationDomains)
    dΓ = _space_measure(dom, b)
    Γ = get_triangulation(dΓ.quad)
    Operation(generalized_normal(b))(get_physical_coordinate(Γ), get_normal_vector(Γ))
end

# ∫ ζ·(A ξ) / |Γ_b| dΓ_b  — assembles the active-mode matrix A exactly.
function _constant_matrix_form(b::RigidBody, dom::IntegrationDomains, A::AbstractMatrix, x, y)
    n = num_modes(b)
    At = TensorValue{n,n}(Tuple(Float64.(A)))
    dΓ = _space_measure(dom, b)
    area = sum(∫(1.0)dΓ)
    ξ = x[variable_symbol(b)]
    ζ = y[variable_symbol(b)]
    ∫(ζ ⋅ (At ⋅ ξ) / area)dΓ
end

# ── Weak forms ───────────────────────────────────────────────

has_mass_form(::RigidBody) = true
has_damping_form(b::RigidBody) = !iszero(b.B_ext)
has_stiffness_form(::RigidBody) = true
has_rhs_form(::RigidBody) = false

"""
    mass(b::RigidBody, dom::IntegrationDomains, x_tt, y)

Rigid-body inertia `ζ·(M̂ ξ̈)`, `M̂ = M_rb / ρw` (see
[`rigid_body_mass_matrix`](@ref)), restricted to the active modes.
"""
mass(b::RigidBody, dom::IntegrationDomains, x_tt, y) =
    _constant_matrix_form(b, dom, _active_block(b, rigid_body_mass_matrix(b)) / b.ρw, x_tt, y)

"""
    damping(b::RigidBody, dom::IntegrationDomains, x_t, y)

External linear damping `ζ·(B̂ ξ̇)`, `B̂ = B_ext / ρw`.
"""
damping(b::RigidBody, dom::IntegrationDomains, x_t, y) =
    _constant_matrix_form(b, dom, _active_block(b, b.B_ext) / b.ρw, x_t, y)

"""
    stiffness(b::RigidBody, dom::IntegrationDomains, x, y)

Hydrostatic and external stiffness.  With `hydrostatics = :zero_draft`, the
restoring term is `∫_Γb g (ζ·N)(ξ·N) dΓ`; with a matrix `C_hs` it is
`ζ·(C_hs ξ)/ρw`.  `K_ext/ρw` is added in all cases.
"""
function stiffness(b::RigidBody, dom::IntegrationDomains, x, y)
    K = b.hydrostatics isa Matrix ? b.hydrostatics + b.K_ext : b.K_ext
    val = _constant_matrix_form(b, dom, _active_block(b, K) / b.ρw, x, y)
    if b.hydrostatics === :zero_draft
        N = _generalized_normal_field(b, dom)
        ξ = x[variable_symbol(b)]
        ζ = y[variable_symbol(b)]
        val = val + ∫(b.g * (ζ ⋅ N) * (ξ ⋅ N))_space_measure(dom, b)
    end
    return val
end

"""
    restoring_matrix(b::RigidBody, dom::IntegrationDomains) -> Matrix{Float64}

Dimensional restoring matrix of the active modes: hydrostatics
(`ρw g ∫_Γb N Nᵀ dΓ` for `:zero_draft`, or the user matrix) plus `K_ext`.
"""
function restoring_matrix(b::RigidBody, dom::IntegrationDomains)
    K = _active_block(b, b.hydrostatics isa Matrix ? b.hydrostatics + b.K_ext : b.K_ext)
    b.hydrostatics === :zero_draft || return K
    N = _generalized_normal_field(b, dom)
    dΓ = _space_measure(dom, b)
    n = num_modes(b)
    e(i) = VectorValue(ntuple(k -> float(k == i), n))
    C = [sum(∫((N ⋅ e(i)) * (N ⋅ e(j)))dΓ) for i in 1:n, j in 1:n]
    return K + b.ρw * b.g * C
end

# ==========================================================================
# PrescribedRigidMotion — forced motion of a rigid body (radiation problem)
# ==========================================================================

"""
    PrescribedRigidMotion(body; mode = nothing, amplitudes = nothing)

Frequency-domain forced motion of `body` with complex displacement
amplitudes `ξ̂` (ordered as `body.modes`; `mode = :heave` is a unit amplitude
in that mode).  The body has no unknowns: the entity only adds the Neumann
data `∂ₙϕ = ∂tξ̂·N` on the wetted surface to the fluid equation, i.e. the
right-hand side `−iω ∫_Γb w (ξ̂·N) dΓ` (time dependence `e^{−iωt}`).

Used to solve radiation problems for added mass and damping (see
`PostProcessing.radiation_coefficients`).  The diffraction problem needs no
entity: leave the body out of the physics and its wetted surface is a fixed
wall.
"""
struct PrescribedRigidMotion{B<:AbstractRigidBody} <: PhysicsParameters
    body::B
    amplitudes::Vector{ComplexF64}
    fe::Nothing
end

function PrescribedRigidMotion(body::AbstractRigidBody; mode = nothing, amplitudes = nothing)
    (isnothing(mode) == isnothing(amplitudes)) &&
        throw(ArgumentError("PrescribedRigidMotion: give exactly one of `mode` and `amplitudes`"))
    if !isnothing(mode)
        j = findfirst(==(mode), body.modes)
        isnothing(j) && throw(ArgumentError("PrescribedRigidMotion: :$mode is not an active mode of the body $(body.modes)"))
        amplitudes = [i == j ? 1.0 + 0im : 0.0 + 0im for i in 1:num_modes(body)]
    end
    length(amplitudes) == num_modes(body) ||
        throw(ArgumentError("PrescribedRigidMotion: need $(num_modes(body)) amplitudes, got $(length(amplitudes))"))
    PrescribedRigidMotion(body, Vector{ComplexF64}(amplitudes), nothing)
end

function print_parameters(pm::PrescribedRigidMotion)
    @printf("\n[MSG] Prescribed motion of rigid body :%s\n", variable_symbol(pm.body))
    println("[VAL] ", NamedTuple{pm.body.modes}(Tuple(pm.amplitudes)))
end

# Field-less entity: no unknowns, only the fluid right-hand side.
variable_symbol(::PrescribedRigidMotion) = nothing
variable_symbols(::PrescribedRigidMotion) = ()
field_fe_configs(::PrescribedRigidMotion) = ()
ambient_dimension(pm::PrescribedRigidMotion) = pm.body.ambient_dim
has_mass_form(::PrescribedRigidMotion) = false
has_damping_form(::PrescribedRigidMotion) = false
has_stiffness_form(::PrescribedRigidMotion) = false
has_rhs_form(::PrescribedRigidMotion) = false
