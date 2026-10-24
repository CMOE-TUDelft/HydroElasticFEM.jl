# ─────────────────────────────────────────────────────────────
# Rigid-body hydrodynamic coefficients
#
# Time dependence e^{−iωt} (the frequency-domain operator substitutes
# ∂t → −iω).  The dynamic pressure is p = −ρ ∂tϕ = iωρ ϕ, and the force on
# mode i is F_i = ∫_Γb p N_i dΓ (N from the fluid normal, which points into
# the body).  The body equation reads
#
#     [−ω²(M + A) − iω(B + B_ext) + C] ξ = X,
#
# where the radiation force of a motion ξ is F_rad = (ω² A + iω B) ξ and X is
# the force on the fixed body (diffraction problem).
# ─────────────────────────────────────────────────────────────

"""
    hydrodynamic_force(problem, result, body; fluid_symbol = :ϕ) -> Vector{ComplexF64}

Dynamic pressure force on the active modes of `body` (dimensional; per unit
width in 2D), `F_i = iωρ ∫_Γb ϕ N_i dΓ`, from a frequency-domain solve.
On a diffraction problem (body left out of the physics) this is the
excitation force; on a radiation problem ([`solve_radiation`](@ref)) the
radiation force.
"""
function hydrodynamic_force(problem, result, body::P.RigidBody; fluid_symbol::Symbol = :ϕ)
    dom = SM.get_integration_domains(problem)
    ω = SM.get_sim_config(problem).ω
    ϕh = result.solution[result.fmap[fluid_symbol]]
    N = P._generalized_normal_field(body, dom)
    dΓ = P._space_measure(dom, body)
    n = P.num_modes(body)
    e(i) = VectorValue(ntuple(k -> float(k == i), n))
    [im * ω * body.ρw * sum(∫(ϕh * (N ⋅ e(i)))dΓ) for i in 1:n]
end

"""
    excitation_force(problem, result, body; fluid_symbol = :ϕ) -> Vector{ComplexF64}

Wave excitation force on `body` from a diffraction problem, i.e. a problem
whose physics does not contain `body` (its wetted surface is then a fixed
wall).  Same as [`hydrodynamic_force`](@ref); divide by the incident
amplitude for the force per unit wave amplitude.
"""
excitation_force(problem, result, body::P.RigidBody; kwargs...) =
    hydrodynamic_force(problem, result, body; kwargs...)

"""
    solve_radiation(domain, fluid, body, config; mode) -> (problem, result)

Radiation problem: `body` forced with unit displacement amplitude in `mode`,
with the fluid entities `fluid` (e.g. `[PotentialFlow, FreeSurface]`, without
incident-wave forcing).
"""
function solve_radiation(domain, fluid::AbstractVector, body::P.RigidBody, config; mode::Symbol)
    physics = P.PhysicsParameters[fluid...; P.PrescribedRigidMotion(body; mode = mode)]
    problem = SM.build_problem(domain, physics, config)
    return problem, SM.simulate(problem)
end

"""
    radiation_coefficients(domain, fluid, body, config) -> (A, B)

Added-mass and radiation-damping matrices of the active modes of `body`
(dimensional; per unit width in 2D), from one radiation problem per mode
([`solve_radiation`](@ref)):

    F_rad,ij = ω² A_ij + iω B_ij     (unit displacement in mode j)
"""
function radiation_coefficients(domain, fluid::AbstractVector, body::P.RigidBody, config)
    ω = config.ω
    n = P.num_modes(body)
    Z = zeros(ComplexF64, n, n)
    for (j, mode) in enumerate(body.modes)
        problem, result = solve_radiation(domain, fluid, body, config; mode = mode)
        Z[:, j] = hydrodynamic_force(problem, result, body)
    end
    return real(Z) / ω^2, imag(Z) / ω
end

"""
    rigid_body_response(body, dom, ω, A, B, X) -> Vector{ComplexF64}

Motion of `body` from its hydrodynamic coefficients:
`ξ = [−ω²(M + A) − iω(B + B_ext) + C]⁻¹ X`, with `M` the rigid-body mass
matrix and `C` the restoring matrix (hydrostatics + `K_ext`) of the active
modes.  Reproduces the coupled solve exactly (linear superposition on the
same mesh).
"""
function rigid_body_response(body::P.RigidBody, dom, ω::Real, A::AbstractMatrix,
                             B::AbstractMatrix, X::AbstractVector)
    idx = P._mode_indices(body)
    M = P.rigid_body_mass_matrix(body)[idx, idx]
    C = P.restoring_matrix(body, dom)
    Bext = body.B_ext[idx, idx]
    return (-ω^2 * (M + A) - im * ω * (B + Bext) + C) \ X
end

"""
    rigid_body_dissipated_power(ξ, ω, B) -> Float64

Time-averaged power absorbed by the linear damping matrix `B` (e.g. the
active block of a body's `B_ext`, a PTO) for the complex motion `ξ`:
`P = ½ ω² Re(ξᴴ B ξ)` (per unit width in 2D).  Use with
[`absorption_from_dissipated_power`](@ref).
"""
rigid_body_dissipated_power(ξ::AbstractVector, ω::Real, B::AbstractMatrix) =
    0.5 * ω^2 * real(dot(ξ, B * ξ))
