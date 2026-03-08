"""
    PotentialFlow <: PhysicsParameters

Parameters for the 2D fluid potential (Laplace equation).

# Fields
- `ρw::Float64` — Density of water
- `g::Float64`  — Gravitational acceleration
"""
@with_kw struct PotentialFlow <: PhysicsParameters
    ρw::Float64 = 1025.0
    g::Float64  = 9.81
end

function print_parameters(f::PotentialFlow)
    @printf("\n[MSG] Fluid Properties:\n")
    @printf("[VAL] ρw = %.2f kg/m3\n", f.ρw)
    @printf("[VAL] g  = %.4f m/s2\n", f.g)
    println()
end

# ── Single-variable weak forms ─────────────────────────────
#    Field access via symbol :ϕ (velocity potential)

function mass(f::PotentialFlow, dom::WeakFormDomains, x_tt, y)
    w = y[:ϕ]
    ∫(0.0 * w)dom[:dΩ]
end

function damping(f::PotentialFlow, dom::WeakFormDomains, x_t, y)
    w = y[:ϕ]
    ∫(0.0 * w)dom[:dΩ]
end

function stiffness(f::PotentialFlow, dom::WeakFormDomains, x, y)
    ϕ = x[:ϕ]
    w = y[:ϕ]
    ∫(∇(w) ⋅ ∇(ϕ))dom[:dΩ]
end

function rhs(pf::PotentialFlow, dom::WeakFormDomains, f, y)
    w = y[:ϕ]
    ∫(w * f[:ϕ])dom[:dΩ]
end
