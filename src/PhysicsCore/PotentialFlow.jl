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

variable_symbol(::PotentialFlow) = :ϕ

# ── Single-variable weak forms ─────────────────────────────
#    Field access via variable_symbol (velocity potential)

function mass(pf::PotentialFlow, dom::WeakFormDomains, x_tt, y)
    sym = variable_symbol(pf)
    w = y[sym]
    ∫(0.0 * w)dom[:dΩ]
end

function damping(pf::PotentialFlow, dom::WeakFormDomains, x_t, y)
    sym = variable_symbol(pf)
    w = y[sym]
    ∫(0.0 * w)dom[:dΩ]
end

function stiffness(pf::PotentialFlow, dom::WeakFormDomains, x, y)
    sym = variable_symbol(pf)
    ϕ = x[sym]
    w = y[sym]
    ∫(∇(w) ⋅ ∇(ϕ))dom[:dΩ]
end

function rhs(pf::PotentialFlow, dom::WeakFormDomains, f, y)
    sym = variable_symbol(pf)
    w = y[sym]
    ∫(w * f[sym])dom[:dΩ]
end
