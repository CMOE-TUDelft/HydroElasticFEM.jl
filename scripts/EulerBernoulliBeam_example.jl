module EulerBernoulliBeamExample

using HydroElasticFEM
using Gridap

# ─────────────────────────────────────────────────────────────
# Without HydroElasticFEM — raw Gridap
# ─────────────────────────────────────────────────────────────
#
# model = CartesianDiscreteModel((0,1),(10,))
# Ω  = Triangulation(model)
# Λ  = Skeleton(Ω)
# dΩ = Measure(Ω, 4)
# dΛ = Measure(Λ, 4)
# nΛ = get_normal_vector(Λ)
# reffe = ReferenceFE(lagrangian, Float64, 2)
# V = TestFESpace(model, reffe, dirichlet_tags="boundary")
# U = TrialFESpace(V, 0.0)
# f(x) = 1.0
# EI = 100.0
# γ_m = 10.0
# h  = 0.1
# a(u,v) = ∫(EI * Δ(v) * Δ(u))dΩ +
#           ∫(EI * (- jump(∇(v) ⋅ nΛ) * mean(Δ(u))
#                   - mean(Δ(v)) * jump(∇(u) ⋅ nΛ)
#                   + (γ_m / h) * jump(∇(v) ⋅ nΛ) * jump(∇(u) ⋅ nΛ)))dΛ
# l(v) = ∫(f * v)dΩ
# op = AffineFEOperator(a, l, U, V)
# uh = solve(LUSolver(), op)

# ─────────────────────────────────────────────────────────────
# With HydroElasticFEM — new weak form structure
# ─────────────────────────────────────────────────────────────

# 1. Define physics parameters
beam = EulerBernoulliBeam(L=1.0, m=1.0, E=100.0, I=1.0, ρw=1.0, g=0.0,
                          bndType=FixedBoundary())
print_parameters(beam)

# 2. Build Gridap model and FE spaces
nel = 40
h   = beam.L / nel
model = CartesianDiscreteModel((0, beam.L), (nel,))
Ω  = Triangulation(model)
Λ  = Skeleton(Ω)
Λb = Boundary(model, tags="boundary")

order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe, dirichlet_tags="boundary")
U = TrialFESpace(V, 0.0)

# Wrap in MultiFieldFESpace so Gridap passes 1-tuples to closures,
# which lets FieldDict index correctly: x[1] → first (only) field.
Y = MultiFieldFESpace([V])
X = MultiFieldFESpace([U])

# 3. Populate WeakFormDomains with measures, normals, Nitsche params
dom = WeakFormDomains(
    dΓ_s   = Measure(Ω, 2 * order + 2),
    dΛ_s   = Measure(Λ, 2 * order + 2),
    n_Λ_s  = get_normal_vector(Λ),
    h_s    = h,
    γ_s    = 10.0 * order^2,
    dΛ_sb  = Measure(Λb, 2 * order + 2),
    n_Λ_sb = get_normal_vector(Λb),
)

# 4. Build the bilinear and linear forms using the weak form interface.
#    The FieldDict maps the beam's symbol (:η_b) to position 1 in the
#    single-field tuple that Gridap passes to our closures.
sym  = variable_symbol(beam)          # :η_b
fmap = Dict(sym => 1)

a((u,), (v,)) = stiffness(beam, dom,
                           FieldDict((u,), fmap), FieldDict((v,), fmap))

f(x) = 1.0
l((v,)) = rhs(beam, dom,
              FieldDict((f,), fmap), FieldDict((v,), fmap))

# 5. Solve
op = AffineFEOperator(a, l, X, Y)
uh = solve(LUSolver(), op)
writevtk(Ω, "scripts/EulerBernoulliBeam_example", nsubcells=10,
         cellfields=["uh" => uh[1]])

end # module
