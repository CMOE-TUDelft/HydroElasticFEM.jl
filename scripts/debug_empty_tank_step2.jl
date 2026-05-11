"""
    debug_empty_tank_step2.jl

Step 2: Compare assembled stiffness matrices and load vectors between
plain and structured implementations to pinpoint the source of the
wave-generation mismatch.

Hypothesis: despite identical math on paper, the assembled systems may differ
due to (a) the SubFacetTriangulation used for Γκ in structured vs the plain
BoundaryTriangulation, (b) DOF ordering differences, or (c) some subtle form
evaluation difference.

What this script checks:
  A) Matrix/RHS norms — confirms whether the assembled systems are numerically equal
  B) Matrix block-by-block — Laplace block, free-surface block, coupling block
  C) RHS vector comparison
  D) Forces the structured problem to use the IDENTICAL triangulations as the plain
     problem and checks if that resolves the mismatch (isolates triangulation type)
"""

using LinearAlgebra
using Printf
using Gridap
using WaveSpec
using Parameters

using HydroElasticFEM: PKG_ROOT, map_vertical_GP_for_const_dep
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S

# ── Import helpers from the example module ────────────────────────────────
push!(LOAD_PATH, joinpath(PKG_ROOT, "examples"))
include(joinpath(PKG_ROOT, "examples", "EmptyTankExample.jl"))
using .EmptyTankExample: EmptyTankTutorialParams, tank_parameters, gp_map,
                         shifted_gp_map, incident_wave, build_regular_wave_state,
                         probe_points

# ── Test parameters ───────────────────────────────────────────────────────
# Use a small mesh so matrices are manageable
const NX = 12
const NY = 4
const ORDER = 1

println("=" ^ 70)
@printf("Step 2 — Matrix comparison (nx=%d, ny=%d, order=%d)\n", NX, NY, ORDER)
println("=" ^ 70)

# ─────────────────────────────────────────────────────────────────────────
# Build plain problem
# ─────────────────────────────────────────────────────────────────────────
@printf("--- Building PLAIN problem ---\n")

p    = EmptyTankTutorialParams(nx=NX, ny=NY, order=ORDER)
tp   = tank_parameters(; H0=p.H0, nx=p.nx, ny=p.ny)
inc  = incident_wave(; H0=p.H0, ω=p.ω, η0=p.η0, α=p.α)

model_plain = CartesianDiscreteModel(tp.domain, tp.partition,
                                     map=gp_map(p.mesh_ry, p.ny, p.H0))

labels_plain = get_face_labeling(model_plain)
add_tag_from_tags!(labels_plain, "surface", [3, 4, 6])
add_tag_from_tags!(labels_plain, "bottom",  [1, 2, 5])
add_tag_from_tags!(labels_plain, "inlet",   [7])
add_tag_from_tags!(labels_plain, "outlet",  [8])
add_tag_from_tags!(labels_plain, "water",   [9])

Ω_plain   = Interior(model_plain)
Γfs_plain = Boundary(model_plain, tags="surface")
Γin_plain = Boundary(model_plain, tags="inlet")
Γout_plain= Boundary(model_plain, tags="outlet")

degree_plain = 2 * p.order
dΩ_p   = Measure(Ω_plain, degree_plain)
dΓfs_p = Measure(Γfs_plain, degree_plain)
dΓin_p = Measure(Γin_plain, degree_plain)
dΓout_p= Measure(Γout_plain, degree_plain)
nΓin_p = get_normal_vector(Γin_plain)

reffe_plain = ReferenceFE(lagrangian, Float64, p.order)
VΩ_p    = TestFESpace(Ω_plain,  reffe_plain; conformity=:H1, vector_type=Vector{ComplexF64})
VΓκ_p   = TestFESpace(Γfs_plain,reffe_plain; conformity=:H1, vector_type=Vector{ComplexF64})
UΩ_p    = TrialFESpace(VΩ_p)
UΓκ_p   = TrialFESpace(VΓκ_p)
X_plain = MultiFieldFESpace([UΩ_p,  UΓκ_p])
Y_plain = MultiFieldFESpace([VΩ_p,  VΓκ_p])

αₕ_plain = -im * p.ω / WaveSpec.PhysicalConstants.g * (1.0 - p.βₕ) / p.βₕ
k_plain  = inc.sea_state.k[1]
g_plain  = WaveSpec.PhysicalConstants.g

a_plain((ϕ, κ), (w, u)) =
    ∫(∇(w) ⋅ ∇(ϕ))dΩ_p +
    ∫(p.βₕ * (u + αₕ_plain * w) * (g_plain * κ - im * p.ω * ϕ) + im * p.ω * w * κ)dΓfs_p +
    ∫(-im * k_plain * w * ϕ)dΓin_p +
    ∫(-im * k_plain * w * ϕ)dΓout_p

l_plain((w, u)) =
    ∫(w * (inc.vin ⋅ nΓin_p))dΓin_p -
    ∫(im * k_plain * w * inc.ϕin)dΓin_p

op_plain = AffineFEOperator(a_plain, l_plain, X_plain, Y_plain)
A_plain  = get_matrix(op_plain)
b_plain  = get_vector(op_plain)

@printf("  Plain DOFs: %d\n", size(A_plain, 1))

# ─────────────────────────────────────────────────────────────────────────
# Build structured problem (default)
# ─────────────────────────────────────────────────────────────────────────
@printf("\n--- Building STRUCTURED problem ---\n")

sea_state = build_regular_wave_state(H=2.0*p.η0, T=2π/p.ω, h=p.H0)
f_in(x)   = (inc.vin(x) ⋅ VectorValue(-1.0, 0.0)) + (-1.0) * im * inc.sea_state.k[1] * inc.ϕin(x)

tank = G.TankDomain(
    L  = tp.LΩ,
    H  = p.H0,
    nx = p.nx,
    ny = p.ny,
    map= shifted_gp_map(tp.x0, p.mesh_ry, p.ny, p.H0),
)

p_flow = P.PotentialFlow(
    ρw  = 1025.0,
    g   = WaveSpec.PhysicalConstants.g,
    sea_state = sea_state,
    boundary_conditions = [
        P.RadiationBC(domain=:dΓin),
        P.RadiationBC(domain=:dΓout),
        P.PrescribedInletPotentialBC(domain=:dΓin, forcing=f_in, quantity=:traction),
    ],
    fe   = PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64}),
    space_domain_symbol = :Ω,
)

free_surface = P.FreeSurface(
    ρw  = 1025.0,
    g   = WaveSpec.PhysicalConstants.g,
    βₕ  = p.βₕ,
    fe  = PH.FESpaceConfig(order=p.order, vector_type=Vector{ComplexF64}),
    space_domain_symbol = :Γκ,
)

config  = S.FreqDomainConfig(ω=p.ω)
problem = S.build_problem(tank, P.PhysicsParameters[p_flow, free_surface], config)

op_struct  = S.get_fe_operator(problem)
A_struct   = get_matrix(op_struct)
b_struct   = get_vector(op_struct)

@printf("  Structured DOFs: %d\n", size(A_struct, 1))

# ─────────────────────────────────────────────────────────────────────────
# A — Global matrix / RHS comparison
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION A — Global system comparison\n")
println("=" ^ 70)

if size(A_plain) != size(A_struct)
    @printf("  [FAIL] Matrix sizes DIFFER: plain=%s  structured=%s\n",
            string(size(A_plain)), string(size(A_struct)))
else
    diff_A    = A_struct - A_plain
    diff_b    = b_struct - b_plain
    relA = norm(diff_A) / norm(A_plain)
    relb = length(b_plain) > 0 ? norm(diff_b) / (norm(b_plain) + eps()) : NaN
    @printf("  ‖A_struct - A_plain‖_F  = %.3e  (relative: %.3e)\n", norm(diff_A), relA)
    @printf("  ‖b_struct - b_plain‖_2  = %.3e  (relative: %.3e)\n", norm(diff_b), relb)

    if relA < 1e-8 && relb < 1e-8
        @printf("  [PASS] Assembled systems are numerically identical — bug is NOT in assembly\n")
    else
        @printf("  [FAIL] Assembled systems DIFFER — bug is in form assembly!\n")
    end
end

# ─────────────────────────────────────────────────────────────────────────
# B — RHS vector element-by-element
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION B — RHS vector inspection\n")
println("=" ^ 70)

n_dof = min(length(b_plain), length(b_struct))
n_print = min(20, n_dof)

@printf("  First %d entries of b_plain vs b_struct:\n", n_print)
@printf("  %5s  %25s  %25s  %12s\n", "DOF", "b_plain", "b_struct", "diff")
for i in 1:n_print
    bp = b_plain[i]
    bs = b_struct[i]
    @printf("  %5d  %25s  %25s  %12.3e\n",
        i, string(round(bp, digits=6)), string(round(bs, digits=6)), abs(bs - bp))
end

# Check which DOFs have nonzero entries in b (these are the inlet DOFs)
nz_plain  = findall(!iszero, b_plain)
nz_struct = findall(!iszero, b_struct)
@printf("\n  Nonzero b_plain  entries: indices = %s\n", string(nz_plain))
@printf("  Nonzero b_struct entries: indices = %s\n", string(nz_struct))

# ─────────────────────────────────────────────────────────────────────────
# C — Solve and compare solutions
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION C — Solution comparison\n")
println("=" ^ 70)

ϕₕ_p, κₕ_p = solve(LUSolver(), op_plain)
result_s    = S.simulate(problem)
ϕₕ_s, κₕ_s = result_s.solution

probes = probe_points(p.probe_x)

κ_plain    = κₕ_p(probes)
κ_struct   = κₕ_s(probes)
ϕ_plain    = ϕₕ_p(probes)
ϕ_struct   = ϕₕ_s(probes)

@printf("  %5s  %25s  %25s  %12s\n", "probe", "κ_plain", "κ_struct", "|diff|")
for (i, (κp, κs)) in enumerate(zip(κ_plain, κ_struct))
    @printf("  %5d  %25s  %25s  %12.4e\n",
            i, string(round(κp, digits=4)), string(round(κs, digits=4)), abs(κp - κs))
end

@printf("\n  %5s  %25s  %25s  %12s\n", "probe", "ϕ_plain", "ϕ_struct", "|diff|")
for (i, (ϕp, ϕs)) in enumerate(zip(ϕ_plain, ϕ_struct))
    @printf("  %5d  %25s  %25s  %12.4e\n",
            i, string(round(ϕp, digits=4)), string(round(ϕs, digits=4)), abs(ϕp - ϕs))
end

# ─────────────────────────────────────────────────────────────────────────
# D — Structured using PLAIN triangulations (isolate SubFacetTriangulation)
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION D — Structured physics on PLAIN triangulations\n")
println("=" ^ 70)
@printf("  (Forces structured assembly to use Boundary(model, tags='surface')\n")
@printf("   instead of the SubFacetTriangulation Γκ = Triangulation(Γ, indices))\n\n")

# Reuse the plain model and its triangulations but assemble via structured physics

# Build integration domains from the PLAIN model  
# We need to make a fake TankTriangulations from the plain model
struct_model = model_plain  # same mesh
trians_from_plain = G.TankTriangulations(Dict{Symbol, Any}(
    :Ω   => Ω_plain,
    :Γ   => Γfs_plain,
    :Γfs => Γfs_plain,
    :Γκ  => Γfs_plain,    # Use plain BoundaryTriangulation directly (not sub-trian)
    :Γη  => Triangulation(Γfs_plain, Int[]),
    :Γin  => Γin_plain,
    :Γout => Γout_plain,
    :Γbot => Boundary(model_plain, tags="bottom"),
    :Γ_structures => Any[],
    :Γ_dampings   => Any[],
))

try
    degrees_d = Dict{Symbol, Int}(
        :dΩ => degree_plain, :dΓκ => degree_plain, :dΓfs => degree_plain,
        :dΓη => degree_plain, :dΓin => degree_plain, :dΓout => degree_plain,
        :dΓbot => degree_plain,
    )
    measures_d = G.get_integration_domains(trians_from_plain; degree=degrees_d)
    ctx_d      = S.build_frequency_context(measures_d, P.PhysicsParameters[p_flow, free_surface], config)

    X_d, Y_d, fmap_d = S.build_fe_spaces(P.PhysicsParameters[p_flow, free_surface], trians_from_plain, config)
    op_d = S.build_frequency_fe_operator(P.PhysicsParameters[p_flow, free_surface], ctx_d, fmap_d, X_d, Y_d)
    A_d  = get_matrix(op_d)
    b_d  = get_vector(op_d)

    @printf("  DOFs: %d\n", size(A_d, 1))

    if size(A_d) == size(A_plain)
        relA_d = norm(A_d - A_plain) / norm(A_plain)
        relb_d = norm(b_d - b_plain) / (norm(b_plain) + eps())
        @printf("  ‖A_plain_trian - A_plain‖_F  = %.3e (relative: %.3e)\n", norm(A_d - A_plain), relA_d)
        @printf("  ‖b_plain_trian - b_plain‖_2  = %.3e (relative: %.3e)\n", norm(b_d - b_plain), relb_d)
        if relA_d < 1e-8 && relb_d < 1e-8
            @printf("  [PASS] Plain triangulation + structured physics = identical to plain!\n")
            @printf("         => Triangulation type is not introducing measurable differences here.\n")
        else
            @printf("  [FAIL] Still differs even with plain triangulations\n")
            @printf("         => Bug is in the structured FORM ASSEMBLY, not the triangulation type\n")
        end
    else
        @printf("  DOF count differs from plain — cannot compare directly\n")
    end
catch e
    @printf("  [ERROR] Section D failed: %s\n", string(e))
    @printf("  (This may fail if TankTriangulations does not accept manual construction)\n")
end

# ─────────────────────────────────────────────────────────────────────────
# E — Diagnose stiffness sub-blocks
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION E — Individual form contributions\n")
println("=" ^ 70)
@printf("  Compute each contribution to plain a(...) separately and compare to structured\n\n")

n_Ω   = num_free_dofs(VΩ_p)
n_Γκ  = num_free_dofs(VΓκ_p)
@printf("  n_Ω (ϕ DOFs) = %d,  n_Γκ (κ DOFs) = %d\n", n_Ω, n_Γκ)

# Laplace block only
a_lap((ϕ, κ), (w, u)) = ∫(∇(w) ⋅ ∇(ϕ))dΩ_p
op_lap   = AffineFEOperator(a_lap, l_plain, X_plain, Y_plain)
A_lap    = get_matrix(op_lap)

# Free-surface block only
a_fs((ϕ, κ), (w, u)) = ∫(p.βₕ * (u + αₕ_plain * w) * (g_plain * κ - im * p.ω * ϕ) + im * p.ω * w * κ)dΓfs_p
op_fs  = AffineFEOperator(a_fs, l_plain, X_plain, Y_plain)
A_fs   = get_matrix(op_fs)

# Radiation BCs only
a_rad((ϕ, κ), (w, u)) = ∫(-im * k_plain * w * ϕ)dΓin_p + ∫(-im * k_plain * w * ϕ)dΓout_p
op_rad = AffineFEOperator(a_rad, l_plain, X_plain, Y_plain)
A_rad  = get_matrix(op_rad)

@printf("  ‖A_lap‖_F  = %.6e\n", norm(A_lap))
@printf("  ‖A_fs‖_F   = %.6e\n", norm(A_fs))
@printf("  ‖A_rad‖_F  = %.6e\n", norm(A_rad))
@printf("  ‖A_lap + A_fs + A_rad - A_plain‖_F = %.3e  (should be ~0)\n",
        norm(A_lap + A_fs + A_rad - A_plain))

# ─────────────────────────────────────────────────────────────────────────
# F — Check DOF ordering between plain and structured FE spaces
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION F — DOF count and diagonal comparison\n")
println("=" ^ 70)

trians_struct = S.get_triangulations(problem)
Ω_s = trians_struct[:Ω]
Γκ_s = trians_struct[:Γκ]

@printf("  Plain:      Ω cells = %d, Γκ cells = %d\n",
        num_cells(Ω_plain), num_cells(Γfs_plain))
@printf("  Structured: Ω cells = %d, Γκ cells = %d\n",
        num_cells(Ω_s), num_cells(Γκ_s))

# Diagonal of A — easiest to compare (permutation-invariant up to DOF relabelling)
diag_plain  = sort(real(diag(A_plain)))
diag_struct = sort(real(diag(A_struct)))

if length(diag_plain) == length(diag_struct)
    diff_diag = norm(diag_struct - diag_plain) / (norm(diag_plain) + eps())
    @printf("  ‖sort(diag(A_struct)) - sort(diag(A_plain))‖ / ‖diag(A_plain)‖ = %.3e\n", diff_diag)
    if diff_diag < 1e-8
        @printf("  [INFO] Diagonal spectra match — DOF count and distribution consistent\n")
    else
        @printf("  [INFO] Diagonal spectra DIFFER — DOF content or distribution different\n")
        n_show = min(15, length(diag_plain))
        @printf("  Plain  diagonal (sorted, first %d): %s\n", n_show,
                string(round.(diag_plain[1:n_show], digits=4)))
        @printf("  Struct diagonal (sorted, first %d): %s\n", n_show,
                string(round.(diag_struct[1:n_show], digits=4)))
    end
else
    @printf("  [FAIL] Diagonal lengths differ: plain=%d, structured=%d\n",
            length(diag_plain), length(diag_struct))
end

# ─────────────────────────────────────────────────────────────────────────
# G — Block-level comparison (ϕ/κ split)
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION G — Block-level matrix mismatch\n")
println("=" ^ 70)

Iϕ = 1:n_Ω
Iκ = (n_Ω + 1):(n_Ω + n_Γκ)

A11_p = A_plain[Iϕ, Iϕ]
A12_p = A_plain[Iϕ, Iκ]
A21_p = A_plain[Iκ, Iϕ]
A22_p = A_plain[Iκ, Iκ]

A11_s = A_struct[Iϕ, Iϕ]
A12_s = A_struct[Iϕ, Iκ]
A21_s = A_struct[Iκ, Iϕ]
A22_s = A_struct[Iκ, Iκ]

rel11 = norm(A11_s - A11_p) / (norm(A11_p) + eps())
rel12 = norm(A12_s - A12_p) / (norm(A12_p) + eps())
rel21 = norm(A21_s - A21_p) / (norm(A21_p) + eps())
rel22 = norm(A22_s - A22_p) / (norm(A22_p) + eps())

@printf("  rel diff A11 (ϕ,ϕ): %.3e\n", rel11)
@printf("  rel diff A12 (ϕ,κ): %.3e\n", rel12)
@printf("  rel diff A21 (κ,ϕ): %.3e\n", rel21)
@printf("  rel diff A22 (κ,κ): %.3e\n", rel22)

max_block = maximum([rel11, rel12, rel21, rel22])
if max_block == rel11
    @printf("  [INFO] Dominant mismatch block: A11 (PotentialFlow core)\n")
elseif max_block == rel12
    @printf("  [INFO] Dominant mismatch block: A12 (PF/FS coupling)\n")
elseif max_block == rel21
    @printf("  [INFO] Dominant mismatch block: A21 (PF/FS coupling)\n")
else
    @printf("  [INFO] Dominant mismatch block: A22 (FreeSurface core)\n")
end

# ─────────────────────────────────────────────────────────────────────────
# H — PotentialFlow-only check (remove FreeSurface field)
# ─────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 70)
@printf("SECTION H — PotentialFlow-only operator comparison\n")
println("=" ^ 70)

# Plain PF-only
VΩ_pf = TestFESpace(Ω_plain, reffe_plain; conformity=:H1, vector_type=Vector{ComplexF64})
UΩ_pf = TrialFESpace(VΩ_pf)

a_pf_plain(ϕ, w) =
    ∫(∇(w) ⋅ ∇(ϕ))dΩ_p +
    ∫(-im * k_plain * w * ϕ)dΓin_p +
    ∫(-im * k_plain * w * ϕ)dΓout_p

l_pf_plain(w) =
    ∫(w * (inc.vin ⋅ nΓin_p))dΓin_p -
    ∫(im * k_plain * w * inc.ϕin)dΓin_p

op_pf_plain = AffineFEOperator(a_pf_plain, l_pf_plain, UΩ_pf, VΩ_pf)
A_pf_plain = get_matrix(op_pf_plain)
b_pf_plain = get_vector(op_pf_plain)

# Structured PF-only
problem_pf = S.build_problem(tank, P.PhysicsParameters[p_flow], config)
op_pf_struct = S.get_fe_operator(problem_pf)
A_pf_struct = get_matrix(op_pf_struct)
b_pf_struct = get_vector(op_pf_struct)

if size(A_pf_struct) == size(A_pf_plain)
    relA_pf = norm(A_pf_struct - A_pf_plain) / (norm(A_pf_plain) + eps())
    relb_pf = norm(b_pf_struct - b_pf_plain) / (norm(b_pf_plain) + eps())
    @printf("  rel diff A_pf_only: %.3e\n", relA_pf)
    @printf("  rel diff b_pf_only: %.3e\n", relb_pf)
    if relA_pf < 1e-8
        @printf("  [INFO] PF-only stiffness matches plain.\n")
    else
        @printf("  [INFO] PF-only stiffness already differs. Mismatch is inside PotentialFlow stiffness assembly.\n")
    end
else
    @printf("  [FAIL] PF-only matrix size mismatch: plain=%s structured=%s\n",
            string(size(A_pf_plain)), string(size(A_pf_struct)))
end

println()
println("=" ^ 70)
@printf("Done.\n")
