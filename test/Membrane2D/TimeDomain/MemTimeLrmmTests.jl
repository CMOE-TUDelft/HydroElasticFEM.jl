using Gridap.ODEs
using Gridap.CellData
using WaveSpec

@testset "Time-domain membrane+LRMM" begin

  # =========================================================================
  # Parameters
  # =========================================================================

  order = 2
  ρw = 1025.0
  g = 9.81

  # Depth
  H0 = 10.0

  # Membrane (same physics as legacy, shorter length)
  Lm = 10.0
  mᵨ = 0.9
  Tᵨ = 0.1/4 * g * 20.0^2  # 98.1 — uses original Lm=20 value
  τ = 0.0

  # Resonator
  rM = 1.0e3
  rK = 5.9e3

  # Wave
  ω = 2.4
  η₀ = 0.10
  k = WaveSpec.AiryWaves.solve_wavenumber(ω, H0)
  λ = 2π / k
  T_wave = 2π / ω
  ph0 = π/2

  # Wave functions
  ηᵢₙ(x, t) = η₀ * cos(k*x[1] - ω*t + ph0)
  ϕᵢₙ(x, t) = (η₀*ω/k) * (cosh(k*(H0 + x[2])) / sinh(k*H0)) *
               sin(k*x[1] - ω*t + ph0)
  vᵢₙ(x, t) = -(η₀*ω) * (cosh(k*(H0 + x[2])) / sinh(k*H0)) *
               cos(k*x[1] - ω*t + ph0)
  vzᵢₙ(x, t) = ω * η₀ * sin(k*x[1] - ω*t + ph0)
  ηᵢₙ(t::Real) = x -> ηᵢₙ(x, t)
  ϕᵢₙ(t::Real) = x -> ϕᵢₙ(x, t)
  vᵢₙ(t::Real) = x -> vᵢₙ(x, t)
  vzᵢₙ(t::Real) = x -> vzᵢₙ(x, t)

  # Domain geometry (half-size)
  Ld = 10.0
  LΩ = 2*Ld + 3*Lm
  x₀ = -Ld
  domain = (x₀, x₀ + LΩ, -H0, 0.0)
  nx = 100
  ny = 5
  partition = (nx, ny)
  mesh_ry = 1.2

  xdᵢₙ = 0.0
  xdₒₜ = x₀ + LΩ - Ld
  xm₀ = xdᵢₙ + Lm
  xm₁ = xm₀ + Lm

  # Time stepping
  γₜ = 0.5
  βₜ = 0.25
  t₀ = 0.0
  Δt = T_wave / 20
  tf = 3 * T_wave

  ∂uₜ_∂u = γₜ / (βₜ * Δt)

  # Numerical constants
  h = LΩ / nx
  βₕ = 0.5
  αₕ = ∂uₜ_∂u / g * (1 - βₕ) / βₕ

  # Damping
  μ₀ = 2.5
  μ₁ᵢₙ(x::VectorValue) = μ₀ * (1.0 - sin(π/2 * (x[1] - x₀) / Ld))
  μ₁ₒᵤₜ(x::VectorValue) = μ₀ * (1.0 - cos(π/2 * (x[1] - xdₒₜ) / Ld))
  μ₂ᵢₙ(x) = μ₁ᵢₙ(x) * k
  μ₂ₒᵤₜ(x) = μ₁ₒᵤₜ(x) * k
  ηd(t) = x -> μ₂ᵢₙ(x) * ηᵢₙ(x, t)
  ∇ₙϕd(t) = x -> μ₁ᵢₙ(x) * vzᵢₙ(x, t)

  # =========================================================================
  # Mesh
  # =========================================================================

  mesh_map(x) = VectorValue(
    x[1],
    HydroElasticFEM.map_vertical_GP_for_const_dep(x[2], mesh_ry, ny, H0)
  )
  model = CartesianDiscreteModel(domain, partition, map=mesh_map)

  # Labelling
  labels_Ω = get_face_labeling(model)
  add_tag_from_tags!(labels_Ω, "surface", [3, 4, 6])
  add_tag_from_tags!(labels_Ω, "bottom",  [1, 2, 5])
  add_tag_from_tags!(labels_Ω, "inlet",   [7])
  add_tag_from_tags!(labels_Ω, "outlet",  [8])
  add_tag_from_tags!(labels_Ω, "water",   [9])

  # Triangulations
  Ω = Interior(model)
  Γ = Boundary(model, tags="surface")
  Γin = Boundary(model, tags="inlet")

  # Surface masking
  function is_mem(xs)
    n = length(xs); x = (1/n) * sum(xs)
    (xm₀ <= x[1] <= xm₁) * (x[2] ≈ 0.0)
  end
  function is_damping1(xs)
    n = length(xs); x = (1/n) * sum(xs)
    (x₀ <= x[1] <= xdᵢₙ) * (x[2] ≈ 0.0)
  end
  function is_damping2(xs)
    n = length(xs); x = (1/n) * sum(xs)
    (xdₒₜ <= x[1]) * (x[2] ≈ 0.0)
  end

  xΓ = get_cell_coordinates(Γ)
  Γm_mask = lazy_map(is_mem, xΓ)
  Γd1_mask = lazy_map(is_damping1, xΓ)
  Γd2_mask = lazy_map(is_damping2, xΓ)

  Γm = Triangulation(Γ, findall(Γm_mask))
  Γd1 = Triangulation(Γ, findall(Γd1_mask))
  Γd2 = Triangulation(Γ, findall(Γd2_mask))
  Γfs = Triangulation(Γ, findall(!, Γm_mask .| Γd1_mask .| Γd2_mask))
  Γη = Triangulation(Γ, findall(Γm_mask))
  Γκ = Triangulation(Γ, findall(!, Γm_mask))

  # Membrane boundary tag
  Λmb = Boundary(Γm)
  xΛmb = get_cell_coordinates(Λmb)
  xΛmb_n1 = findall(model.grid_topology.vertex_coordinates .== xΛmb[1])
  xΛmb_n2 = findall(model.grid_topology.vertex_coordinates .== xΛmb[2])
  new_entity = num_entities(labels_Ω) + 1
  labels_Ω.d_to_dface_to_entity[1][xΛmb_n1[1]] = new_entity
  labels_Ω.d_to_dface_to_entity[1][xΛmb_n2[1]] = new_entity
  add_tag!(labels_Ω, "mem_bnd", [new_entity])

  # =========================================================================
  # Measures and normals
  # =========================================================================

  degree = 2 * order
  dΩ = Measure(Ω, degree)
  dΓm = Measure(Γm, degree)
  dΓd1 = Measure(Γd1, degree)
  dΓd2 = Measure(Γd2, degree)
  dΓfs = Measure(Γfs, degree)
  dΓin = Measure(Γin, degree)

  # =========================================================================
  # FE spaces (4-field: phi, kappa, eta, q)
  # =========================================================================

  reffe = ReferenceFE(lagrangian, Float64, order)
  V_Ω = TestFESpace(Ω, reffe, conformity=:H1)
  V_Γκ = TestFESpace(Γκ, reffe, conformity=:H1)
  V_Γη = TestFESpace(Γη, reffe, conformity=:H1)
  V_Γq = ConstantFESpace(Ω, vector_type=Vector{Float64},
    field_type=VectorValue{1, Float64})

  U_Ω = TransientTrialFESpace(V_Ω)
  U_Γκ = TransientTrialFESpace(V_Γκ)
  U_Γη = TransientTrialFESpace(V_Γη)
  U_Γq = TransientTrialFESpace(V_Γq)

  î1 = VectorValue(1.0)

  X = TransientMultiFieldFESpace([U_Ω, U_Γκ, U_Γη, U_Γq])
  Y = MultiFieldFESpace([V_Ω, V_Γκ, V_Γη, V_Γq])

  δ_p = DiracDelta(Γ, [Point(15.0, 0.0)])

  # =========================================================================
  # Weak form (free BC, tau=0 branch)
  # =========================================================================

  ∇ₙ(ϕ) = ∇(ϕ) ⋅ VectorValue(0.0, 1.0)

  m(t, (ϕₜₜ, κₜₜ, ηₜₜ, qₜₜ), (w, u, v, ξ)) =
    ∫(mᵨ * v * ηₜₜ)dΓm +
    rM * δ_p(qₜₜ ⋅ ξ)

  c(t, (ϕₜ, κₜ, ηₜ, qₜ), (w, u, v, ξ)) =
    ∫(βₕ*(u + αₕ*w)*ϕₜ - w*κₜ)dΓfs +
    ∫(βₕ*(u + αₕ*w)*ϕₜ - w*κₜ)dΓd1 +
    ∫(βₕ*(u + αₕ*w)*ϕₜ - w*κₜ)dΓd2 +
    ∫(v*ϕₜ - w*ηₜ)dΓm

  a(t, (ϕ, κ, η, q), (w, u, v, ξ)) =
    ∫(∇(w) ⋅ ∇(ϕ))dΩ +
    ∫(βₕ*(u + αₕ*w)*g*κ)dΓfs +
    ∫(βₕ*(u + αₕ*w)*g*κ - μ₂ᵢₙ*κ*w + μ₁ᵢₙ*∇ₙ(ϕ)*(u + αₕ*w))dΓd1 +
    ∫(βₕ*(u + αₕ*w)*g*κ - μ₂ₒᵤₜ*κ*w + μ₁ₒᵤₜ*∇ₙ(ϕ)*(u + αₕ*w))dΓd2 +
    ∫(v*(g*η) + Tᵨ*∇(v) ⋅ ∇(η))dΓm +
    (-rK/ρw) * δ_p(v * (q ⋅ î1 - η)) +
    rK * δ_p(ξ ⋅ q - (ξ ⋅ î1) * η)

  l(t, (w, u, v, ξ)) =
    ∫(w * vᵢₙ(t))dΓin -
    ∫(ηd(t)*w - ∇ₙϕd(t)*(u + αₕ*w))dΓd1

  # =========================================================================
  # Operator and solver
  # =========================================================================

  constant_forms = (true, true, true)
  op = TransientLinearFEOperator((a, c, m), l, X, Y;
    constant_forms=constant_forms)

  @testset "FEM assembly" begin
    @test length(X.spaces) == 4
    @test length(Y.spaces) == 4
  end

  ls = LUSolver()
  ode_solver = GeneralizedAlpha2(ls, Δt, 1.0)

  # Initial conditions
  u0   = interpolate_everywhere([0.0, 0.0, 0.0, VectorValue(0.0)], X(0.0))
  u0t  = interpolate_everywhere([0.0, 0.0, 0.0, VectorValue(0.0)], X(0.0))
  u0tt = interpolate_everywhere([0.0, 0.0, 0.0, VectorValue(0.0)], X(0.0))

  uht = solve(ode_solver, op, t₀, tf, (u0, u0t, u0tt))

  # =========================================================================
  # Time loop — collect diagnostics
  # =========================================================================

  step_count = 0
  ϕ_l2 = Float64[]
  η_l2 = Float64[]
  q_vals = Float64[]

  for (t, uh) in uht
    step_count += 1
    ϕₕ, κₕ, ηₕ, qₕ = uh

    push!(ϕ_l2, sqrt(abs(sum(∫(ϕₕ * ϕₕ)dΩ))))
    push!(η_l2, sqrt(abs(sum(∫(ηₕ * ηₕ)dΓm))))
    push!(q_vals, abs(get_free_dof_values(qₕ)[1]))
  end

  expected_steps = round(Int, tf / Δt)

  # =========================================================================
  # Assertions
  # =========================================================================

  @testset "Time integration" begin
    @test step_count == expected_steps
    @test all(isfinite, ϕ_l2)
    @test all(isfinite, η_l2)
    @test all(isfinite, q_vals)
  end

  @testset "Physical reasonableness" begin
    # Wave has entered the domain
    @test ϕ_l2[end] > 0.0
    # Membrane responds
    @test η_l2[end] > 0.0
    # Resonator responds
    @test any(q_vals .> 0.0)
  end

  @testset "Energy boundedness" begin
    # Solution grows from quiescent start
    @test ϕ_l2[end] > ϕ_l2[1]
    # No blow-up
    @test maximum(ϕ_l2)  < 1e6
    @test maximum(η_l2)  < 1e4
    @test maximum(q_vals) < 1e4
  end

end
