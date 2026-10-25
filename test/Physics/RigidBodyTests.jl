using Test
using Gridap
using LinearAlgebra
using WaveSpec
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO

_rb_cross(a, b) = [a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]]
_rb_skew(r) = [0.0 -r[3] r[2]; r[3] 0.0 -r[1]; -r[2] r[1] 0.0]

@testset "RigidBody — constructor and validation" begin
  b = P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], modes = (:heave, :pitch))
  @test b.ambient_dim == 2 && P.num_modes(b) == 2 && b.x_ref == b.x_G
  @test P.variable_symbol(b) == :ξ && b.space_domain_symbol == :Γ_b
  @test !P.has_damping_form(b) && P.has_mass_form(b) && P.has_stiffness_form(b)
  @test P.has_damping_form(P.PotentialFlow(), b)
  @test_nowarn P.print_parameters(b)
  @test P.rigid_body_modes(3) == (:surge, :sway, :heave, :roll, :pitch, :yaw)

  # zero-draft surge has no restoring and no fluid coupling
  @test_throws ArgumentError P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0])
  K = zeros(3, 3); K[1, 1] = 5.0
  @test P.num_modes(P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], K_ext = K)) == 3
  @test P.num_modes(P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], hydrostatics = :none)) == 3
  @test_throws ArgumentError P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], modes = (:roll,))
  @test_throws ArgumentError P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], modes = (:heave, :heave))
  @test_throws ArgumentError P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [0.0, 0.0, 0.0], modes = (:heave,))
  @test_throws ArgumentError P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], modes = (:heave,),
                                         K_ext = zeros(2, 2))
  @test_throws ArgumentError P.RigidBody(mass = -1.0, inertia = 2.0, x_G = [1.0, 0.0], modes = (:heave,))
  @test_throws ArgumentError P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], modes = (:heave,),
                                         hydrostatics = :from_somewhere)
end

@testset "RigidBody — kinematics" begin
  # 2D: u = (ξ₁ + α r_z, ξ₃ − α r_x)
  b2 = P.RigidBody(mass = 1.0, inertia = 1.0, x_G = [1.0, 2.0], hydrostatics = :none)
  @test P.rigid_body_displacement(b2, [0.0, 0.0, 0.1], [3.0, 2.5]) ≈ [0.1 * 0.5, -0.1 * 2.0]
  @test P.rigid_body_displacement(b2, [1.0, 2.0, 0.0], [7.0, -3.0]) ≈ [1.0, 2.0]

  # ξ·N(x, n) = u(x)·n and ξ·U_z(x) = u_z(x), 2D and 3D, all modes and subsets
  b3 = P.RigidBody(mass = 1.0, inertia = Matrix(1.0I, 3, 3), x_G = [1.0, -2.0, 0.5],
                   x_ref = [0.3, 0.2, -0.1], hydrostatics = :none)
  b2s = P.RigidBody(mass = 1.0, inertia = 1.0, x_G = [1.0, 2.0], modes = (:pitch, :heave))
  b3s = P.RigidBody(mass = 1.0, inertia = Matrix(1.0I, 3, 3), x_G = [1.0, -2.0, 0.5],
                    modes = (:roll, :heave, :pitch))
  for b in (b2, b3, b2s, b3s)
    D = b.ambient_dim
    N = P.generalized_normal(b)
    U = P.vertical_shape(b)
    for k in 1:5
      ξ = randn(P.num_modes(b))
      x = randn(D)
      n = normalize(randn(D))
      u = P.rigid_body_displacement(b, ξ, x)
      @test dot(ξ, collect(N(VectorValue(x...), VectorValue(n...)).data)) ≈ dot(u, n)
      @test dot(ξ, collect(U(VectorValue(x...)).data)) ≈ u[end]
    end
  end
  # 3D rotation matches α × r
  α = [0.1, -0.2, 0.3]; x = [2.0, 1.0, -1.0]
  @test P.rigid_body_displacement(b3, [0, 0, 0, α...], x) ≈ _rb_cross(α, x .- b3.x_ref)
end

@testset "RigidBody — mass matrix from point masses" begin
  # M = Σ m_k G_kᵀ G_k with u_k = G_k ξ, for a cloud of point masses
  pts3 = [randn(3) for _ in 1:6]; ms = rand(6) .+ 0.5
  m = sum(ms); xG = sum(ms .* pts3) / m
  IG = sum(mk * (sum(abs2, p - xG) * I(3) - (p - xG) * (p - xG)') for (mk, p) in zip(ms, pts3))
  xr = [0.4, -0.3, 0.2]
  M3 = sum(mk * (G = [I(3) -_rb_skew(p - xr)]; G' * G) for (mk, p) in zip(ms, pts3))
  b3 = P.RigidBody(mass = m, inertia = IG, x_G = xG, x_ref = xr, hydrostatics = :none)
  @test P.rigid_body_mass_matrix(b3) ≈ M3

  pts2 = [randn(2) for _ in 1:5]; ms2 = rand(5) .+ 0.5
  m2 = sum(ms2); xG2 = sum(ms2 .* pts2) / m2
  Iyy = sum(mk * sum(abs2, p - xG2) for (mk, p) in zip(ms2, pts2))
  xr2 = [0.2, 0.7]
  M2 = sum(mk * (r = p - xr2; G = [1.0 0.0 r[2]; 0.0 1.0 -r[1]]; G' * G) for (mk, p) in zip(ms2, pts2))
  b2 = P.RigidBody(mass = m2, inertia = Iyy, x_G = xG2, x_ref = xr2, hydrostatics = :none)
  @test P.rigid_body_mass_matrix(b2) ≈ M2
end

@testset "RigidBody — assembled blocks" begin
  ρw, g, L = 1025.0, 9.81, 2.0
  tank = G.TankDomain(L = 8.0, H = 1.0, nx = 16, ny = 2,
    structure_domains = [G.StructureDomain(L = L, x₀ = [3.0, 1.0], domain_symbol = :Γ_b)])
  tr = G.build_triangulations(tank, G.build_model(tank))
  fe = PH.FESpaceConfig(vector_type = Vector{ComplexF64})
  Kx = [4.0e3 1.0e2 -3.0e2; 2.0e2 7.0e3 5.0e1; 0.0 -1.0e2 9.0e3]       # not symmetric
  Bx = [1.0e2 0.0 0.0; 0.0 2.0e2 3.0e1; 0.0 3.0e1 5.0e2]
  body = P.RigidBody(mass = 2.0e3, inertia = 7.0e2, x_G = [4.0, 1.2], x_ref = [3.5, 1.0],
                     K_ext = Kx, B_ext = Bx, ρw = ρw, g = g, fe = fe)
  physics = P.PhysicsParameters[P.PotentialFlow(fe = fe), P.FreeSurface(fe = fe), body]
  dom = G.get_integration_domains(tr; degree = SM.get_integration_degrees(tr, physics))
  X, Y, fmap = SM.build_fe_spaces(physics, tr, PH.FreqDomainConfig(ω = 1.0))
  @test fmap[:ξ] == 3 && num_free_dofs(Y.spaces[3]) == 3

  offs = cumsum([0; [num_free_dofs(V) for V in Y.spaces]])
  blk(A, i, j) = Matrix(A[offs[i]+1:offs[i+1], offs[j]+1:offs[j+1]])
  asm(f) = assemble_matrix((x, y) -> f(FO.FieldMap(x, fmap), FO.FieldMap(y, fmap)), X, Y)
  iξ = fmap[:ξ]

  @test blk(asm((x, y) -> P.mass(body, dom, x, y)), iξ, iξ) ≈ P.rigid_body_mass_matrix(body) / ρw
  @test blk(asm((x, y) -> P.damping(body, dom, x, y)), iξ, iξ) ≈ Bx / ρw
  # zero-draft restoring about x_ref: C_heave = gL, C_heave,pitch = −g∫r_x,
  # C_pitch = g∫r_x², with r_x ∈ (3, 5) − 3.5
  C = zeros(3, 3)
  C[2, 2] = g * L
  C[2, 3] = C[3, 2] = -g * ((5.0 - 3.5)^2 - (3.0 - 3.5)^2) / 2
  C[3, 3] = g * ((5.0 - 3.5)^3 - (3.0 - 3.5)^3) / 3
  @test blk(asm((x, y) -> P.stiffness(body, dom, x, y)), iξ, iξ) ≈ Kx / ρw + C

  # Mode subset
  bs = P.RigidBody(mass = 2.0e3, inertia = 7.0e2, x_G = [4.0, 1.2], x_ref = [3.5, 1.0],
                   modes = (:pitch, :heave), ρw = ρw, g = g, fe = fe)
  physics_s = P.PhysicsParameters[P.PotentialFlow(fe = fe), P.FreeSurface(fe = fe), bs]
  Xs, Ys, fmaps = SM.build_fe_spaces(physics_s, tr, PH.FreqDomainConfig(ω = 1.0))
  offs_s = cumsum([0; [num_free_dofs(V) for V in Ys.spaces]])
  Ms = assemble_matrix((x, y) -> P.mass(bs, dom, FO.FieldMap(x, fmaps), FO.FieldMap(y, fmaps)), Xs, Ys)
  M = P.rigid_body_mass_matrix(bs)
  @test Matrix(Ms[offs_s[3]+1:end, offs_s[3]+1:end]) ≈ M[[3, 2], [3, 2]] / ρw

  # Fluid coupling: body row = −(fluid column)ᵀ, heave column sums to −L
  Cc = asm((x, y) -> P.damping(physics[1], body, dom, x, y))
  iϕ = fmap[:ϕ]
  @test blk(Cc, iξ, iϕ) ≈ -transpose(blk(Cc, iϕ, iξ))
  @test sum(blk(Cc, iϕ, iξ)[:, 2]) ≈ -L
  @test norm(blk(Cc, iϕ, iξ)[:, 1]) < 1e-12
end

# --------------------------------------------------------------------------
# Rigid limit: a free-free Euler-Bernoulli beam with very large EI responds
# like a zero-draft rigid body in heave and pitch (same mass and inertia).
# --------------------------------------------------------------------------

function _rb_incident_wave(; H0, ω, η0)
  spec = WaveSpec.ContinuousSpectrums.RegularWave(2η0, 2π / ω)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess = false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
  k = WaveSpec.AiryWaves.solve_wavenumber(ω, H0)
  sea = WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, [ω], [k], [0.0], H0, 1)
  wave(x) = WaveSpec.AiryWaves.generate_sea(sea, [x[1]], [0.0], [x[2]], [0.0], vars = [:ϕ, :u, :w])
  ϕin(x) = wave(x)[:ϕ][1]
  vin(x) = VectorValue(wave(x)[:u][1], wave(x)[:w][1])
  f_in(x) = (vin(x) ⋅ VectorValue(-1.0, 0.0)) - im * k * ϕin(x)
  return sea, f_in
end

function _rb_solve(structure; H0, LΩ, nx, ny, xb0, Lb, ω, η0, order)
  tank = G.TankDomain(L = LΩ, H = H0, nx = nx, ny = ny,
    map = x -> VectorValue(x[1], x[2] - H0),
    structure_domains = [G.StructureDomain(L = Lb, x₀ = [xb0, 0.0], domain_symbol = :Γ_b)])
  sea, f_in = _rb_incident_wave(; H0, ω, η0)
  fe = PH.FESpaceConfig(order = order, vector_type = Vector{ComplexF64})
  pf = P.PotentialFlow(ρw = 1025.0, g = 9.81, sea_state = sea,
    boundary_conditions = [P.RadiationBC(domain = :dΓin), P.RadiationBC(domain = :dΓout),
                           P.PrescribedInletPotentialBC(domain = :dΓin, forcing = f_in, quantity = :traction)],
    fe = fe, space_domain_symbol = :Ω)
  fs = P.FreeSurface(ρw = 1025.0, g = 9.81, βₕ = 0.5, fe = fe, space_domain_symbol = :Γκ)
  physics = P.PhysicsParameters[pf, fs, structure]
  SM.simulate(SM.build_problem(tank, physics, SM.FreqDomainConfig(ω = ω)))
end

@testset "RigidBody — rigid limit of a stiff Euler-Bernoulli beam" begin
  ρw, g = 1025.0, 9.81
  prm = (H0 = 10.0, LΩ = 60.0, nx = 60, ny = 4, xb0 = 25.0, Lb = 10.0, ω = 1.5, η0 = 0.1, order = 2)
  mᵨ = 0.9
  fe = PH.FESpaceConfig(order = prm.order, vector_type = Vector{ComplexF64})
  beam = P.EulerBernoulliBeam(L = prm.Lb, mᵨ = mᵨ, EIᵨ = 1.0e9, g = g, fe = fe)
  m = ρw * mᵨ * prm.Lb
  xG = [prm.xb0 + prm.Lb / 2, 0.0]
  body = P.RigidBody(mass = m, inertia = m * prm.Lb^2 / 12, x_G = xG, modes = (:heave, :pitch),
                     ρw = ρw, g = g, fe = fe)

  res_b = _rb_solve(beam; prm...)
  res_r = _rb_solve(body; prm...)
  ϕb, _, ηb = res_b.solution
  ϕr = res_r.solution[1]

  q = SM.body_motion(res_r, body)
  @test keys(q) == (:heave, :pitch)
  @test abs(q.heave) > 1e-3 * prm.η0     # the body actually moves

  xs = range(prm.xb0 + 0.5, prm.xb0 + prm.Lb - 0.5, length = 7)
  η_beam = ηb(Point.(xs, 0.0))
  η_body = [q.heave - q.pitch * (x - xG[1]) for x in xs]
  # the difference scales like 1/EI (≈ 2e-2 at EIᵨ = 1e3, 2e-4 at 1e5) down to
  # a discretization floor (≈ 8e-5 here)
  @test maximum(abs, η_beam .- η_body) < 1e-3 * maximum(abs, η_body)
  # …and the comparison discriminates: a flexible beam does not match
  soft = P.EulerBernoulliBeam(L = prm.Lb, mᵨ = mᵨ, EIᵨ = 1.0e3, g = g, fe = fe)
  η_soft = _rb_solve(soft; prm...).solution[3](Point.(xs, 0.0))
  @test maximum(abs, η_soft .- η_body) > 1e-2 * maximum(abs, η_body)

  probes = [Point(x, z) for x in (5.0, 20.0, 45.0, 55.0) for z in (0.0, -5.0)]
  @test maximum(abs, ϕb(probes) .- ϕr(probes)) < 1e-3 * maximum(abs, ϕr(probes))
end

# --------------------------------------------------------------------------
# hydrostatics = :from_mesh on finite-draft boxes (BodyCutout) against the
# closed-form restoring matrix of a box (Newman 1977, §6.16)
# --------------------------------------------------------------------------

@testset "RigidBody — hydrostatics from the mesh" begin
  ρw, g = 1025.0, 9.81
  ρg = ρw * g
  @test P.RigidBody(mass = 1.0, inertia = 2.0, x_G = [1.0, 0.0], hydrostatics = :from_mesh).z_surface == 0.0

  # 2D box x ∈ [3, 5], draft 0.5 (submerged area V = 1), free surface z = 0
  B, d = 2.0, 0.5
  tank = G.TankDomain(L = 8.0, H = 2.0, nx = 32, ny = 8, map = x -> VectorValue(x[1], x[2] - 2.0),
                      body_cutouts = [G.BodyCutout(L = B, d = d, x₀ = [3.0, 0.0])])
  tr = G.build_triangulations(tank, G.build_model(tank))
  dom = G.get_integration_domains(tr)
  fe = PH.FESpaceConfig(vector_type = Vector{ComplexF64})
  xG, xr = [4.0, -0.1], [4.5, 0.2]              # reference point off the body centre
  body = P.RigidBody(mass = ρw * B * d, inertia = 1.0, x_G = xG, x_ref = xr,
                     hydrostatics = :from_mesh, ρw = ρw, g = g, fe = fe)
  rx = 4.0 - xr[1]                             # waterplane centre relative to x_ref
  C = zeros(3, 3)
  C[2, 2] = ρg * B
  C[2, 3] = C[3, 2] = -ρg * B * rx
  C[3, 3] = ρg * (B^3 / 12 + B * rx^2) + ρg * B * d * (-d / 2 - xr[2]) - ρg * B * d * (xG[2] - xr[2])
  Ch = @test_logs P.hydrostatic_matrix(body, dom)   # in equilibrium: no warning
  @test Ch ≈ C rtol = 1e-12
  @test P.restoring_matrix(body, dom) ≈ C rtol = 1e-12

  # assembled stiffness block = (C + K_ext)/ρw; surge is coupled to the fluid
  # through the walls, so it needs no restraint
  Kx = zeros(3, 3); Kx[1, 1] = 1.0e3
  bk = P.RigidBody(mass = ρw * B * d, inertia = 1.0, x_G = xG, x_ref = xr, K_ext = Kx,
                   hydrostatics = :from_mesh, ρw = ρw, g = g, fe = fe)
  physics = P.PhysicsParameters[P.PotentialFlow(fe = fe), P.FreeSurface(fe = fe), bk]
  X, Y, fmap = SM.build_fe_spaces(physics, tr, PH.FreqDomainConfig(ω = 1.0))
  K = assemble_matrix((x, y) -> P.stiffness(bk, dom, FO.FieldMap(x, fmap), FO.FieldMap(y, fmap)), X, Y)
  n0 = sum(num_free_dofs(Y.spaces[i]) for i in 1:fmap[:ξ]-1)
  @test Matrix(K[n0+1:end, n0+1:end]) ≈ (C + Kx) / ρw rtol = 1e-12
  Cc = assemble_matrix((x, y) -> P.damping(physics[1], bk, dom, FO.FieldMap(x, fmap), FO.FieldMap(y, fmap)), X, Y)
  @test norm(Cc[1:num_free_dofs(Y.spaces[1]), n0+1]) > 0   # surge column of the fluid coupling

  # zero draft: :from_mesh reduces to :zero_draft when z_G = z_ref
  tz = G.TankDomain(L = 8.0, H = 1.0, nx = 16, ny = 2,
    structure_domains = [G.StructureDomain(L = 2.0, x₀ = [3.0, 1.0], domain_symbol = :Γ_b)])
  domz = G.get_integration_domains(G.build_triangulations(tz, G.build_model(tz)); degree = Dict(:dΩ => 4, :Γ_b => 4))
  kw = (mass = 1.0, inertia = 1.0, x_G = [4.0, 1.0], x_ref = [3.5, 1.0], modes = (:heave, :pitch), ρw = ρw, g = g)
  @test P.restoring_matrix(P.RigidBody(; kw..., hydrostatics = :from_mesh, z_surface = 1.0), domz) ≈
        P.restoring_matrix(P.RigidBody(; kw...), domz) rtol = 1e-12

  # 3D box x ∈ [3, 5], y ∈ [1, 3], draft 0.5 (V = 2)
  t3 = G.TankDomain(L = 8.0, W = 4.0, H = 2.0, nx = 16, ny = 8, nz = 4,
                    map = x -> VectorValue(x[1], x[2], x[3] - 2.0),
                    body_cutouts = [G.BodyCutout(L = 2.0, W = 2.0, d = 0.5, x₀ = [3.0, 1.0, 0.0])])
  dom3 = G.get_integration_domains(G.build_triangulations(t3, G.build_model(t3)))
  xG3, xr3, V = [4.0, 2.0, -0.1], [4.5, 1.75, 0.2], 2.0
  b3 = P.RigidBody(mass = ρw * V, inertia = Matrix(1.0I, 3, 3), x_G = xG3, x_ref = xr3,
                   hydrostatics = :from_mesh, ρw = ρw, g = g)
  cx, cy = 4.0 - xr3[1], 2.0 - xr3[2]          # waterplane centre, 2 × 2 square
  A, Ixx = 4.0, 2.0^4 / 12
  rB = ρg * V * (-0.25 - xr3[3]) - ρg * V * (xG3[3] - xr3[3])
  C3 = zeros(6, 6)
  C3[3, 3] = ρg * A
  C3[3, 4] = C3[4, 3] = ρg * A * cy
  C3[3, 5] = C3[5, 3] = -ρg * A * cx
  C3[4, 5] = C3[5, 4] = -ρg * A * cx * cy
  C3[4, 4] = ρg * (Ixx + A * cy^2) + rB
  C3[5, 5] = ρg * (Ixx + A * cx^2) + rB
  C3[4, 6] = -ρg * V * (4.0 - xr3[1]) + ρg * V * (xG3[1] - xr3[1])   # = 0: x_B = x_G
  C3[5, 6] = -ρg * V * (2.0 - xr3[2]) + ρg * V * (xG3[2] - xr3[2])
  @test P.hydrostatic_matrix(b3, dom3) ≈ C3 rtol = 1e-12
end
