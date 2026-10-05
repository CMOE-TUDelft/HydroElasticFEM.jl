using Test
using Gridap
using Gridap.ODEs
using LinearAlgebra
using WaveSpec
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Geometry as G

const _g = WaveSpec.PhysicalConstants.g

function _airy_ref(ω, kx, ky, A, φ, k, h, x, y, z, t)
  η = ϕ = u = v = w = 0.0
  for i in eachindex(ω)
    θ = kx[i] * x + ky[i] * y - ω[i] * t + φ[i]
    ch = cosh(k[i] * (z + h)) / cosh(k[i] * h)
    chs = cosh(k[i] * (z + h)) / sinh(k[i] * h)
    shs = sinh(k[i] * (z + h)) / sinh(k[i] * h)
    η += A[i] * cos(θ)
    ϕ += A[i] * _g / ω[i] * ch * sin(θ)
    u += A[i] * ω[i] * kx[i] / k[i] * chs * cos(θ)
    v += A[i] * ω[i] * ky[i] / k[i] * chs * cos(θ)
    w += A[i] * ω[i] * shs * sin(θ)
  end
  (; η, ϕ, u, v, w)
end

@testset "IncidentSea fields" begin
  h = 12.0
  ω = [0.8, 0.8, 0.8, 1.3, 1.3, 1.3]
  θd = repeat([-0.3, 0.0, 0.4], outer = 2)
  k = [WaveSpec.AiryWaves.solve_wavenumber(w, h) for w in ω]
  kx, ky = k .* cos.(θd), k .* sin.(θd)
  A = [0.1, 0.2, 0.05, 0.07, 0.03, 0.02]
  φ = [0.3, 1.1, 2.0, 4.0, 5.5, 0.9]
  ramp(t) = min(t / 5, 1.0)
  sea = P.IncidentSea(ω, kx, ky, A, φ, k, h; ramp, z_surface = 2.0)
  @test P.num_frequency_groups(sea) == 2

  η, ϕ, w, vel = P.incident_elevation(sea), P.incident_potential(sea),
                 P.incident_vertical_velocity(sea), P.incident_velocity(sea)
  for (x, t) in ((Point(3.0, -4.0, 1.5), 0.7), (Point(-10.0, 2.5, -6.0), 12.3))
    ref = _airy_ref(ω, kx, ky, A, φ, k, h, x[1], x[2], x[3] - 2.0, t)
    r = ramp(t)
    @test η(x, t) ≈ r * ref.η
    @test ϕ(x, t) ≈ r * ref.ϕ
    @test w(x, t) ≈ r * ref.w
    @test vel(x, t) ≈ r * VectorValue(ref.u, ref.v, ref.w)

    for f in (η, ϕ, w, vel)
      rec = sum(1:2) do g_
        c = P.with_sea_mode(() -> f(x, t), sea, :cos, g_)
        s = P.with_sea_mode(() -> f(x, t), sea, :sin, g_)
        cos(sea.ω_groups[g_] * t) * c + sin(sea.ω_groups[g_] * t) * s
      end
      @test r * rec ≈ f(x, t)
      @test iszero(P.with_sea_mode(() -> f(x, t), sea, :off))
    end
  end
  @test sea.mode[] == (:time, 0)

  sea2 = P.IncidentSea(ω[[1, 4]], k[[1, 4]], zeros(2), A[[1, 4]], φ[[1, 4]], k[[1, 4]], h; dim = 2)
  ref2 = _airy_ref(ω[[1, 4]], k[[1, 4]], zeros(2), A[[1, 4]], φ[[1, 4]], k[[1, 4]], h, 5.0, 0.0, -3.0, 2.0)
  @test P.incident_velocity(sea2)(Point(5.0, -3.0), 2.0) ≈ VectorValue(ref2.u, ref2.w)

  spec = WaveSpec.ContinuousSpectrums.RegularWave(0.4, 6.0)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess = false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
  ωs = [2π / 6.0]
  st = WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, ωs,
    [WaveSpec.AiryWaves.solve_wavenumber(ωs[1], h)], [0.0], h, 7)
  seaS = P.IncidentSea(st)
  @test seaS.A ≈ vec(permutedims(WaveSpec.AiryWaves.get_amplitudes(st)))
  @test seaS.φ ≈ vec(permutedims(WaveSpec.AiryWaves.get_random_phases(st)))
  @test seaS.kx ≈ st.k && all(iszero, seaS.ky)
end

@testset "Spectral forcing of an IncidentSea-driven problem" begin
  H, L, Ld, Lm = 10.0, 60.0, 10.0, 10.0
  ωc = [0.9, 0.9, 1.4]
  kc = [WaveSpec.AiryWaves.solve_wavenumber(w, H) for w in ωc]
  sea = P.IncidentSea(ωc, kc, zeros(3), [0.05, 0.03, 0.02], [0.2, 2.5, 1.0], kc, H;
                      dim = 2, z_surface = H, ramp = t -> min(t / 2, 1.0))
  tank = G.TankDomain(L = L, H = H, nx = 30, ny = 4,
    structure_domains = [G.StructureDomain(L = Lm, x₀ = [25.0, H])],
    damping_zones = [G.DampingZone(L = Ld, x₀ = [0.0, H], domain_symbol = :Γd_in),
                     G.DampingZone(L = Ld, x₀ = [L - Ld, H], domain_symbol = :Γd_out)])
  fe = PH.FESpaceConfig(vector_type = Vector{Float64})
  μ₁(x) = 2.0 * (1 - x[1] / Ld)
  bcs = [P.PrescribedInletPotentialBC(domain = :dΓin, quantity = :velocity,
                                      forcing = P.incident_velocity(sea)),
         P.DampingZoneBC(domain = :dΓd_1, μ₁ = μ₁, μ₂ = (x -> 1.0),
                         η_in = P.incident_elevation(sea), vz_in = P.incident_vertical_velocity(sea)),
         P.DampingZoneBC(domain = :dΓd_2, μ₁ = (x -> 2.0), μ₂ = (x -> 1.0),
                         η_in = (x -> 0.0), vz_in = (x -> 0.0))]
  physics = P.PhysicsParameters[P.PotentialFlow(fe = fe, boundary_conditions = bcs),
                                P.FreeSurface(fe = fe),
                                P.Membrane(L = Lm, mᵨ = 0.9, Tᵨ = 98.1, fe = fe)]
  Δt = 0.1
  tc = PH.TimeConfig(Δt = Δt, tf = 1.0, αₕ = 0.5 / (0.25 * Δt) / 9.81,
                     u0 = zeros(3), u0t = zeros(3), u0tt = zeros(3))
  prob = SM.build_problem(tank, physics, PH.TimeDomainConfig(tf = 1.0); tconfig = tc)
  op = SM.get_fe_operator(prob)
  @test op isa FO.SpectralForcingTFEOperator
  @test op.forcing isa FO.SpectralForcing
  @test length(op.forcing.ω) == 2

  X, Y, fmap = SM.get_trial_fe_space(prob), SM.get_test_fe_space(prob), SM.get_field_map(prob)
  ctx = SM.get_assembly_context(prob)
  op_std = FO.build_time_fe_operator(physics, ctx, fmap, X, Y;
                                     rhs_fn = (t, y) -> zeros(length(fmap)))
  V = Gridap.FESpaces.get_test(op_std)
  res = Gridap.ODEs.get_res(op_std)
  X0 = X(0.0)
  @testset "forcing equals the assembled right-hand side" begin
    for t in (0.0, 0.35, 1.0, 7.3)
      F = assemble_vector(v -> res(t, zero(X0), v), V)
      r = zeros(length(F))
      FO._subtract_forcing!(r, op.forcing, t)
      @test norm(-r - F) <= 1e-11 * max(norm(F), 1.0)
    end
  end

  @testset "time integration matches the stock operator" begin
    u0 = (zero(X0), zero(X0), zero(X0))
    solver = GeneralizedAlpha2(LUSolver(), Δt, 0.9)
    s_alg = [copy(get_free_dof_values(u)) for (_, u) in solve(solver, op, 0.0, 1.0, u0)]
    s_std = [copy(get_free_dof_values(u)) for (_, u) in solve(solver, op_std, 0.0, 1.0, u0)]
    @test maximum(abs, s_std[end]) > 0
    for (a, b) in zip(s_alg, s_std)
      @test norm(a - b) <= 1e-10 * norm(b)
    end
  end

  @testset "falls back for other time-dependent inputs" begin
    bcs2 = [bcs[1], bcs[2],
            P.PrescribedInletPotentialBC(domain = :dΓbot, quantity = :normal_gradient,
                                         forcing = (x, t) -> 0.01 * sin(t))]
    phys2 = P.PhysicsParameters[P.PotentialFlow(fe = fe, boundary_conditions = bcs2), physics[2], physics[3]]
    @test FO._incident_seas(phys2) === nothing
    op2 = FO.build_time_fe_operator(phys2, ctx, fmap, X, Y)
    @test op2 isa Gridap.ODEs.TransientFEOperator
    @test !(op2 isa FO.SpectralForcingTFEOperator)
    op3 = FO.build_time_fe_operator(physics, ctx, fmap, X, Y;
                                    rhs_fn = (t, y) -> zeros(length(fmap)))
    @test !(op3 isa FO.SpectralForcingTFEOperator)
  end
end
