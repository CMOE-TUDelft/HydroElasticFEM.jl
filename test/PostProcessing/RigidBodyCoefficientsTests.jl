using Test
using Gridap
using LinearAlgebra
using WaveSpec
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO
import HydroElasticFEM.PostProcessing as PP

# --------------------------------------------------------------------------
# 2D zero-draft rigid pontoon (heave, pitch) in a tank with absorbing ends.
# Free surface at z = 0; incident wave from the left (inlet at x = 0).
# --------------------------------------------------------------------------

const _HC = (H0 = 10.0, LΩ = 80.0, nx = 80, ny = 4, xb0 = 35.0, Lb = 10.0,
             ω = 1.5, η0 = 0.1, order = 2, ρw = 1025.0, g = 9.81, mᵨ = 0.9)

function _hc_sea()
  p = _HC
  spec = WaveSpec.ContinuousSpectrums.RegularWave(2p.η0, 2π / p.ω)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess = false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
  k = WaveSpec.AiryWaves.solve_wavenumber(p.ω, p.H0)
  sea = WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, [p.ω], [k], [0.0], p.H0, 1)
  # Complex amplitude of the incident Airy potential (time dependence
  # e^{−iωt}, η = iωϕ/g at z = 0):  ϕ̂_I = −i g η0/ω cosh k(z+H)/cosh kH e^{ikx}.
  # (Real WaveSpec snapshots used as complex amplitudes would only feed half
  # of the wave into the tank.)  On the inlet, n = −e_x and ∂ₙϕ̂_I = −ik ϕ̂_I.
  ϕI(x) = -im * p.g * p.η0 / p.ω * cosh(k * (x[2] + p.H0)) / cosh(k * p.H0) * exp(im * k * x[1])
  f_in(x) = -2im * k * ϕI(x)
  return sea, f_in, k
end

const _HC_SEA, _HC_FIN, _HC_K = _hc_sea()
const _HC_FE = PH.FESpaceConfig(order = _HC.order, vector_type = Vector{ComplexF64})

_hc_tank() = G.TankDomain(L = _HC.LΩ, H = _HC.H0, nx = _HC.nx, ny = _HC.ny,
  map = x -> VectorValue(x[1], x[2] - _HC.H0),
  structure_domains = [G.StructureDomain(L = _HC.Lb, x₀ = [_HC.xb0, 0.0], domain_symbol = :Γ_b)])

function _hc_fluid(; incident::Bool)
  bcs = Any[P.RadiationBC(domain = :dΓin), P.RadiationBC(domain = :dΓout)]
  incident && push!(bcs, P.PrescribedInletPotentialBC(domain = :dΓin, forcing = _HC_FIN, quantity = :traction))
  pf = P.PotentialFlow(ρw = _HC.ρw, g = _HC.g, sea_state = _HC_SEA, boundary_conditions = [bcs...],
                       fe = _HC_FE, space_domain_symbol = :Ω)
  fs = P.FreeSurface(ρw = _HC.ρw, g = _HC.g, βₕ = 0.5, fe = _HC_FE, space_domain_symbol = :Γκ)
  P.PhysicsParameters[pf, fs]
end

function _hc_body(; B_ext = nothing)
  m = _HC.ρw * _HC.mᵨ * _HC.Lb
  P.RigidBody(mass = m, inertia = m * _HC.Lb^2 / 12, x_G = [_HC.xb0 + _HC.Lb / 2, 0.0],
              modes = (:heave, :pitch), B_ext = B_ext, ρw = _HC.ρw, g = _HC.g, fe = _HC_FE)
end

const _HC_LW = 2π / _HC_K
const _HC_XS_IN = PP.probe_positions(_HC.xb0, :upwave, _HC_LW, _HC.H0)
const _HC_XS_OUT = PP.probe_positions(_HC.xb0 + _HC.Lb, :downwave, _HC_LW, _HC.H0)
_hc_fit(κh, xs) = PP.fit_wave_components(xs, PP.sample_probe_line(κh, xs, 0.0), _HC_K)

@testset "PrescribedRigidMotion" begin
  body = _hc_body()
  pm = P.PrescribedRigidMotion(body; mode = :pitch)
  @test pm.amplitudes == [0.0, 1.0]
  @test P.variable_symbols(pm) == () && !P.has_rhs_form(pm)
  @test P.has_rhs_form(P.PotentialFlow(), pm)
  @test_throws ArgumentError P.PrescribedRigidMotion(body)
  @test_throws ArgumentError P.PrescribedRigidMotion(body; mode = :heave, amplitudes = [1.0, 0.0])
  @test_throws ArgumentError P.PrescribedRigidMotion(body; mode = :surge)
  @test_throws ArgumentError P.PrescribedRigidMotion(body; amplitudes = [1.0])
  @test_nowarn P.print_parameters(pm)

  # rhs = −iω ∫ w (ξ̂·N): for ξ̂ = (a, b) on the top boundary, N = (1, −r_x)
  tank = _hc_tank()
  tr = G.build_triangulations(tank, G.build_model(tank))
  physics = P.PhysicsParameters[_hc_fluid(incident = false)...]
  dom = G.get_integration_domains(tr; degree = SM.get_integration_degrees(tr, physics))
  X, Y, fmap = SM.build_fe_spaces(physics, tr, PH.FreqDomainConfig(ω = _HC.ω))
  ctx = SM.build_frequency_context(dom, physics, PH.FreqDomainConfig(ω = _HC.ω))
  a, b = 0.3 + 0.1im, -0.02 + 0.05im
  pm2 = P.PrescribedRigidMotion(body; amplitudes = [a, b])
  l = assemble_vector(y -> P.rhs(physics[1], pm2, ctx, nothing, FO.FieldMap(y, fmap)), Y)
  xG = body.x_G[1]
  Vϕ = Y.spaces[fmap[:ϕ]]
  ref = assemble_vector(w -> ∫(-im * _HC.ω * w * (a + b * ((x -> -(x[1] - xG)) ∘ get_physical_coordinate(tr[:Γ_b]))))dom[:dΓ_b], Vϕ)
  @test l[1:num_free_dofs(Vϕ)] ≈ ref
  @test norm(l[num_free_dofs(Vϕ)+1:end]) < 1e-14 * norm(ref)
end

@testset "Rigid-body hydrodynamic coefficients (2D zero-draft pontoon)" begin
  p = _HC
  tank = _hc_tank()
  body = _hc_body()
  cfg = PH.FreqDomainConfig(ω = p.ω)
  Cg = PP.group_velocity(_HC_K, p.H0; g = p.g)

  # ── Radiation: A, B and radiated wave energy ──
  A, B = PP.radiation_coefficients(tank, _hc_fluid(incident = false), body, cfg)
  @test size(A) == size(B) == (2, 2)
  @test A ≈ transpose(A) rtol = 1e-10        # reciprocity (exact for the discrete operator)
  @test B ≈ transpose(B) rtol = 1e-10
  @test all(diag(A) .> 0)
  @test all(eigvals(Symmetric((B + B') / 2)) .> 0)

  # Radiated power ½ω²B_jj (unit amplitude) leaves as outgoing waves on both
  # sides: ½ρg C_g (|a_L|² + |a_R|²).  With e^{−iωt}, e^{+ikx} travels to the
  # right, so the outgoing amplitudes are fit.B on the left, fit.A on the right.
  for (j, mode) in enumerate(body.modes)
    _, res = PP.solve_radiation(tank, _hc_fluid(incident = false), body, cfg; mode = mode)
    κh = res.solution[2]
    fl, fr = _hc_fit(κh, _HC_XS_IN), _hc_fit(κh, _HC_XS_OUT)
    @test abs(fl.A) < 0.01 * abs(fl.B) && abs(fr.B) < 0.01 * abs(fr.A)  # outgoing only (≈ 7e-4)
    P_rad = 0.5 * p.ρw * p.g * Cg * (abs2(fl.B) + abs2(fr.A))
    @test 0.5 * p.ω^2 * B[j, j] ≈ P_rad rtol = 5e-3       # observed ≈ 4e-4
  end

  # ── Diffraction: excitation force on the fixed body ──
  prob_d = SM.build_problem(tank, _hc_fluid(incident = true), cfg)
  res_d = SM.simulate(prob_d)
  η_inc = _hc_fit(res_d.solution[2], _HC_XS_IN).A
  @test abs(η_inc) ≈ p.η0 rtol = 0.05
  X = PP.excitation_force(prob_d, res_d, body)
  # Haskind (2D, body symmetric about x_G, so |X(β=0)| = |X(β=π)|):
  #   B_jj = |X_j|² / (2ρ g C_g), X per unit incident amplitude
  for j in 1:2
    @test B[j, j] ≈ abs2(X[j] / η_inc) / (2 * p.ρw * p.g * Cg) rtol = 5e-3   # observed ≈ 2.5e-4
  end

  # ── Coupled solve = response from (A, B, X), and energy balance ──
  for Bpto in (0.0, B[1, 1])
    Bext = zeros(3, 3); Bext[2, 2] = Bpto              # heave PTO (full-mode index 2)
    bodyc = _hc_body(B_ext = Bext)
    prob_c = SM.build_problem(tank, P.PhysicsParameters[_hc_fluid(incident = true)..., bodyc], cfg)
    res_c = SM.simulate(prob_c)
    ξ = collect(SM.body_motion(res_c, bodyc))
    ξ_ref = PP.rigid_body_response(bodyc, SM.get_integration_domains(prob_c), p.ω, A, B, X)
    @test ξ ≈ ξ_ref rtol = 1e-8

    κh = res_c.solution[2]
    rt = PP.reflection_transmission_coefficients(κh, _HC_XS_IN, _HC_XS_OUT, _HC_K; y = 0.0)
    P_pto = PP.rigid_body_dissipated_power(ξ, p.ω, Bext[2:3, 2:3])
    A_pto = PP.absorption_from_dissipated_power(P_pto, abs(rt.incident_fit.A), _HC_K, p.H0; ρw = p.ρw, g = p.g)
    @test abs2(rt.R) + abs2(rt.T) + A_pto ≈ 1.0 atol = 2e-3   # observed ≤ 2e-4
    # At this ω the pontoon is at heave resonance (−ω²(M + A₃₃) + C₃₃ ≈ 0), so
    # the matched PTO B_pto = B₃₃ reaches the maximum absorption of a
    # symmetric 2D body in one mode, 1/2 (Evans 1976; Mei 1976).
    Bpto > 0 && @test A_pto ≈ 0.5 rtol = 0.01
  end
end
