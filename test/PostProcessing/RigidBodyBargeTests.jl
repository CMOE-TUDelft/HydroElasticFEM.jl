using Test
using Gridap
using LinearAlgebra
using WaveSpec
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.PostProcessing as PP

include("BoxEigenfunctionReference.jl")

# --------------------------------------------------------------------------
# 2D rectangular barge (breadth 10, draft 2) in water of depth 10: surge,
# heave and pitch about the waterline centre.  Literature benchmark: the
# matched eigenfunction expansion of Zheng, You & Shen (2004), reimplemented
# independently in BoxEigenfunctionReference.jl (no open tabulated data).
# --------------------------------------------------------------------------

const _BG = (H0 = 10.0, LΩ = 80.0, xb0 = 35.0, Lb = 10.0, d = 2.0, ω = 1.5, η0 = 0.1,
             ρw = 1025.0, g = 9.81)

function _bg_sea()
  p = _BG
  spec = WaveSpec.ContinuousSpectrums.RegularWave(2p.η0, 2π / p.ω)
  ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess = false)
  spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
  k = WaveSpec.AiryWaves.solve_wavenumber(p.ω, p.H0)
  sea = WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, [p.ω], [k], [0.0], p.H0, 1)
  # complex incident amplitude (e^{−iωt}), inlet traction −2ik ϕ̂_I
  ϕI(x) = -im * p.g * p.η0 / p.ω * cosh(k * (x[2] + p.H0)) / cosh(k * p.H0) * exp(im * k * x[1])
  return sea, x -> -2im * k * ϕI(x), k
end

const _BG_SEA, _BG_FIN, _BG_K = _bg_sea()
const _BG_FE = PH.FESpaceConfig(order = 2, vector_type = Vector{ComplexF64})

_bg_tank(nx, ny) = G.TankDomain(L = _BG.LΩ, H = _BG.H0, nx = nx, ny = ny,
  map = x -> VectorValue(x[1], x[2] - _BG.H0),
  body_cutouts = [G.BodyCutout(L = _BG.Lb, d = _BG.d, x₀ = [_BG.xb0, 0.0])])

function _bg_fluid(; incident::Bool)
  bcs = Any[P.RadiationBC(domain = :dΓin), P.RadiationBC(domain = :dΓout)]
  incident && push!(bcs, P.PrescribedInletPotentialBC(domain = :dΓin, forcing = _BG_FIN, quantity = :traction))
  pf = P.PotentialFlow(ρw = _BG.ρw, g = _BG.g, sea_state = _BG_SEA, boundary_conditions = [bcs...],
                       fe = _BG_FE, space_domain_symbol = :Ω)
  fs = P.FreeSurface(ρw = _BG.ρw, g = _BG.g, βₕ = 0.5, fe = _BG_FE, space_domain_symbol = :Γκ)
  P.PhysicsParameters[pf, fs]
end

function _bg_body(; B_ext = nothing)
  p = _BG
  m = p.ρw * p.Lb * p.d                         # freely floating: ρV = m
  P.RigidBody(mass = m, inertia = m * (p.Lb^2 + p.d^2) / 12,
              x_G = [p.xb0 + p.Lb / 2, -p.d / 2], x_ref = [p.xb0 + p.Lb / 2, 0.0],
              hydrostatics = :from_mesh, B_ext = B_ext, ρw = p.ρw, g = p.g, fe = _BG_FE)
end

const _BG_LW = 2π / _BG_K
const _BG_XS_IN = PP.probe_positions(_BG.xb0, :upwave, _BG_LW, _BG.H0)
const _BG_XS_OUT = PP.probe_positions(_BG.xb0 + _BG.Lb, :downwave, _BG_LW, _BG.H0)

@testset "Box eigenfunction reference" begin
  p = _BG
  ref(M) = box_radiation_reference(b = p.Lb / 2, d = p.d, h = p.H0, ω = p.ω, ρ = p.ρw, g = p.g, M = M)
  A40, B40 = ref(40)
  A80, B80 = ref(80)
  @test A40 ≈ A80 rtol = 1e-3                                     # observed ≈ 2e-4
  @test B40 ≈ B80 rtol = 1e-3
  @test A80 ≈ transpose(A80) rtol = 1e-3
  @test B80 ≈ transpose(B80) rtol = 1e-3
  @test iszero(A80[1, 2]) && iszero(B80[2, 3])                   # heave decoupled (symmetry)
  @test all(diag(B80) .> 0)
  # one antisymmetric outgoing wave: the surge–pitch damping block has rank one
  @test abs(det(B80[[1, 3], [1, 3]])) < 1e-3 * B80[1, 1] * B80[3, 3]
end

@testset "Finite-draft barge: radiation against eigenfunction matching" begin
  p = _BG
  cfg = PH.FreqDomainConfig(ω = p.ω)
  Aref, Bref = box_radiation_reference(b = p.Lb / 2, d = p.d, h = p.H0, ω = p.ω, ρ = p.ρw, g = p.g, M = 60)
  err(A, Aref) = abs.(A .- Aref) ./ maximum(abs, Aref)
  A1, B1 = PP.radiation_coefficients(_bg_tank(160, 10), _bg_fluid(incident = false), _bg_body(), cfg)
  A2, B2 = PP.radiation_coefficients(_bg_tank(320, 20), _bg_fluid(incident = false), _bg_body(), cfg)
  @test A2 ≈ transpose(A2) rtol = 1e-10
  @test B2 ≈ transpose(B2) rtol = 1e-10
  # observed on the 320 × 20 mesh: diag A ≤ 0.65 %, diag B ≤ 1.5 %, off-diagonal ≤ 0.5 %
  @test all(abs.(diag(A2) .- diag(Aref)) .< 0.015 .* diag(Aref))
  @test all(abs.(diag(B2) .- diag(Bref)) .< 0.03 .* diag(Bref))
  @test abs(A2[1, 3] - Aref[1, 3]) < 0.01 * maximum(abs, Aref)
  @test abs(B2[1, 3] - Bref[1, 3]) < 0.01 * maximum(abs, Bref)
  # mesh convergence (corner singularity: rate ≈ h^{4/3}, ratio ≈ 2.5)
  @test maximum(err(A2, Aref)) < 0.5 * maximum(err(A1, Aref))
  @test maximum(err(B2, Bref)) < 0.5 * maximum(err(B1, Bref))
end

@testset "Finite-draft barge: Haskind, coupled response and energy balance" begin
  p = _BG
  tank = _bg_tank(320, 20)
  body = _bg_body()
  cfg = PH.FreqDomainConfig(ω = p.ω)
  Cg = PP.group_velocity(_BG_K, p.H0; g = p.g)
  A, B = PP.radiation_coefficients(tank, _bg_fluid(incident = false), body, cfg)

  prob_d = SM.build_problem(tank, _bg_fluid(incident = true), cfg)
  res_d = SM.simulate(prob_d)
  η_inc = PP.fit_wave_components(_BG_XS_IN, PP.sample_probe_line(res_d.solution[2], _BG_XS_IN, 0.0), _BG_K).A
  @test abs(η_inc) ≈ p.η0 rtol = 0.05
  X = PP.excitation_force(prob_d, res_d, body)
  # Haskind (2D body symmetric about x_G): B_jj = |X_j|² / (2ρ g C_g)
  for j in 1:3
    @test B[j, j] ≈ abs2(X[j] / η_inc) / (2 * p.ρw * p.g * Cg) rtol = 5e-3   # observed ≈ 2.7e-4
  end

  for Bpto in (0.0, B[2, 2])
    Bext = zeros(3, 3); Bext[2, 2] = Bpto                   # heave PTO
    bodyc = _bg_body(B_ext = Bext)
    prob_c = SM.build_problem(tank, P.PhysicsParameters[_bg_fluid(incident = true)..., bodyc], cfg)
    res_c = SM.simulate(prob_c)
    ξ = collect(SM.body_motion(res_c, bodyc))
    @test ξ ≈ PP.rigid_body_response(bodyc, SM.get_integration_domains(prob_c), p.ω, A, B, X) rtol = 1e-8
    rt = PP.reflection_transmission_coefficients(res_c.solution[2], _BG_XS_IN, _BG_XS_OUT, _BG_K; y = 0.0)
    P_pto = PP.rigid_body_dissipated_power(ξ, p.ω, Bext)
    A_pto = PP.absorption_from_dissipated_power(P_pto, abs(rt.incident_fit.A), _BG_K, p.H0; ρw = p.ρw, g = p.g)
    @test abs2(rt.R) + abs2(rt.T) + A_pto ≈ 1.0 atol = 2e-3   # observed ≤ 9e-4
  end
end
