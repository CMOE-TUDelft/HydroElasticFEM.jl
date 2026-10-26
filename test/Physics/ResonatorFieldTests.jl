using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM
import HydroElasticFEM.Simulation.FEOperators as FO

# =========================================================================
# ResonatorField: homogenised resonators on a structure.
#
# (1) Dispersion. A beam covering a whole channel carries a uniform field.
#     Eliminating q gives the beam the effective mass density (normalised,
#     time factor e^{−iωt})
#       m_eff = mᵨ + (nM/ρw)(K − iωC)/(K − iωC − ω²M),
#     so k(ω) is the root of (g + EIᵨk⁴ − m_eff ω²) k tanh(kH) = ω²; without
#     damping, ω must lie on one of the closed-form branches ω±(k)
#     (dispersion_theory.tex, eq:omega_k_x_closed_form). The FE wavenumber
#     is read from equally spaced deflection samples (as in
#     StructuralDampingZoneDispersionTests.jl).
# (2) Consistency with a ResonatorArray of spacing a = 1/n (O(a²)), and
#     matrix-level checks (symmetry, rigid mode, energy sign, total mass).
# (3) 3D: plate + field on a one-cell-wide channel (ν = 0) == 2D beam +
#     field to round-off.
# =========================================================================

const _RF = (g=9.81, H=5.0, m=0.01, EI=49.0, ω=2.0, ρw=1025.0)

_rf_fe(; kw...) = PH.FESpaceConfig(; order=2, vector_type=Vector{ComplexF64}, kw...)

function _rf_field(; μ=0.5, ωr=2.5, ζ=0.0, n=1.0, kw...)
  M = μ * _RF.ρw / n
  P.ResonatorField(; M=M, K=M * ωr^2, C=2ζ * M * ωr, n=n, ρw=_RF.ρw,
                   fe=_rf_fe(), kw...)
end

function _rf_exact_k(rf; k=0.4 + 0im)
  (; g, H, m, EI, ω) = _RF
  madd = isnothing(rf) ? 0.0 : P.effective_mass_density(rf, ω)
  f(k) = (g + EI * k^4 - (m + madd) * ω^2) * k * tanh(k * H) - ω^2
  for _ in 1:60
    k -= f(k) / ((f(k + 1e-7) - f(k - 1e-7)) / 2e-7)
  end
  abs(f(k)) < 1e-10 || error("dispersion root not converged")
  return k
end

# Closed-form undamped branches ω±(k) of the paper (with M_f = 1/(k tanh kH)).
function _rf_branches(k, rf)
  (; g, H, m, EI) = _RF
  ωr2 = rf.K / rf.M
  Kk = g + EI * k^2 * k^2
  Mt = m + 1 / (k * tanh(k * H))
  Mr = rf.n * rf.M / rf.ρw
  Σ = Kk + ωr2 * (Mt + Mr)
  d = sqrt(Σ^2 - 4 * ωr2 * Kk * Mt)
  return sqrt((Σ - d) / (2Mt)), sqrt((Σ + d) / (2Mt))
end

_fit_k(η, Δ) = begin
  a = η[2:end-1]
  z = dot(a, η[3:end] .+ η[1:end-2]) / (2 * dot(a, a))
  k = acos(z) / Δ
  imag(k) < 0 ? -k : k
end

# Beam over the whole channel, inlet flux, wall at the end. `extra` is the
# resonator entity (field or array) or nothing; `rdoms` its point domains.
function _rf_channel(extra; λ, n_wavelengths, per_wavelength, ny=8,
                     rdoms=G.ResonatorDomain[])
  (; g, H, m, EI, ω) = _RF
  L = n_wavelengths * λ
  tank = G.TankDomain(L=L, H=H, nx=round(Int, n_wavelengths * per_wavelength),
    ny=ny, structure_domains=[G.StructureDomain(L=L, x₀=[0.0, H])],
    resonator_domains=rdoms)
  inlet = P.PrescribedInletPotentialBC(domain=:dΓin, quantity=:normal_gradient,
                                       forcing=x -> 1.0 + 0im)
  physics = P.PhysicsParameters[
    P.PotentialFlow(g=g, fe=_rf_fe(), boundary_conditions=[inlet]),
    P.EulerBernoulliBeam(L=L, mᵨ=m, EIᵨ=EI, g=g, fe=_rf_fe(γ=6.0))]
  isnothing(extra) || push!(physics, extra)
  result = SM.simulate(SM.build_problem(tank, physics, SM.FreqDomainConfig(ω=ω)))
  return result, L
end

function _rf_fe_k(rf; per_wavelength=24, n_wavelengths=8)
  λ = 2π / real(_rf_exact_k(rf))
  result, L = _rf_channel(rf; λ, n_wavelengths, per_wavelength)
  ηh = result.solution[result.fmap[:η_b]]
  Δ = λ / 8
  η = ComplexF64[ηh(Point(x, _RF.H)) for x in (2λ):Δ:(L - 2λ)]
  return _fit_k(η, Δ), result
end

@testset "ResonatorField" begin

  @testset "parameters and effective mass" begin
    rf = _rf_field(ζ=0.05)
    @test P.variable_symbol(rf) == :q
    @test P.has_damping_form(rf) && !P.has_rhs_form(rf)
    @test !P.has_damping_form(_rf_field())
    @test P.effective_mass_density(rf, 0.0) ≈ 0.5            # ω → 0: rigid
    @test abs(P.effective_mass_density(rf, 1e5)) < 1e-5      # ω → ∞: decoupled
    beam = P.EulerBernoulliBeam(L=1.0, mᵨ=1.0, EIᵨ=1.0)
    @test P.hosts(rf, beam) && P.has_stiffness_form(rf, beam)
    @test !P.hosts(P.ResonatorField(M=1.0, K=1.0, space_domain_symbol=:Γ_a), beam)
    @test !P.has_stiffness_form(beam, rf)                     # ordered pair
    # graded population
    rg = P.ResonatorField(M=1.0, K=x -> 4.0, n=x -> 2.0, ρw=1.0)
    @test P.effective_mass_density(rg, 1.0, Point(0.0, 0.0)) ≈ 2.0 * 4 / 3
  end

  @testset "dispersion of a beam with a uniform field" begin
    k0_h, _ = _rf_fe_k(nothing)
    k0 = _rf_exact_k(nothing)

    # undamped: real k, matches the root and the closed-form branch
    rf = _rf_field()
    k_h, res = _rf_fe_k(rf)
    k = _rf_exact_k(rf)
    @test abs(imag(k_h)) < 1e-6
    err = abs(k_h - k) / abs(k)
    @info "ResonatorField dispersion (undamped)" k0 k real(k_h) err
    @test err < 2e-3
    @test real(k) > 1.3 * real(k0)                   # the field matters
    ωm, ωp = _rf_branches(real(k_h), rf)
    @test min(abs(ωm - _RF.ω), abs(ωp - _RF.ω)) / _RF.ω < 1e-3
    # q and η share the FE space: the q equation gives q = K/(K − ω²M) η
    qh, ηh = res.solution[res.fmap[:q]], res.solution[res.fmap[:η_b]]
    T = rf.K / (rf.K - _RF.ω^2 * rf.M)
    @test norm(get_free_dof_values(qh) - T * get_free_dof_values(ηh)) <
          1e-8 * norm(get_free_dof_values(qh))

    # damped: compare the shift caused by the field (removes mesh dispersion)
    rd = _rf_field(ζ=0.05)
    kd_h, _ = _rf_fe_k(rd)
    kd = _rf_exact_k(rd)
    δ_err = abs((kd_h - k0_h) - (kd - k0)) / abs(kd - k0)
    @info "ResonatorField dispersion (ζ = 0.05)" kd kd_h δ_err
    @test imag(kd_h) > 0
    @test abs(kd_h - kd) / abs(kd) < 2e-3
    @test δ_err < 0.01
  end

  @testset "matrix checks: symmetry, rigid mode, energy, mass" begin
    (; H) = _RF
    L, nx = 10.0, 20
    rf = _rf_field(ζ=0.05, n=0.7)
    tank = G.TankDomain(L=L, H=H, nx=nx, ny=2,
      structure_domains=[G.StructureDomain(L=L, x₀=[0.0, H])])
    fe = PH.FESpaceConfig(order=2, vector_type=Vector{Float64})
    beam = P.EulerBernoulliBeam(L=L, mᵨ=0.01, EIᵨ=49.0, fe=PH.FESpaceConfig(
      order=2, vector_type=Vector{Float64}, γ=6.0))
    rf = P.ResonatorField(M=rf.M, K=rf.K, C=rf.C, n=rf.n, ρw=rf.ρw, fe=fe)
    tr = G.build_triangulations(tank, G.build_model(tank))
    dom = G.get_integration_domains(tr; degree=4)
    X, Y, fmap = SM.build_fe_spaces(P.PhysicsParameters[beam, rf], tr,
                                    PH.TimeDomainConfig())
    U = X(0.0)
    fm(u) = FO.FieldMap(u, fmap)
    Kf = assemble_matrix((x, y) -> P.stiffness(rf, dom, fm(x), fm(y)) +
                                   P.stiffness(rf, beam, dom, fm(x), fm(y)), U, Y)
    Cf = assemble_matrix((x, y) -> P.damping(rf, dom, fm(x), fm(y)) +
                                   P.damping(rf, beam, dom, fm(x), fm(y)), U, Y)
    Mf = assemble_matrix((x, y) -> P.mass(rf, dom, fm(x), fm(y)), U, Y)
    @test norm(Kf - Kf') < 1e-10 * norm(Kf)
    @test norm(Cf - Cf') < 1e-10 * norm(Cf)
    one_ = ones(size(Kf, 1))                 # q = η = 1: springs unloaded
    @test norm(Kf * one_) < 1e-10 * norm(Kf)
    @test norm(Cf * one_) < 1e-10 * norm(Cf)
    for _ in 1:5
      z = randn(size(Kf, 1))
      @test dot(z, Kf * z) ≥ -1e-10 * norm(Kf) * dot(z, z)
      @test dot(z, Cf * z) ≥ -1e-10 * norm(Cf) * dot(z, z)
    end
    # total mass ∫ nM/ρw dΓ = (number of resonators)·M/ρw
    @test dot(one_, Mf * one_) ≈ rf.n * L * rf.M / rf.ρw
    # Cf = (C/K) Kf for constant parameters
    @test norm(Cf - (rf.C / rf.K) * Kf) < 1e-10 * norm(Kf)
  end

  @testset "consistency with a ResonatorArray of spacing a = 1/n" begin
    # Undamped, below resonance, 8 and 16 resonators per wavelength (ka ≈
    # 0.8, 0.4), placed inside cells (25 cells per wavelength). The array →
    # field difference must decay like (ka)².
    μ, ωr = 0.2, 3.0
    rf = _rf_field(; μ, ωr)
    λ = 2π / real(_rf_exact_k(rf))
    nw, pw = 3, 25
    xs = range(0.5λ, nw * λ - 0.5λ, length=41)
    sample(res) = (η = res.solution[res.fmap[:η_b]];
                   ComplexF64[η(Point(x, _RF.H)) for x in xs])
    rfield, L = _rf_channel(rf; λ, n_wavelengths=nw, per_wavelength=pw)
    ηf = sample(rfield)
    diffs = Float64[]
    for per in (8, 16)
      a = λ / per
      N = nw * per
      ra = P.resonator_array(N, μ * _RF.ρw * a, μ * _RF.ρw * a * ωr^2, 0.0;
                             ρw=_RF.ρw, fe=_rf_fe(order=1))
      rdoms = [G.ResonatorDomain(location=[(i - 0.5) * a, _RF.H],
                                 delta_symbol=Symbol("δ_p$i")) for i in 1:N]
      rarr, _ = _rf_channel(ra; λ, n_wavelengths=nw, per_wavelength=pw, rdoms)
      push!(diffs, norm(sample(rarr) - ηf) / norm(ηf))
    end
    rate = log2(diffs[1] / diffs[2])
    @info "ResonatorField vs ResonatorArray" diffs rate
    @test diffs[2] < diffs[1]
    @test rate > 1.5
    @test diffs[2] < 0.03
  end

  @testset "time domain: floating membrane with a field" begin
    tank = G.TankDomain(L=20.0, H=2.0, nx=40, ny=2,
      structure_domains=[G.StructureDomain(L=4.0, x₀=[4.0, 2.0])])
    fe = PH.FESpaceConfig(vector_type=Vector{Float64})
    mem = P.Membrane(L=4.0, mᵨ=0.9, Tᵨ=98.1, fe=fe)
    rf = P.ResonatorField(M=100.0, K=500.0, C=5.0, n=2.0, fe=fe)
    inlet = P.PrescribedInletPotentialBC(domain=:dΓin, quantity=:traction,
                                         forcing=(t -> (x -> 0.1 * sin(2t))))
    physics = P.PhysicsParameters[P.PotentialFlow(fe=fe, boundary_conditions=[inlet]),
                                  P.FreeSurface(fe=fe), mem, rf]
    @test count(p -> p[1] === rf, SM.detect_couplings(physics)) == 1
    cfg = PH.TimeDomainConfig(t₀=0.0, tf=0.5)
    tc = PH.TimeConfig(Δt=0.1, tf=0.5, u0=zeros(4), u0t=zeros(4), u0tt=zeros(4))
    prob = SM.build_problem(tank, physics, cfg; tconfig=tc)
    qmax = 0.0
    res = SM.simulate(prob, tc)
    for (_, uh) in res.solution
      q = get_free_dof_values(uh[res.fmap[:q]])
      @test all(isfinite, q)
      qmax = max(qmax, maximum(abs, q))
    end
    @test qmax > 0
  end

  @testset "3D plate + field == 2D beam + field (ν = 0)" begin
    (; g, H, m, EI, ω) = _RF
    λ = 2π / real(_rf_exact_k(_rf_field()))
    L, nx, nz = 3λ, 24, 4
    dx = L / nx
    hb = 0.01
    rf = _rf_field(ζ=0.05)
    inlet = P.PrescribedInletPotentialBC(domain=:dΓin, quantity=:normal_gradient,
                                         forcing=x -> 1.0 + 0im)
    xs = collect(range(0.3λ, L, length=15))

    tank2 = G.TankDomain(L=L, H=H, nx=nx, ny=nz,
      structure_domains=[G.StructureDomain(L=L, x₀=[0.0, H])])
    beam = P.EulerBernoulliBeam(L=L, mᵨ=m, EIᵨ=EI, g=g, fe=_rf_fe(γ=6.0))
    r2 = SM.simulate(SM.build_problem(tank2,
      P.PhysicsParameters[P.PotentialFlow(g=g, fe=_rf_fe(),
                                          boundary_conditions=[inlet]), beam, rf],
      SM.FreqDomainConfig(ω=ω)))
    v2 = ComplexF64[r2.solution[r2.fmap[:q]](Point(x, H)) for x in xs]

    tank3 = G.TankDomain(L=L, W=dx, H=H, nx=nx, ny=1, nz=nz,
      structure_domains=[G.StructureDomain(L=L, W=dx, x₀=[0.0, 0.0, H])])
    plate = P.KirchhoffLovePlate(E=12 * EI / hb^3, ν=0.0, hb=hb, ρ=1.0,
                                 ρb=m / hb, g=g, fe=_rf_fe(γ=6.0))
    r3 = SM.simulate(SM.build_problem(tank3,
      P.PhysicsParameters[P.PotentialFlow(g=g, dim=3, fe=_rf_fe(),
                                          boundary_conditions=[inlet]), plate, rf],
      SM.FreqDomainConfig(ω=ω)))
    q3 = r3.solution[r3.fmap[:q]]
    v3 = ComplexF64[q3(Point(x, dx / 2, H)) for x in xs]
    @test all(isfinite, v3) && norm(v3) > 0
    @test norm(v3 - v2) / norm(v2) < 1e-10
  end
end
