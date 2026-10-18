using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM

# =========================================================================
# StructuralDampingZone validation: damped dispersion of a plate-covered
# channel.
#
# A beam covering the whole surface, with a uniform zone (σ ≡ 1), carries
# waves whose complex wavenumber is the root of the damped host symbol
# (density-normalised, time factor e^{−iωt}):
#
#   (g + C(1 − iωτ)k⁴ − mω² − iωc) k tanh(kH) = ω².
#
# The FE wavenumber is read from equally spaced samples ηₙ of the beam
# deflection: any combination of e^{±ikx} satisfies
# ηₙ₊₁ + ηₙ₋₁ = 2cos(kΔ) ηₙ, so a linear least-squares fit gives k without
# assuming anything about reflections. The samples stay two wavelengths away
# from the inlet and the wall (evanescent fields).
#
# The observable is the shift caused by the damping, δk = k(c, τ) − k(0, 0),
# FE against exact, which removes the (second-order, C/DG P2) mesh
# dispersion error of the undamped wave. Tolerance 1 % of |δk|; the error at
# 24 elements per wavelength is 0.16–0.26 % and halves with h (rate 2).
# =========================================================================

const _SDD = (g=9.81, H=5.0, m=0.01, C=49.0, ω=2.0)

function _sdd_exact_k(c, τ; k=0.4 + 0im)
  (; g, H, m, C, ω) = _SDD
  f(k) = (g + C * (1 - im * ω * τ) * k^4 - m * ω^2 - im * ω * c) *
         k * tanh(k * H) - ω^2
  for _ in 1:50
    df = (f(k + 1e-7) - f(k - 1e-7)) / 2e-7
    k -= f(k) / df
  end
  abs(f(k)) < 1e-10 || error("dispersion root not converged")
  return k
end

function _sdd_fe_k(c, τ; per_wavelength=24, n_wavelengths=8, ny=8)
  (; g, H, m, C, ω) = _SDD
  λ = 2π / real(_sdd_exact_k(0.0, 0.0))
  L = n_wavelengths * λ
  tank = G.TankDomain(L=L, H=H, nx=n_wavelengths * per_wavelength, ny=ny,
    structure_domains=[G.StructureDomain(L=L, x₀=[0.0, H])])
  fe(; kw...) = PH.FESpaceConfig(; order=2, vector_type=Vector{ComplexF64},
                                 kw...)
  zones = iszero(c) && iszero(τ) ? P.StructuralDampingZone[] :
    [P.StructuralDampingZone(profile=x -> 1.0, cᵨ=c, τ=τ)]
  inlet = P.PrescribedInletPotentialBC(domain=:dΓin, quantity=:normal_gradient,
                                       forcing=x -> 1.0 + 0im)
  physics = P.PhysicsParameters[
    P.PotentialFlow(g=g, fe=fe(), boundary_conditions=[inlet]),
    P.EulerBernoulliBeam(L=L, mᵨ=m, EIᵨ=C, g=g, fe=fe(γ=6.0),
                         damping_zones=zones)]
  result = SM.simulate(SM.build_problem(tank, physics, SM.FreqDomainConfig(ω=ω)))
  ηh = result.solution[result.fmap[:η_b]]

  Δ = λ / 8
  η = ComplexF64[ηh(Point(x, H)) for x in (2λ):Δ:(L - 2λ)]
  a = η[2:end-1]
  z = dot(a, η[3:end] .+ η[1:end-2]) / (2 * dot(a, a))
  k = acos(z) / Δ
  return imag(k) < 0 ? -k : k
end

@testset "Damped dispersion of a plate-covered channel" begin
  k0_h, k0 = _sdd_fe_k(0.0, 0.0), _sdd_exact_k(0.0, 0.0)
  @test abs(imag(k0_h)) < 1e-6
  @test abs(k0_h - k0) / abs(k0) < 1e-3

  # dashpot, Kelvin-Voigt, and both together
  for (c, τ) in ((0.5, 0.0), (0.0, 0.1), (0.3, 0.05))
    @testset "cᵨ = $c, τ = $τ" begin
      δk_h = _sdd_fe_k(c, τ) - k0_h
      δk = _sdd_exact_k(c, τ) - k0
      @test imag(δk_h) > 0                          # outgoing wave decays
      @test abs(δk_h - δk) / abs(δk) < 0.01
    end
  end
end
