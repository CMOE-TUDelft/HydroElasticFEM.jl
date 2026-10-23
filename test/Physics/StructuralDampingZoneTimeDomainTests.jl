using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.Physics as P
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Simulation as SM

# =========================================================================
# StructuralDampingZone in the time domain: a plate-covered channel driven
# by an inlet flux cos(ωt), with a ramped zone (dashpot + Kelvin-Voigt) at
# the far end.
#
# The time-domain run starts on the frequency-domain steady state
# (η = Re η̂, ηₜ = Re(−iωη̂), ηₜₜ = Re(−ω²η̂) for every field), so no
# start-up transient is excited: if the time-domain damping form equals the
# frequency-domain one, the solution stays on Re(η̂ e^{−iωt}) up to the
# time-integration error. Both solves share the mesh, so the spatial error
# cancels. Checked over one period: the error is second order in Δt
# (trapezoidal generalised-α, ρ∞ = 1) and below 1 % at 60 steps per period
# (≈ 0.4 %).
# =========================================================================

const _SDT = (g=9.81, H=5.0, m=0.01, C=49.0, ω=2.0, λ=2π / 0.3849)

function _sdt_setup(vector_type, forcing)
  (; g, H, m, C, λ) = _SDT
  L = 3λ
  tank = G.TankDomain(L=L, H=H, nx=24, ny=4,
    structure_domains=[G.StructureDomain(L=L, x₀=[0.0, H])])
  fe(; kw...) = PH.FESpaceConfig(; order=2, vector_type=vector_type, kw...)
  zone = P.StructuralDampingZone(
    profile=P.ramp_profile(x_start=L - λ, x_end=L), cᵨ=3.0, τ=0.1)
  inlet = P.PrescribedInletPotentialBC(domain=:dΓin, quantity=:normal_gradient,
                                       forcing=forcing)
  physics = P.PhysicsParameters[
    P.PotentialFlow(g=g, fe=fe(), boundary_conditions=[inlet]),
    P.EulerBernoulliBeam(L=L, mᵨ=m, EIᵨ=C, g=g, fe=fe(γ=6.0),
                         damping_zones=[zone])]
  return tank, physics, L
end

# Relative deviation of the time-domain run from Re(η̂ e^{−iωt}), sampled on
# the beam over `periods` periods with `steps` steps per period.
function _sdt_deviation(û, η̂, xs, steps; periods=1)
  (; H, ω) = _SDT
  T = 2π / ω
  tank, physics, _ = _sdt_setup(Vector{Float64}, t -> (x -> cos(ω * t)))
  tf = periods * T
  prob = SM.build_problem(tank, physics, PH.TimeDomainConfig(t₀=0.0, tf=tf))
  fmap = SM.get_field_map(prob)
  X0 = SM.get_test_fe_space(prob)(0.0)
  init(a) = (u = FEFunction(X0, real.(a * û)); [u[i] for i in 1:length(fmap)])
  tc = PH.TimeConfig(Δt=T / steps, t₀=0.0, tf=tf, ρ∞=1.0,
                     u0=init(1.0), u0t=init(-im * ω), u0tt=init(-ω^2))
  sym = fmap[:η_b]
  deviation = 0.0
  for (t, uh) in SM.simulate(prob, tc).solution
    ηt = [uh[sym](Point(x, H)) for x in xs]
    deviation = max(deviation, norm(ηt .- real.(η̂ .* exp(-im * ω * t))))
  end
  return deviation / norm(η̂)
end

@testset "Time domain stays on the frequency-domain steady state" begin
  (; H, ω, λ) = _SDT
  tank, physics, L = _sdt_setup(Vector{ComplexF64}, x -> 1.0 + 0im)
  result = SM.simulate(SM.build_problem(tank, physics, SM.FreqDomainConfig(ω=ω)))
  ηh = result.solution[result.fmap[:η_b]]
  xs = collect(range(0.3λ, L, length=15))   # includes the zone
  η̂ = ComplexF64[ηh(Point(x, H)) for x in xs]
  û = get_free_dof_values(result.solution)

  e30 = _sdt_deviation(û, η̂, xs, 30)
  e60 = _sdt_deviation(û, η̂, xs, 60)
  @test e60 < 0.01
  @test log2(e30 / e60) > 1.8                  # second order in Δt
end
