module KhabakhpashevaBeamJointExample

using Gridap
using Parameters
using Printf
using Plots

using WaveSpec
using HydroElasticFEM: PKG_ROOT
import HydroElasticFEM.Geometry as G
import HydroElasticFEM.ParameterHandler as PH
import HydroElasticFEM.Physics as P
import HydroElasticFEM.Simulation as S

export KhabakhpashevaCaseParams
export run_khabakhpasheva_case
export run_khabakhpasheva_two_cases

@with_kw struct KhabakhpashevaCaseParams
    name::String = "KhabakhpashevaFreqDomain"
    nx::Int = 20
    ny::Int = 5
    order::Int = 4
    ξ::Float64 = 0.0
    vtk_output::Bool = true
    make_plot::Bool = true
end

function _constants()
    Lb = 12.5
    m = 8.36
    EI1 = 47100.0
    EI2 = 471.0
    β = 0.2
    H = 1.1
    α = 0.249

    Ld = Lb
    LΩ = 2 * Ld + 2 * Lb

    x0 = 0.0
    xd_in = x0 + Ld
    xb0 = xd_in + Lb / 2
    xbj = xb0 + β * Lb
    xb1 = xb0 + Lb
    xd_out = LΩ - Ld

    return (; Lb, m, EI1, EI2, β, H, α, Ld, LΩ, x0, xd_in, xb0, xbj, xb1, xd_out)
end

function _physics_constants(c, ξ)
    g = 9.81
    ρ = 1025.0
    d0 = c.m / ρ
    a1 = c.EI1 / ρ
    a2 = c.EI2 / ρ
    kr = ξ * a1 / c.Lb
    return (; g, ρ, d0, a1, a2, kr)
end

function _incident_functions(c, pconst)
    λ = c.α * c.Lb
    k = 2 * pi / λ
    ω = sqrt(pconst.g * k * tanh(k * c.H))
    η0 = 0.01

    H_wave = 2 * η0
    T = 2π / ω
    spec = WaveSpec.ContinuousSpectrums.RegularWave(H_wave, T)
    ds = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(spec; mess=false)
    spread = WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0)
    θ_vec = [0.0]
    sea_state = WaveSpec.AiryWaves.AiryState(ds, spread, 1, 1, [ω], [k], θ_vec, c.H, 1)

    ηin(x) = η0 * exp(im * k * x[1])
    ϕin(x) = -im * (η0 * ω / k) * (cosh(k * x[2]) / sinh(k * c.H)) * exp(im * k * x[1])
    vin(x) = (η0 * ω) * (cosh(k * x[2]) / sinh(k * c.H)) * exp(im * k * x[1])
    vzin(x) = -im * ω * η0 * exp(im * k * x[1])

    return (; λ, k, ω, η0, sea_state, ηin, ϕin, vin, vzin)
end

function _numerics(c, p::KhabakhpashevaCaseParams, ω)
    nx_total = Int(ceil(p.nx / c.β) * ceil(c.LΩ / c.Lb))
    h = c.LΩ / nx_total
    γ = 1.0 * p.order * (p.order - 1) / h
    βh = 0.5
    αh = -im * ω / 9.81 * (1 - βh) / βh
    return (; nx_total, h, γ, βh, αh)
end

function _damping(c, wave)
    μ0 = 2.5
    μ1_in(x) = μ0 * (1.0 - sin(pi / 2 * (x[1]) / c.Ld))
    μ1_out(x) = μ0 * (1.0 - cos(pi / 2 * (x[1] - c.xd_out) / c.Ld))
    μ2_in(x) = μ1_in(x) * wave.k
    μ2_out(x) = μ1_out(x) * wave.k
    return (; μ1_in, μ1_out, μ2_in, μ2_out)
end

function _graded_map(H, ny)
    function f_y(y)
        if y == H
            return H
        end
        i = y / (H / ny)
        H - H / (2.5^i)
    end
    x -> VectorValue(x[1], f_y(x[2]))
end

"""
    run_khabakhpasheva_case(params::KhabakhpashevaCaseParams)

Run one HydroElasticFEM frequency-domain case for the Khabakhpasheva beam+joint setup.
Returns `(xs, η_rel_xs, meta)` where:
- `xs` are normalized coordinates `(x - xb0)/Lb` on `[0,1]`
- `η_rel_xs` is `|η|/η0` sampled along the beam
- `meta` includes derived constants (ω, k, kᵣ, etc.)

Note: current `EulerBernoulliBeam` uses a single `EIᵨ` on `:Γη`.
This example therefore uses `EI1` globally for the beam and represents the
joint via `JointRotationalSpring` at `xbj`.
"""
function run_khabakhpasheva_case(params::KhabakhpashevaCaseParams)
    c = _constants()
    pconst = _physics_constants(c, params.ξ)
    wave = _incident_functions(c, pconst)
    num = _numerics(c, params, wave.ω)
    damp = _damping(c, wave)

    tank = G.TankDomain2D(
        L = c.LΩ,
        H = c.H,
        nx = num.nx_total,
        ny = params.ny,
        map = _graded_map(c.H, params.ny),
        structure_domains = [
            G.StructureDomain1D(L = c.Lb, x₀ = [c.xb0, c.H], domain_symbol = :Γ_beam),
        ],
        damping_zones = [
            G.DampingZone1D(L = c.Ld, x₀ = [0.0, c.H], domain_symbol = :Γ_d_in),
            G.DampingZone1D(L = c.Ld, x₀ = [c.xd_out, c.H], domain_symbol = :Γ_d_out),
        ],
        joint_domains = [
            G.JointDomain1D(location = [c.xbj, c.H], domain_symbol = :dΛj_1, normal_symbol = :n_Λ_j_1),
        ],
    )

    f_in(x) = -wave.vin(x) - im * wave.k * wave.ϕin(x)

    potential = P.PotentialFlow(
        ρw = pconst.ρ,
        g = pconst.g,
        boundary_conditions = [
            P.RadiationBC(domain = :dΓin),
            P.RadiationBC(domain = :dΓout),
            P.PrescribedInletPotentialBC(domain = :dΓin, forcing = f_in, quantity = :traction),
            P.DampingZoneBC(
                domain = :dΓd_1,
                μ₁ = damp.μ1_in,
                μ₂ = damp.μ2_in,
                η_in = wave.ηin,
                vz_in = wave.vzin,
            ),
            P.DampingZoneBC(
                domain = :dΓd_2,
                μ₁ = damp.μ1_out,
                μ₂ = damp.μ2_out,
                η_in = (x -> 0.0 + 0.0im),
                vz_in = (x -> 0.0 + 0.0im),
            ),
        ],
        sea_state = wave.sea_state,
        fe = PH.FESpaceConfig(order = params.order, vector_type = Vector{ComplexF64}),
        space_domain_symbol = :Ω,
    )

    free_surface = P.FreeSurface(
        ρw = pconst.ρ,
        g = pconst.g,
        βₕ = num.βh,
        fe = PH.FESpaceConfig(order = params.order, vector_type = Vector{ComplexF64}),
        space_domain_symbol = :Γκ,
    )

    beam = P.EulerBernoulliBeam(
        L = c.Lb,
        mᵨ = pconst.d0,
        EIᵨ = x -> pconst.a1*(x[1]<c.xbj) + pconst.a2*(x[1]>=c.xbj),
        g = pconst.g,
        joints = [P.JointRotationalSpring(:dΛj_1, :n_Λ_j_1, pconst.kr)],
        fe = PH.FESpaceConfig(order = params.order, vector_type = Vector{ComplexF64}, γ = num.γ),
        space_domain_symbol = :Γη,
    )

    physics = P.PhysicsParameters[potential, free_surface, beam]
    config = S.FreqDomainConfig(ω = wave.ω)

    problem = S.build_problem(tank, physics, config)
    result = S.simulate(problem)

    ϕh, κh, ηh = result.solution

    if params.vtk_output
        trians = S.get_triangulations(problem)
        outdir = joinpath(PKG_ROOT, "data", "VTK", "examples", "KhabakhpashevaBeamJointExample", params.name)
        isdir(outdir) || mkpath(outdir)
        writevtk(trians[:Ω], joinpath(outdir, "omega_field"), cellfields = ["phi_re" => real(ϕh), "phi_im" => imag(ϕh)])
        writevtk(trians[:Γη], joinpath(outdir, "beam_field"), cellfields = ["eta_re" => real(ηh), "eta_im" => imag(ηh)])
        writevtk(trians[:Γκ], joinpath(outdir, "free_surface_field"), cellfields = ["kappa_re" => real(κh), "kappa_im" => imag(κh)])
    end

    ξs = collect(range(0.0, 1.0, length = 400))
    probes = [Point(c.xb0 + ξi * c.Lb, c.H) for ξi in ξs]
    ηvals = ηh(probes)

    xs = ξs
    η_rel_xs = abs.(ηvals) ./ wave.η0

    meta = (
        ω = wave.ω,
        k = wave.k,
        λ = wave.λ,
        kᵣ = pconst.kr,
        EIᵨ = pconst.a1,
        EI2_over_ρ = pconst.a2,
        nx_total = num.nx_total,
        h = num.h,
        γ = num.γ,
        xb0 = c.xb0,
        xbj = c.xbj,
        xb1 = c.xb1,
    )

    return xs, η_rel_xs, meta
end

"""
    run_khabakhpasheva_two_cases(; kwargs...)

Run the two benchmark-style cases:
- ξ = 0.0
- ξ = 625.0

Returns a named tuple with both result curves and metadata.
"""
function run_khabakhpasheva_two_cases(; nx=20, ny=5, order=4, vtk_output=true, make_plot=true)
    case_with = KhabakhpashevaCaseParams(name = "xi_0", nx = nx, ny = ny, order = order, ξ = 0.0, vtk_output = vtk_output, make_plot = make_plot)
    case_without = KhabakhpashevaCaseParams(name = "xi_625", nx = nx, ny = ny, order = order, ξ = 625.0, vtk_output = vtk_output, make_plot = make_plot)

    xs_with, η_with, meta_with = run_khabakhpasheva_case(case_with)
    xs_without, η_without, meta_without = run_khabakhpasheva_case(case_without)

    plt = nothing
    if make_plot
        plt = plot(
            xs_with, η_with,
            lw = 2,
            label = "ξ = 0 (kᵣ = $(round(meta_with.kᵣ, digits=4)))",
            xlabel = "x/L",
            ylabel = "|η|/η₀",
            xlims = (0.0, 1.0),
            legend = :topright,
            title = "HydroElasticFEM: Khabakhpasheva beam with joint",
        )
        plot!(plt, xs_without, η_without, lw = 2, ls = :dash, label = "ξ = 625 (kᵣ = $(round(meta_without.kᵣ, digits=4)))")

        outdir = joinpath(PKG_ROOT, "data", "VTK", "examples", "KhabakhpashevaBeamJointExample")
        isdir(outdir) || mkpath(outdir)
        savefig(plt, joinpath(outdir, "khabakhpasheva_joint_cases.png"))
    end

    return (
        with_joint = (xs = xs_with, η_rel = η_with, meta = meta_with),
        without_joint = (xs = xs_without, η_rel = η_without, meta = meta_without),
        plot = plt,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    @printf("Running HydroElasticFEM Khabakhpasheva beam-joint cases...\n")
    @printf("Warm-up run (nx=2, ny=1, ξ=0)...\n")
    warmup = KhabakhpashevaCaseParams(
        name = "warmup_nx2_ny1",
        nx = 2,
        ny = 1,
        order = 2,
        ξ = 0.0,
        vtk_output = false,
        make_plot = false,
    )
    run_khabakhpasheva_case(warmup)
    results = run_khabakhpasheva_two_cases()
    @printf("Done. Generated %d points per case.\n", length(results.with_joint.xs))
end

end # module
