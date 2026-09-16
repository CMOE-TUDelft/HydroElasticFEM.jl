using Test
using Gridap
import HydroElasticFEM.PostProcessing as PP

@testset "absorption_coefficient" begin
    @test PP.absorption_coefficient(0.0 + 0im, 1.0 + 0im) ≈ 0.0       # fully transmitted
    @test PP.absorption_coefficient(1.0 + 0im, 0.0 + 0im) ≈ 0.0       # fully reflected
    @test PP.absorption_coefficient(0.0 + 0im, 0.0 + 0im) ≈ 1.0       # fully absorbed
    R, T = 0.6 + 0.1im, 0.5 - 0.2im
    @test PP.absorption_coefficient(R, T) ≈ 1 - abs2(R) - abs2(T)
end

@testset "reflection_transmission_coefficients — synthetic scattering field" begin
    # Build a fabricated (but physically-shaped) scattering field: incident +
    # reflected wave upstream of x_split, purely transmitted wave downstream.
    # R, T here are chosen freely (not required to satisfy energy balance) —
    # this test only checks that the fit recovers the prescribed values, not
    # that A comes out to any particular number.
    k = 0.22
    η0 = 0.4
    R_true = 0.35 - 0.05im
    T_true = 0.80 + 0.10im
    x_split = 50.0

    function κ_synthetic(p::Point)
        x = p[1]
        if x < x_split
            return η0 * exp(im * k * x) + R_true * η0 * exp(-im * k * x)
        else
            return T_true * η0 * exp(im * k * x)
        end
    end

    Lwave = 2π / k
    xs_in = PP.probe_positions(x_split, :upwave, Lwave, 10.0; n=4)
    xs_out = PP.probe_positions(x_split, :downwave, Lwave, 10.0; n=4)

    rt = PP.reflection_transmission_coefficients(κ_synthetic, xs_in, xs_out, k; y=0.0, η0=η0)

    @test isapprox(rt.R, R_true; atol=1e-8)
    @test isapprox(rt.T, T_true; atol=1e-8)
    @test isapprox(rt.A, 1 - abs2(R_true) - abs2(T_true); atol=1e-8)
    @test isapprox(abs(rt.incident_fit.A), η0; atol=1e-8)
    @test isapprox(abs(rt.transmitted_fit.B), 0.0; atol=1e-8)   # no spurious back-reflection by construction

    # Raw-vector method should agree exactly with the field method.
    vals_in = PP.sample_probe_line(κ_synthetic, xs_in, 0.0)
    vals_out = PP.sample_probe_line(κ_synthetic, xs_out, 0.0)
    rt2 = PP.reflection_transmission_coefficients(xs_in, vals_in, xs_out, vals_out, k)
    @test rt2.R ≈ rt.R
    @test rt2.T ≈ rt.T
end

@testset "reflection_transmission_coefficients — energy-conserving case" begin
    # |R|^2 + |T|^2 = 1 exactly by construction ⇒ A should be ≈ 0.
    k = 0.15
    η0 = 1.0
    Rmag, Tmag = 0.6, sqrt(1 - 0.6^2)
    R_true = Rmag * cis(0.7)
    T_true = Tmag * cis(-1.1)
    x_split = 20.0

    function κ_synthetic(p::Point)
        x = p[1]
        x < x_split ?
            η0 * exp(im * k * x) + R_true * η0 * exp(-im * k * x) :
            T_true * η0 * exp(im * k * x)
    end

    Lwave = 2π / k
    xs_in = PP.probe_positions(x_split, :upwave, Lwave, 8.0; n=5)
    xs_out = PP.probe_positions(x_split, :downwave, Lwave, 8.0; n=5)

    rt = PP.reflection_transmission_coefficients(κ_synthetic, xs_in, xs_out, k; y=0.0)
    @test isapprox(rt.A, 0.0; atol=1e-8)
    @test isapprox(abs2(rt.R) + abs2(rt.T), 1.0; atol=1e-8)
end

@testset "reflection_transmission_coefficients — ill-conditioned probes raise an error" begin
    k = 0.3
    Lhalf = π / k
    xs_in = [0.0, Lhalf, 2Lhalf]     # degenerate spacing on purpose
    xs_out = [100.0, 100.0 + Lhalf, 100.0 + 2Lhalf]
    vals_in = fill(1.0 + 0im, 3)
    vals_out = fill(1.0 + 0im, 3)

    @test_throws ErrorException PP.reflection_transmission_coefficients(
        xs_in, vals_in, xs_out, vals_out, k,
    )
end

@testset "group_velocity limits" begin
    g = 9.81
    H = 50.0

    # Deep water: kH ≫ 1 ⇒ cg → ω / (2k)  (half the phase speed)
    k_deep = 2.0
    ω_deep = sqrt(g * k_deep * tanh(k_deep * H))
    @test isapprox(PP.group_velocity(k_deep, H; g=g), ω_deep / (2k_deep); rtol=1e-3)

    # Shallow water: kH ≪ 1 ⇒ cg → sqrt(g H)
    k_shallow = 1e-4
    @test isapprox(PP.group_velocity(k_shallow, H; g=g), sqrt(g * H); rtol=1e-3)
end

@testset "energy_flux and dissipation-based absorption" begin
    g, ρw, H = 9.81, 1025.0, 10.0
    k = 0.2
    η0 = 0.3
    cg = PP.group_velocity(k, H; g=g)
    P_in = PP.energy_flux(η0, k, H; ρw=ρw, g=g)
    @test isapprox(P_in, 0.5 * ρw * g * η0^2 * cg)

    # A resonator dissipating exactly P_in should give A ≈ 1.
    ω = sqrt(g * k * tanh(k * H))
    C = 5.0e3
    qh = 2.0 + 0.0im
    # Choose qh magnitude so that resonator_dissipated_power(qh, ω, C) == P_in
    qh_scaled = sqrt(2 * P_in / (C * ω^2)) + 0im
    P_diss = PP.resonator_dissipated_power(qh_scaled, ω, C)
    @test isapprox(P_diss, P_in; rtol=1e-10)
    @test isapprox(PP.absorption_from_dissipated_power(P_diss, η0, k, H; ρw=ρw, g=g), 1.0; rtol=1e-10)
end
