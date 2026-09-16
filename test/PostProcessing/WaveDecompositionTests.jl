using Test
using Gridap
using LinearAlgebra
import HydroElasticFEM.PostProcessing as PP

@testset "fit_wave_components — exact recovery" begin
    k = 0.35
    A0 = 1.3 + 0.4im
    B0 = 0.25 - 0.1im
    xs = [2.0, 5.7, 9.3, 14.1, 21.8]   # arbitrary, irregular spacing
    vals = [A0 * exp(im * k * x) + B0 * exp(-im * k * x) for x in xs]

    fit = PP.fit_wave_components(xs, vals, k)
    @test isapprox(fit.A, A0; atol=1e-9)
    @test isapprox(fit.B, B0; atol=1e-9)
    @test fit.residual < 1e-10
    @test fit.cond < 50.0
    @test fit.k == k
    @test length(fit.xs) == length(xs)

    # 2-probe (minimal, exactly-determined) case should also recover exactly,
    # provided the two probes aren't at a degenerate spacing.
    fit2 = PP.fit_wave_components(xs[1:2], vals[1:2], k)
    @test isapprox(fit2.A, A0; atol=1e-9)
    @test isapprox(fit2.B, B0; atol=1e-9)
end

@testset "fit_wave_components — degenerate spacing is flagged via cond" begin
    k = 0.35
    Lhalf = π / k   # spacing at exactly Lwave/2 is the classic degeneracy
    xs = [0.0, Lhalf, 2Lhalf]
    A0, B0 = 1.0 + 0im, 0.3 + 0im
    vals = [A0 * exp(im * k * x) + B0 * exp(-im * k * x) for x in xs]

    fit = PP.fit_wave_components(xs, vals, k)
    # fit_wave_components itself never throws — it just reports a large
    # condition number for the caller (or a higher-level function) to act on.
    @test fit.cond > 1e6
end

@testset "fit_wave_components — argument validation" begin
    @test_throws ErrorException PP.fit_wave_components([1.0], [1.0 + 0im], 0.1)          # N < 2
    @test_throws ErrorException PP.fit_wave_components([1.0, 2.0], [1.0 + 0im], 0.1)      # length mismatch
end

@testset "suggest_probe_offsets" begin
    Lwave = 12.5
    for n in 2:6
        offs = PP.suggest_probe_offsets(Lwave; n=n)
        @test length(offs) == n
        @test issorted(offs)
        @test all(0 .< offs .< Lwave)
        # No offset should land within 2% of a half-wavelength multiple.
        @test all(offs) do o
            r = o / (Lwave / 2)
            abs(r - round(r)) > 0.02
        end
    end
    @test_throws ErrorException PP.suggest_probe_offsets(Lwave; n=1)
    @test_throws ErrorException PP.suggest_probe_offsets(Lwave; n=7)
end

@testset "probe_positions" begin
    Lwave, H = 20.0, 5.0
    edge = 30.0
    xs_up = PP.probe_positions(edge, :upwave, Lwave, H; margin=1.5, n=4)
    xs_down = PP.probe_positions(edge, :downwave, Lwave, H; margin=1.5, n=4)

    @test length(xs_up) == 4
    @test all(xs_up .< edge)
    @test all(xs_down .> edge)
    # Closest probe on each side should clear the margin*H standoff.
    @test edge - maximum(xs_up) >= 1.5 * H
    @test minimum(xs_down) - edge >= 1.5 * H

    @test_throws ErrorException PP.probe_positions(edge, :sideways, Lwave, H)
end

@testset "sample_probe_line" begin
    k = 0.2
    field(p::Point) = ComplexF64(exp(im * k * p[1]))   # depends only on x, ignores y
    xs = [1.0, 4.0, 9.5]
    y = 7.3   # arbitrary — this synthetic field ignores it, exercising the call path
    vals = PP.sample_probe_line(field, xs, y)
    @test vals ≈ [exp(im * k * x) for x in xs]

    # Real-valued fields are promoted to ComplexF64.
    real_field(p::Point) = cos(p[1])
    rvals = PP.sample_probe_line(real_field, xs, y)
    @test rvals ≈ ComplexF64.(cos.(xs))
    @test eltype(rvals) == ComplexF64
end
