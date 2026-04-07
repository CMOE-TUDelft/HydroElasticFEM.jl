"""
    benchmark/runbenchmarks.jl

BenchmarkTools suite for HydroElasticFEM.jl.

Benchmarks three representative problems at reduced mesh sizes so they run
quickly in CI while still exercising the full assembly + solve path:

  1. Empty tank — frequency domain  (plain Gridap formulation)
  2. Floating membrane — frequency domain (plain Gridap formulation)
  3. Empty tank — time domain with sponge-layer damping zones

Usage:
    julia --project=benchmark benchmark/runbenchmarks.jl [--save-baseline]

With --save-baseline the results are written to benchmark/baseline.json.
Without it the results are compared against the baseline (if present) and a
summary is printed to stdout (non-blocking: exit code 0 regardless of deltas).
"""

using BenchmarkTools
using JSON
using Printf

# ── load package + example modules ───────────────────────────────────────────
# Develop-install from the parent directory so that `using HydroElasticFEM`
# picks up the local checkout rather than any registered version.
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using HydroElasticFEM

const PKG_EXAMPLES = joinpath(pkgdir(HydroElasticFEM), "examples")

include(joinpath(PKG_EXAMPLES, "EmptyTankExample.jl"))
include(joinpath(PKG_EXAMPLES, "FloatingMembraneExample.jl"))
include(joinpath(PKG_EXAMPLES, "EmptyTankTimeDomainDampingExample.jl"))

# ── benchmark suite ───────────────────────────────────────────────────────────

const SUITE = BenchmarkGroup()

# Benchmark 1 — empty tank, frequency domain.
# nx=20, ny=4 gives a coarse but non-trivial mesh; order=1 keeps compile time low.
SUITE["empty_tank_freq"] = @benchmarkable(
    EmptyTankExample.run_plain_implementation(; nx=20, ny=4, order=1),
)

# Benchmark 2 — floating membrane, frequency domain.
SUITE["floating_membrane_freq"] = @benchmarkable(
    FloatingMembraneExample.run_plain_implementation(; nx=20, ny=4, order=1),
)

# Benchmark 3 — empty tank, time domain with damping zones.
# Short simulation (tf=0.10, Δt=0.05) with a coarse mesh.
SUITE["empty_tank_time"] = @benchmarkable(
    EmptyTankTimeDomainDampingExample.run_example(;
        nx=8, ny=3, tf=0.10, Δt=0.05, write_vtk=false,
    ),
)

# ── configure and run ─────────────────────────────────────────────────────────

# These are end-to-end solves — each sample takes several seconds; cap at 3.
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60
BenchmarkTools.DEFAULT_PARAMETERS.samples = 3
BenchmarkTools.DEFAULT_PARAMETERS.evals   = 1

# Warm-up / auto-tune each benchmark (skips expensive calibration for evals=1)
println("Tuning benchmarks...")
tune!(SUITE)

println("Running benchmarks...")
results = run(SUITE; verbose=true)

# ── persist or compare ────────────────────────────────────────────────────────

const BASELINE_PATH = joinpath(@__DIR__, "baseline.json")
const SAVE_BASELINE = "--save-baseline" in ARGS

if SAVE_BASELINE
    println("\nSaving baseline → $BASELINE_PATH")
    BenchmarkTools.save(BASELINE_PATH, results)
    println("Baseline saved.")
else
    println("\n── Results ─────────────────────────────────────────────────────")
    for (name, trial) in results
        t = median(trial)
        @printf("  %-35s  %s   mem=%s\n", name,
                BenchmarkTools.prettytime(time(t)),
                BenchmarkTools.prettymemory(memory(t)))
    end

    if isfile(BASELINE_PATH)
        println("\n── Delta vs baseline ────────────────────────────────────────")
        saved   = first(BenchmarkTools.load(BASELINE_PATH))
        verdict = judge(minimum(results), minimum(saved))
        for (name, j) in verdict
            ratio  = j.ratio.time
            symbol = ratio > 1.20 ? "REGRESSION ⚠" :
                     ratio < 0.80 ? "IMPROVEMENT ✓" : "unchanged"
            @printf("  %-35s  ratio=%.2fx  %s\n", name, ratio, symbol)
        end
        println("\n(Performance policy: soft-warn only — no CI gate)")
    else
        println("\nNo baseline at $BASELINE_PATH — rerun with --save-baseline to create one.")
    end
end
