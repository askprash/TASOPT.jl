"""
    generate_engine_benchmark_baseline.jl

Standalone script to regenerate the engine-sweep performance baseline.

## Usage

    julia --project=. test/generate_engine_benchmark_baseline.jl

## When to regenerate

Run this script after intentional algorithm changes that shift allocations
or runtime (e.g. reducing heap allocations, restructuring loops).  Commit
the updated baseline alongside the code change so reviewers can audit the
improvement.

## What it produces

    test/fixtures/engine_benchmark_baseline.toml

A TOML file with the allocation count (deterministic) and median elapsed
time (hardware-dependent reference) for `run_engine_sweep` on the default
aircraft model.
"""

using TASOPT, TOML, Dates

@warn """
*** REGENERATING ENGINE SWEEP BENCHMARK BASELINE ***
This script rewrites test/fixtures/engine_benchmark_baseline.toml.
Run it only when a deliberate performance change is intended.
Commit the updated baseline with an explanatory git message.
"""

println("Loading and sizing default aircraft...")
ac = TASOPT.load_default_model()
TASOPT.size_aircraft!(ac; printiter=false)

# Warmup: two calls to reach stable post-JIT allocation count.
# engine_state_to_pare! now writes 22 design-constant fields to bare pare
# (tasopt-j9l.45.14.4); Julia needs two passes to fully specialise all paths.
println("Warmup run 1 (JIT)...")
_ = TASOPT.engine.run_engine_sweep(ac)
println("Warmup run 2 (stable)...")
_ = TASOPT.engine.run_engine_sweep(ac)

# Allocation count — deterministic post-JIT.
println("Measuring allocations...")
alloc = @allocations TASOPT.engine.run_engine_sweep(ac)

# Elapsed time — median of 7 post-warmup runs.
println("Measuring elapsed time (7 runs)...")
times = [(@elapsed TASOPT.engine.run_engine_sweep(ac)) for _ in 1:7]
sort!(times)
med_time = times[4]  # median of 7

println("alloc_count = $(alloc)")
println("elapsed_med_s = $(med_time)")

baseline_path = joinpath(@__DIR__, "fixtures", "engine_benchmark_baseline.toml")

# Build the new baseline dict and write it.
d = Dict{String,Any}(
    "meta" => Dict{String,Any}(
        "julia_version" => string(VERSION),
        "date"          => string(Dates.today()),
        "description"   => "Baseline captured post-warmup on ralph-bd-afk Docker sandbox",
    ),
    "benchmark" => Dict{String,Any}(
        "alloc_count" => alloc,
        "elapsed_s"   => med_time,
    ),
)

open(baseline_path, "w") do io
    println(io, "# Engine sweep performance baseline — tasopt-j9l.7")
    println(io, "#")
    println(io, "# To regenerate after intentional performance changes, run:")
    println(io, "#   julia --project=. test/generate_engine_benchmark_baseline.jl")
    println(io, "#")
    println(io, "# Allocation check tolerance: 5%  — hard @test failure on regression.")
    println(io, "# Runtime check tolerance:   10%  — informational log only (hardware-dependent).")
    println(io)
    TOML.print(io, d)
end

println("Baseline written to $(baseline_path)")
