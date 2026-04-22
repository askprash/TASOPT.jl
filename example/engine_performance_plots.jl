# Engine performance plots — standalone script
#
# Generates two PNG figures from the default TASOPT aircraft model:
#
#   engine_station_profiles.png
#       Two-panel figure: total-temperature (K) and total-pressure (Pa) along
#       the gas path for all 16 mission points.
#
#   engine_performance.png
#       Seven-panel figure: net thrust, TSFC, bypass ratio, core mass flow,
#       and corrected fan/LPC/HPC spool speeds versus mission point.
#
# Usage:
#   julia --project=. example/engine_performance_plots.jl [OUTPUT_DIR]
#
# OUTPUT_DIR defaults to the current directory if not provided.
# The directory is created if it does not exist.
#
# Environment:
#   GKSwstype=nul  — set this to suppress window creation in headless mode
#                    (already the default in the ralph-bd-afk sandbox).

using TASOPT

outdir = isempty(ARGS) ? "." : ARGS[1]
mkpath(outdir)

println("Loading default aircraft model...")
ac = TASOPT.load_default_model()
size_aircraft!(ac; printiter=false)
println("Aircraft sized successfully.")

println("Running engine sweep (all 16 mission points)...")
sweep = TASOPT.engine.run_engine_sweep(ac)
println("  $(length(sweep.ip_indices)) points: $(join(sweep.ip_labels, ", "))")

# ── Figure 1: station Tt / pt profiles ────────────────────────────────────
println("Generating station profile plots...")
p_profiles = TASOPT.plot_engine_station_profiles(sweep;
    size=(900, 600), dpi=150)
path_profiles = joinpath(outdir, "engine_station_profiles.png")
savefig(p_profiles, path_profiles)
println("  Saved: $path_profiles")

# ── Figure 2: performance + spool speeds vs mission point ─────────────────
println("Generating performance plots...")
p_perf = TASOPT.plot_engine_performance(sweep; ac=ac,
    size=(1100, 900), dpi=150)
path_perf = joinpath(outdir, "engine_performance.png")
savefig(p_perf, path_perf)
println("  Saved: $path_perf")

println("Done.")
