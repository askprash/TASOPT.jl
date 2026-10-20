"""
benchmark_optimizers_fair.jl

Fair comparison of derivative-free optimizers on the TASOPT 737-800 problem
defined by tasopt2.16_clean/runs/737/737_opt_3000.tas (Fortran simpop).

What is matched to Fortran (fobj.f / voptset):
  * the same 11 design variables and the same initial perturbation steps
  * altitude steps in metres (10 m), not feet
  * warm-started sizing: `initwgt=true` after the first evaluation
  * the same exterior penalties: lBF (5 Wpay), fuel volume (5 Wpay), span (25 Wpay),
    top-of-climb gradient (1 Wpay), plus the ad hoc sweep, lambdas, lambdat penalties
  * the lambdat >= 0.1 clamp (applied to the aircraft; see note below)
  * Tt4CR on ip cruise1:cruisen, Tt4TO on static/rotate/takeoff, CL/rcls/rclt on climb1+1:descentn-1,
    wing sweep also written to the horizontal tail sweep (igsweeph)

Optimizers: NLopt LN_NELDERMEAD, LN_SBPLX, LN_BOBYQA, and Optim.jl NelderMead with
FixedParameters(1, 1.4, 0.7, 0.7) (simpop's coefficients) inside a restart wrapper.

Starting points: the baseline point and 3 perturbed points (fixed seed).

Run (from the repository root):
  julia --project=test/optimizer_bench test/benchmark_optimizers_fair.jl
Environment variables:
  TASOPT_BENCH_MAXEVAL   evaluation cap per run (default 6000)
  TASOPT_BENCH_OPTS      comma list of optimizers to run (default all)
  TASOPT_BENCH_STARTS    comma list of start ids to run (default 0,1,2,3)
  TASOPT_BENCH_OUT       output directory (default test/benchmark_results)

Outputs (in the output directory): runs.json, runs_summary.csv, evals_<run>.csv (per-evaluation trace).
Summarize with test/summarize_optimizer_benchmark.jl.

Note on the lambdat clamp: Fortran writes the clamped value back into the simplex vertex. The
optimizer libraries own their vertices, so here the clamp is applied to the aircraft only; below
lambdat = 0.1 the objective is flat in that coordinate.
"""

using TASOPT
using NLopt
using Optim
using Printf
using Random

include(__TASOPTindices__)

const ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUT_TOML = joinpath(ROOT, "example", "737_compare_input.toml")
const OUTDIR = get(ENV, "TASOPT_BENCH_OUT", joinpath(@__DIR__, "benchmark_results"))
const MAXEVAL = parse(Int, get(ENV, "TASOPT_BENCH_MAXEVAL", "6000"))
const FT = 0.3048

# ── Problem definition ───────────────────────────────────────────
const VAR_NAMES = ["CL", "AR", "sweep_deg", "hboxo", "hboxs", "lambdat",
                   "rcls", "rclt", "alt_m", "Tt4CR_K", "Tt4TO_K"]
# initial perturbations (= initial simplex axis steps), from 737_opt_3000.tas
const DX = [0.001, 0.05, 0.025, 0.00025, 0.00025, 0.001, 0.001, 0.001, 10.0, 0.5, 0.25]
# half-widths used to perturb the extra start points
const START_SCALE = [0.03, 0.5, 2.0, 0.008, 0.008, 0.04, 0.06, 0.1, 1500.0, 30.0, 40.0]
# bounds, used only by BOBYQA (which needs a scaled box); the Nelder-Mead variants are unbounded like simpop
const LB = [0.40, 6.0, 15.0, 0.10, 0.10, 0.10, 0.8, 0.5, 7600.0, 1300.0, 1500.0]
const UB = [0.75, 18.0, 35.0, 0.17, 0.17, 0.45, 1.5, 1.5, 13700.0, 1750.0, 2000.0]

function set_design!(ac, x)
    ac.para[iaCL, ipclimb1+1:ipdescentn-1, 1]   .= x[1]
    ac.wing.layout.AR                            = x[2]
    ac.wing.layout.sweep                         = x[3]
    ac.htail.layout.sweep                        = x[3]
    ac.wing.inboard.cross_section.thickness_to_chord  = x[4]
    ac.wing.outboard.cross_section.thickness_to_chord = x[5]
    ac.wing.outboard.λ                           = max(x[6], 0.1)
    ac.para[iarcls, ipclimb1+1:ipdescentn-1, 1] .= x[7]
    ac.para[iarclt, ipclimb1+1:ipdescentn-1, 1] .= x[8]
    ac.para[iaalt, ipcruise1, 1]                 = x[9]
    ac.pare[ieTt4, ipcruise1:ipcruisen, 1]      .= x[10]
    ac.pare[ieTt4, ipstatic:iptakeoff, 1]       .= x[11]
    return nothing
end

function get_design(ac)
    return [ac.para[iaCL, ipcruise1, 1], ac.wing.layout.AR, ac.wing.layout.sweep,
            ac.wing.inboard.cross_section.thickness_to_chord,
            ac.wing.outboard.cross_section.thickness_to_chord, ac.wing.outboard.λ,
            ac.para[iarcls, ipcruise1, 1], ac.para[iarclt, ipcruise1, 1],
            ac.para[iaalt, ipcruise1, 1], ac.pare[ieTt4, ipcruise1, 1],
            ac.pare[ieTt4, ipstatic, 1]]
end

"Constraint values (positive = violated), matching the `con` array of fobj.f."
function constraints(ac)
    return (span = ac.wing.span / ac.wing.layout.max_span - 1.0,
            lBF  = ac.parm[imlBF, 1] / ac.parg[iglBFmax] - 1.0,
            fuel = ac.parg[igWfuel] / (ac.parg[igrWfmax] * ac.parg[igWfmax]) - 1.0,
            gtoc = 1.0 - ac.para[iagamV, ipclimbn, 1] / ac.parg[iggtocmin])
end

function penalized_objective(ac)
    Wpay = ac.parg[igWpay]
    f = ac.parg[igWfuel]
    # ad hoc penalties from fobj.f
    f += Wpay * max((1.5 - ac.wing.layout.sweep) / 100.0, 0.0)^2
    f += Wpay * max((ac.wing.inboard.λ - 1.0) / 0.25, 0.0)^2
    f += Wpay * max((0.15 - ac.wing.outboard.λ) / 0.25, 0.0)^2
    c = constraints(ac)
    f += 5.0 * Wpay * max(c.lBF, 0.0)^2
    f += 5.0 * Wpay * max(c.fuel, 0.0)^2
    f += 25.0 * Wpay * max(c.span, 0.0)^2
    f += 1.0 * Wpay * max(c.gtoc, 0.0)^2
    return f
end

# ── Evaluation wrapper: warm start, counting, trace, best tracking ─
mutable struct Problem
    ac::Any
    neval::Int
    warm::Bool
    t0::UInt64
    tsizing::Float64
    nfail::Int
    first_err::String
    best_f::Float64
    best_x::Vector{Float64}
    best_Wfuel::Float64
    best_con::NamedTuple
    trace::Vector{NTuple{4,Float64}}   # (eval, f, best_f, elapsed_s)
end

function Problem()
    ac = read_aircraft_model(INPUT_TOML)
    size_aircraft!(ac, printiter=false)   # baseline sizing; also initializes arrays
    return Problem(ac, 0, false, time_ns(), 0.0, 0, "", Inf, Float64[], NaN,
                   (span = NaN, lBF = NaN, fuel = NaN, gtoc = NaN), NTuple{4,Float64}[])
end

function fail!(p::Problem)
    # A failed sizing can leave NaNs in the aircraft state (weights, engine guesses) that make every
    # later evaluation fail too, so rebuild the aircraft from the input file; the next sizing is cold.
    p.nfail += 1
    p.warm = false
    p.ac = read_aircraft_model(INPUT_TOML)
    return 1e10
end

function evaluate!(p::Problem, x)
    p.neval += 1
    set_design!(p.ac, x)
    f = 1e10
    try
        ts = time_ns()
        size_aircraft!(p.ac, iter=50, printiter=false, initwgt=p.warm)
        p.tsizing += (time_ns() - ts) * 1e-9
        f = penalized_objective(p.ac)
        if isfinite(f)
            p.warm = true
        else
            f = fail!(p)
        end
    catch e
        isempty(p.first_err) && (p.first_err = sprint(showerror, e)[1:min(end, 300)])
        f = fail!(p)
    end
    if f < p.best_f
        p.best_f = f
        p.best_x = collect(Float64, x)
        if f < 1e10
            p.best_Wfuel = p.ac.parg[igWfuel]
            p.best_con = constraints(p.ac)
        end
    end
    push!(p.trace, (p.neval, f, p.best_f, (time_ns() - p.t0) * 1e-9))
    return f
end

# ── Optimizer runners: each returns (x, f, reason) ───────────────
function run_nlopt(alg::Symbol, p::Problem, x0; scaled=false, ftol_rel=1e-7, xtol_rel=0.0, maxeval=MAXEVAL)
    n = length(x0)
    opt = NLopt.Opt(alg, n)
    if scaled
        to_x(u) = LB .+ u .* (UB .- LB)
        opt.min_objective = (u, g) -> evaluate!(p, to_x(u))
        opt.lower_bounds = zeros(n)
        opt.upper_bounds = ones(n)
        opt.initial_step = fill(0.02, n)
        u0 = clamp.((x0 .- LB) ./ (UB .- LB), 0.0, 1.0)
        start = u0
    else
        opt.min_objective = (x, g) -> evaluate!(p, x)
        opt.initial_step = DX
        start = collect(Float64, x0)
    end
    opt.ftol_rel = ftol_rel
    opt.xtol_rel = xtol_rel
    opt.maxeval = maxeval
    (f, xs, ret) = NLopt.optimize(opt, start)
    return (scaled ? LB .+ xs .* (UB .- LB) : xs, f, string(ret))
end

struct AxisSimplexer <: Optim.Simplexer
    steps::Vector{Float64}
end
function Optim.simplexer(S::AxisSimplexer, x0::AbstractVector)
    n = length(x0)
    simplex = [collect(Float64, x0) for _ in 1:n+1]
    for i in 1:n
        simplex[i+1][i] += S.steps[i]
    end
    return simplex
end

"""
Optim.jl Nelder-Mead with simpop's coefficients (reflect 1, expand 1.4, contract 0.7, shrink 0.7).
Restart wrapper: after each convergence rebuild a fresh axis-step simplex at the best point and
continue until a restart improves the objective by less than `restart_rtol` (relative) or the
evaluation cap is hit.
"""
function run_optim_restart(p::Problem, x0; restart_rtol=1e-7, maxrestarts=30, maxeval=MAXEVAL)
    method = Optim.NelderMead(parameters = Optim.FixedParameters(1.0, 1.4, 0.7, 0.7),
                              initial_simplex = AxisSimplexer(DX))
    x = collect(Float64, x0)
    nrestart = 0
    reason = "max_restarts"
    for r in 1:maxrestarts
        remaining = maxeval - p.neval
        remaining <= 0 && (reason = "maxeval"; break)
        f_before = p.best_f
        Optim.optimize(z -> evaluate!(p, z), x, method,
                       Optim.Options(g_tol = 0.01, f_calls_limit = remaining, iterations = 10^6))
        nrestart = r
        if r > 1 && f_before - p.best_f < restart_rtol * abs(p.best_f)
            reason = "no_improvement_after_restart"
            break
        end
        x = copy(p.best_x)    # fresh simplex at the best point
    end
    return (p.best_x, p.best_f, "$(reason) (restarts=$(nrestart))")
end

const OPTIMIZERS = Dict(
    "NLopt_NelderMead" => (p, x0) -> run_nlopt(:LN_NELDERMEAD, p, x0),
    "NLopt_SBPLX"      => (p, x0) -> run_nlopt(:LN_SBPLX, p, x0),
    "NLopt_BOBYQA"     => (p, x0) -> run_nlopt(:LN_BOBYQA, p, x0; scaled=true, ftol_rel=1e-9, xtol_rel=1e-7),
    "Optim_NM_simpop_restart" => (p, x0) -> run_optim_restart(p, x0),
)
const OPT_ORDER = ["NLopt_NelderMead", "Optim_NM_simpop_restart", "NLopt_SBPLX", "NLopt_BOBYQA"]

# ── Starting points ──────────────────────────────────────────────
function start_points()
    base = get_design(read_aircraft_model(INPUT_TOML))
    starts = Dict(0 => base)
    rng = MersenneTwister(20261020)
    for k in 1:3
        starts[k] = base .+ START_SCALE .* (2 .* rand(rng, length(base)) .- 1)
        starts[k][6] = max(starts[k][6], 0.12)
    end
    return starts
end

# ── JSON / CSV output (no extra dependencies) ────────────────────
jval(v::AbstractString) = "\"" * replace(v, "\"" => "\\\"") * "\""
jval(v::Bool) = v ? "true" : "false"
jval(v::Integer) = string(v)
jval(v::Real) = isfinite(v) ? repr(Float64(v)) : "null"
jval(v::AbstractVector) = "[" * join(jval.(v), ",") * "]"
jval(v::NamedTuple) = "{" * join(("\"$k\":" * jval(x) for (k, x) in pairs(v)), ",") * "}"
jval(v::AbstractDict) = "{" * join(("\"$k\":" * jval(x) for (k, x) in v), ",") * "}"

function main()
    mkpath(OUTDIR)
    opts = split(get(ENV, "TASOPT_BENCH_OPTS", join(OPT_ORDER, ",")), ",")
    ids = parse.(Int, split(get(ENV, "TASOPT_BENCH_STARTS", "0,1,2,3"), ","))
    starts = start_points()

    # compile everything once so the first run's wall time is not dominated by JIT
    let p = Problem(); evaluate!(p, starts[0]); evaluate!(p, starts[0]) end

    records = String[]
    csv = ["run,optimizer,start,wall_s,sizing_s,nevals,nfail,f_final,Wfuel_N,span_viol,lBF_viol,fuel_viol,gtoc_viol,max_viol,feasible_1e-3,reason"]
    for oname in opts, id in ids
        runid = "$(oname)_s$(id)"
        p = Problem()               # fresh aircraft; first evaluation is cold, the rest warm
        p.t0 = time_ns()
        local x, f, reason
        wall = @elapsed begin
            try
                x, f, reason = OPTIMIZERS[oname](p, starts[id])
            catch e
                x, f, reason = p.best_x, p.best_f, "error: " * sprint(showerror, e)
            end
        end
        c = p.best_con
        maxv = maximum(max.(values(c), 0.0))
        feasible = maxv <= 1e-3
        @printf("%-26s start %d  wall %7.1fs  evals %5d  f* %.2f N  Wfuel %.2f  maxviol %.2e  %s\n",
                oname, id, wall, p.neval, p.best_f, p.best_Wfuel, maxv, reason)
        push!(records, "{" * join([
            "\"run\":" * jval(runid), "\"optimizer\":" * jval(oname), "\"start\":" * jval(id),
            "\"x0\":" * jval(starts[id]), "\"wall_s\":" * jval(wall), "\"sizing_s\":" * jval(p.tsizing),
            "\"nevals\":" * jval(p.neval), "\"nfail\":" * jval(p.nfail), "\"first_error\":" * jval(p.first_err),
            "\"f_best\":" * jval(p.best_f), "\"Wfuel_N\":" * jval(p.best_Wfuel),
            "\"x_best\":" * jval(p.best_x), "\"constraints\":" * jval(c),
            "\"max_violation\":" * jval(maxv), "\"feasible_1e-3\":" * jval(feasible),
            "\"reason\":" * jval(reason)], ",") * "}")
        push!(csv, join([runid, oname, id, wall, p.tsizing, p.neval, p.nfail, p.best_f, p.best_Wfuel,
                         c.span, c.lBF, c.fuel, c.gtoc, maxv, feasible, "\"$reason\""], ","))
        open(joinpath(OUTDIR, "evals_$(runid).csv"), "w") do io
            println(io, "eval,f,best_f,elapsed_s")
            for t in p.trace
                @printf(io, "%d,%.6f,%.6f,%.4f\n", t[1], t[2], t[3], t[4])
            end
        end
        # write incrementally so a crash keeps earlier results
        open(joinpath(OUTDIR, "runs.json"), "w") do io
            println(io, "{\"var_names\":", jval(VAR_NAMES), ",\"runs\":[\n", join(records, ",\n"), "\n]}")
        end
        open(joinpath(OUTDIR, "runs_summary.csv"), "w") do io
            println(io, join(csv, "\n"))
        end
    end
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
