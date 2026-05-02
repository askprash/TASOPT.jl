"""
engine_harness.jl — Engine-standalone runners (single-point and multi-point sweep).

Provides:
- `run_engine_design_point`: single-point design sizing runner — accepts an
  aircraft object, initialises ambient conditions, runs `tfsize!` via
  `tfwrap!`, and returns an `EngineState`.
- `run_engine_sweep`: multi-point off-design runner — iterates the engine
  across a range of mission points, writes results into `ac.missions[imission]`,
  and returns that `Mission{Float64}`.
- `write_sweep_csv`: serialise a `Mission` sweep to CSV.
- `write_sweep_toml`: serialise a `Mission` sweep to TOML regression baseline.
- `SweepResult`: **deprecated** tabular container (use `Mission{T}` instead).
"""

# ---------------------------------------------------------------------------
# Legacy ducted-fan compatibility bridge — pare-slice → FlowStation helpers
#
# These helpers exist only because ducted_fan_harness.jl still stages engine
# state through a pare-style vector before constructing FlowStation objects.
# They are NOT used by the turbofan path (tfcalc!/tfoper!/tfsize!) and should
# be deleted once ducted_fan_harness.jl is migrated to EngineState.
# ---------------------------------------------------------------------------

"""
    _fill_total!(fs, pare, iTt, iht, ipt, icpt, iRt)

Write the five total-state scalars from a `pare` slice into `fs`.

**Ducted-fan bridge only.** Called by `ducted_fan_harness.jl`; not used by
the turbofan Newton solver.
"""
@inline function _fill_total!(fs::FlowStation, pare,
                               iTt::Int, iht::Int, ipt::Int, icpt::Int, iRt::Int)
    fs.Tt  = pare[iTt]
    fs.ht  = pare[iht]
    fs.pt  = pare[ipt]
    fs.cpt = pare[icpt]
    fs.Rt  = pare[iRt]
end

"""
    _fill_static!(fs, pare, ip, iT, iR, icp, iu)

Write the five static-state scalars (p, T, R, cp, u) from a `pare` slice
into `fs`.  Static enthalpy `hs` and entropy complement `ss` are not in
`pare`; they remain at their prior value.

**Ducted-fan bridge only.** Called by `ducted_fan_harness.jl`; not used by
the turbofan Newton solver.
"""
@inline function _fill_static!(fs::FlowStation, pare,
                                ip::Int, iT::Int, iR::Int, icp::Int, iu::Int)
    fs.ps  = pare[ip]
    fs.Ts  = pare[iT]
    fs.Rs  = pare[iR]
    fs.cps = pare[icp]
    fs.u   = pare[iu]
end

# engine_state_to_pare_vec was deleted (tasopt-j9l.45.14.7.7.4): ie* index constants
# removed from index.inc; NPSS integration paths in odperformance.jl that called
# this function have been stubbed with an empty matrix placeholder.

# ---------------------------------------------------------------------------
# run_engine_design_point
# ---------------------------------------------------------------------------

"""
    run_engine_design_point(ac; imission=1, ip=ipcruise1) -> EngineState

Run the turbofan design-point thermodynamic sizing for one mission point and
return the result as a populated `EngineState`.

The full structural/aerodynamic aircraft weight-convergence loop is **not**
run — only the engine thermodynamics (`tfsize!` via `tfcalc!`/`tfwrap!`).
Ambient conditions are computed from the ISA standard atmosphere at the
design altitude and Mach number declared in the input file.

## Arguments
- `ac`: an `aircraft` object loaded with `load_default_model()` or
  `read_aircraft_model()`.  The aircraft should be sized with `size_aircraft!`
  first so that the design thrust is stored in the typed EngineState.
- `imission::Int`: mission index (default `1`).
- `ip::Int`: mission-point index (default `ipcruise1`).

## Returns
An `EngineState{Float64}` populated from the converged design-point typed state
at `ip`.  Station fields not stored in the typed state (total state for 19c,
25c, 4a, 49c; `hs` and `ss` for all stations) have value zero.

## Notes
For a fully-sized aircraft, `run_engine_design_point` reproduces the same
thermodynamic result as the final iteration of `size_aircraft!`.

## Example
```julia
using TASOPT
ac = TASOPT.load_default_model()
size_aircraft!(ac; printiter=false)
eng = TASOPT.engine.run_engine_design_point(ac)
TASOPT.engine.dump_stations(eng)
```
"""
function run_engine_design_point(ac; imission::Int=1, ip::Int=ipcruise1)

    # -----------------------------------------------------------------------
    # Initialise atmosphere ISA offset — mirrors step 1 of _size_aircraft!
    # so that the standard-atmosphere call uses the correct ΔT.
    # -----------------------------------------------------------------------
    altTO = ac.parm[imaltTO]
    T_std = atmos(altTO).T
    ac.parm[imDeltaTatm] = ac.parm[imT0TO] - T_std

    # -----------------------------------------------------------------------
    # Set ambient condition scalars (T0, p0, a0, M0, u0, ρ0, μ0) from design altitude + Mach.
    # This is the same logic as set_ambient_conditions! in size_aircraft.jl,
    # reproduced here so the harness does not depend on that unexported helper.
    # -----------------------------------------------------------------------
    ΔTatmos = ac.parm[imDeltaTatm]
    alt_m   = ac.para[iaalt, ip, imission]
    as      = atmos(alt_m, ΔTatmos)
    Mach    = ac.para[iaMach, ip, imission]

    ac.para[iaReunit, ip, imission] = Mach * as.a * as.ρ / as.μ

    # -----------------------------------------------------------------------
    # Run engine design-point sizing (no aircraft weight / aero loop)
    # -----------------------------------------------------------------------
    tfwrap!(ac, "design", imission, ip, true)

    return deepcopy(ac.missions[imission].points[ip].engine)
end

# ---------------------------------------------------------------------------
# SweepResult
# ---------------------------------------------------------------------------

"""
    SweepResult{T<:AbstractFloat}

!!! warning "Deprecated"
    `SweepResult` is deprecated.  [`run_engine_sweep`](@ref) now returns a
    `Mission{T}` directly and writes results into `ac.missions[imission]`.
    Use `mission.points[ip].engine` to access per-point engine state.
    `SweepResult` will be removed in a future release.

Tabular container for the output of [`run_engine_sweep`](@ref).  Holds one
[`EngineState`](@ref) per mission point together with the key engine
performance scalars extracted from the typed `EngineState`.

## Fields

| Field       | Size  | Unit    | Description                              |
|:------------|:------|:--------|:-----------------------------------------|
| `ip_indices`| `n`   | —       | Mission-point indices (`ip` values)       |
| `ip_labels` | `n`   | —       | Human-readable labels (e.g. `"cruise1"`) |
| `alt`       | `n`   | m       | Altitude at each point                   |
| `Mach`      | `n`   | —       | Mach number at each point                |
| `engines`   | `n`   | —       | `EngineState` per point                  |
| `Fe`        | `n`   | N       | Net thrust per engine                    |
| `TSFC`      | `n`   | kg/N/s  | Thrust-specific fuel consumption         |
| `BPR`       | `n`   | —       | Bypass ratio                             |
| `mcore`     | `n`   | kg/s    | Core mass flow rate (per engine)         |
| `mdotf`     | `n`   | kg/s    | Fuel mass flow rate (total, all engines)                       |
"""
struct SweepResult{T<:AbstractFloat}
    ip_indices ::Vector{Int}
    ip_labels  ::Vector{String}
    alt        ::Vector{T}
    Mach       ::Vector{T}
    engines    ::Vector{EngineState{T}}
    Fe         ::Vector{T}
    TSFC       ::Vector{T}
    BPR        ::Vector{T}
    mcore      ::Vector{T}
    mdotf      ::Vector{T}
end

# ---------------------------------------------------------------------------
# run_engine_sweep
# ---------------------------------------------------------------------------

"""
    run_engine_sweep(ac; imission=1, ip_range=ipstatic:ipdescentn,
                    initializes_engine=false) -> Mission{Float64}

Run the turbofan off-design performance routine across every mission point in
`ip_range`, write the converged engine state into `ac.missions[imission]`, and
return that `Mission{Float64}` directly.

The aircraft must already be sized with `size_aircraft!` (or equivalent)
before calling this function.  In particular:
- The required thrust at each mission point must be stored in the typed
  `EngineState` (set by `size_aircraft!` ↔ `fly_mission!`).
- Ambient conditions in typed state (T0, p0, a0, M0, …) must be populated.
- `ac.missions[imission]` must be pre-allocated (done by `read_input.jl`).

For climb and ground-roll points (`ipstatic:ipclimbn`) the engine is operated
at maximum turbine-inlet temperature (`CalcMode.FixedTt4OffDes`); for all
other points it operates at the thrust target (`CalcMode.FixedFeOffDes`).

## Arguments
- `ac`: a sized `aircraft` object.
- `imission::Int`: mission index (default `1`).
- `ip_range`: iterable of mission-point indices (default `ipstatic:ipdescentn`,
  i.e., all 16 regular mission points).
- `initializes_engine::Bool`: passed to `tfwrap!` for each point.
  Default `false` — use the converged `EngineState` as the Newton initial guess.
  Pass `true` to reinitialise the Newton iteration from scratch at every point
  (slower, but independent of the prior converged state).

## Returns
`ac.missions[imission]::Mission{Float64}` — the mission whose `points[ip].engine`
fields have been populated for every `ip` in `ip_range`.

## Example
```julia
using TASOPT
ac = TASOPT.load_default_model()
size_aircraft!(ac; printiter=false)
mission = TASOPT.engine.run_engine_sweep(ac)
TASOPT.engine.write_sweep_csv("engine_sweep.csv", mission, ipstatic:ipdescentn, ac)
```
"""
function run_engine_sweep(ac;
                          imission::Int        = 1,
                          ip_range             = ipstatic:ipdescentn,
                          initializes_engine::Bool = false)

    mission = ac.missions[imission]

    for ip in ip_range
        # Run engine off-design at this mission point
        tfwrap!(ac, "off_design", imission, ip, initializes_engine)
    end

    return mission
end

# ---------------------------------------------------------------------------
# write_sweep_csv
# ---------------------------------------------------------------------------

const _SWEEP_CSV_STATIONS = (
    ("0",  :st0),  ("2",  :st2),  ("3",   :st3),
    ("4",  :st4),  ("41", :st41), ("45",  :st45),
    ("5",  :st5),  ("8",  :st8),  ("18",  :st18))

"""
    write_sweep_csv(io::IO, mission::Mission, ip_range, ac; imission=1)
    write_sweep_csv(path::AbstractString, mission::Mission, ip_range, ac; imission=1)

Write a [`Mission`](@ref) sweep to CSV format for the mission points in `ip_range`.

Each row corresponds to one mission point.  Columns:
- Metadata: `ip`, `label`, `alt_m`, `Mach`
- Engine performance: `Fe_N`, `TSFC_kg_Ns`, `BPR`, `mcore_kg_s`, `mfuel_kg_s`
- Station totals at key stations (0, 2, 3, 4, 41, 45, 5, 8, 18):
  `Tt<N>_K`, `pt<N>_Pa`, `ht<N>_J_kg`
- Station exit velocities at nozzle stations 8 and 18:
  `u8_m_s`, `u18_m_s`

The `path` overload opens the file, writes, and closes it.
"""
function write_sweep_csv(io::IO, mission, ip_range, ac; imission::Int=1)
    # ----- Header -----
    header_parts = [
        "ip", "label", "alt_m", "Mach",
        "Fe_N", "TSFC_kg_Ns", "BPR", "mcore_kg_s", "mfuel_kg_s",
    ]
    for (sname, _) in _SWEEP_CSV_STATIONS
        push!(header_parts, "Tt$(sname)_K", "pt$(sname)_Pa", "ht$(sname)_J_kg")
    end
    push!(header_parts, "u8_m_s", "u18_m_s")
    println(io, join(header_parts, ","))

    # ----- Rows -----
    for ip in ip_range
        eng = mission.points[ip].engine
        lbl = (ip <= length(ip_labels)) ? ip_labels[ip] : string(ip)
        row = [
            string(ip),
            lbl,
            @sprintf("%.2f", ac.para[iaalt,  ip, imission]),
            @sprintf("%.6g", ac.para[iaMach, ip, imission]),
            @sprintf("%.6g", eng.Fe),
            @sprintf("%.6g", eng.TSFC),
            @sprintf("%.6g", eng.BPR),
            @sprintf("%.6g", eng.st2.mdot),
            @sprintf("%.6g", eng.mfuel),
        ]
        for (_, stfield) in _SWEEP_CSV_STATIONS
            fs = getfield(eng, stfield)
            push!(row,
                  @sprintf("%.6g", fs.Tt),
                  @sprintf("%.6g", fs.pt),
                  @sprintf("%.6g", fs.ht))
        end
        push!(row,
              @sprintf("%.6g", eng.st8.u),
              @sprintf("%.6g", eng.st18.u))
        println(io, join(row, ","))
    end
    return nothing
end

function write_sweep_csv(path::AbstractString, mission, ip_range, ac; imission::Int=1)
    open(path, "w") do io
        write_sweep_csv(io, mission, ip_range, ac; imission=imission)
    end
    return nothing
end

"""
    write_sweep_csv(io::IO, result::SweepResult)
    write_sweep_csv(path::AbstractString, result::SweepResult)

!!! warning "Deprecated"
    This overload is deprecated.  Use `write_sweep_csv(io, mission, ip_range, ac)`
    with the `Mission{T}` returned by [`run_engine_sweep`](@ref).

Write a [`SweepResult`](@ref) to CSV format.
"""
function write_sweep_csv(io::IO, result::SweepResult)
    # ----- Header -----
    header_parts = [
        "ip", "label", "alt_m", "Mach",
        "Fe_N", "TSFC_kg_Ns", "BPR", "mcore_kg_s", "mfuel_kg_s",
    ]
    for (sname, _) in _SWEEP_CSV_STATIONS
        push!(header_parts, "Tt$(sname)_K", "pt$(sname)_Pa", "ht$(sname)_J_kg")
    end
    push!(header_parts, "u8_m_s", "u18_m_s")
    println(io, join(header_parts, ","))

    # ----- Rows -----
    n = length(result.ip_indices)
    for k in 1:n
        eng = result.engines[k]
        row = [
            string(result.ip_indices[k]),
            result.ip_labels[k],
            @sprintf("%.2f", result.alt[k]),
            @sprintf("%.6g", result.Mach[k]),
            @sprintf("%.6g", result.Fe[k]),
            @sprintf("%.6g", result.TSFC[k]),
            @sprintf("%.6g", result.BPR[k]),
            @sprintf("%.6g", result.mcore[k]),
            @sprintf("%.6g", result.mdotf[k]),
        ]
        for (_, stfield) in _SWEEP_CSV_STATIONS
            fs = getfield(eng, stfield)
            push!(row,
                  @sprintf("%.6g", fs.Tt),
                  @sprintf("%.6g", fs.pt),
                  @sprintf("%.6g", fs.ht))
        end
        push!(row,
              @sprintf("%.6g", eng.st8.u),
              @sprintf("%.6g", eng.st18.u))
        println(io, join(row, ","))
    end
    return nothing
end

function write_sweep_csv(path::AbstractString, result::SweepResult)
    open(path, "w") do io
        write_sweep_csv(io, result)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# write_sweep_toml — diffable full-fidelity station baseline
# ---------------------------------------------------------------------------

import TOML

# Stations serialised in physical flow-path order (mirrors _STATION_DUMP_ORDER).
const _TOML_STATION_ORDER = (
    ("0",    "Freestream",     :st0),
    ("2",    "FanFaceFan",     :st2),
    ("12",   "FanFaceOuter",   :st12),
    ("2a",   "FanFaceLPC",     :st2a),
    ("2ac",  "PreCoolerOut",   :st2ac),
    ("13",   "FanExit",        :st13),
    ("25",   "LPCExit",        :st25),
    ("25c",  "InterCoolerOut", :st25c),
    ("3",    "HPCExit",        :st3),
    ("4",    "CombustorExit",  :st4),
    ("4a",   "CoolMixInlet",   :st4a),
    ("41",   "TurbineInlet",   :st41),
    ("45",   "HPTExit",        :st45),
    ("5",    "LPTExit",        :st5),
    ("5c",   "RegenCoolerOut", :st5c),
    ("8",    "CoreNozzle",     :st8),
    ("9",    "CoreNozzleExit", :st9),
    ("18",   "FanNozzle",      :st18),
    ("19",   "FanNozzleExit",  :st19),
    ("25off","OfftakeDisch",   :st25off),
)

"""
    _station_to_dict(fs::FlowStation) -> Dict{String,Float64}

Serialise the populated scalar fields of a `FlowStation` to a plain
`Dict{String,Float64}`.  Fields not written by `tfcalc!`
(hs, ss, st, alpha) are excluded to keep the baseline compact and to
avoid noisy diffs on zero-valued entries.
"""
function _station_to_dict(fs::FlowStation)
    Dict{String,Float64}(
        "Tt"   => Float64(fs.Tt),
        "ht"   => Float64(fs.ht),
        "pt"   => Float64(fs.pt),
        "cpt"  => Float64(fs.cpt),
        "Rt"   => Float64(fs.Rt),
        "Ts"   => Float64(fs.Ts),
        "ps"   => Float64(fs.ps),
        "cps"  => Float64(fs.cps),
        "Rs"   => Float64(fs.Rs),
        "u"    => Float64(fs.u),
        "A"    => Float64(fs.A),
        "mdot" => Float64(fs.mdot),
    )
end

"""
    write_sweep_toml(io::IO, mission::Mission, ip_range, ac; imission=1)
    write_sweep_toml(path::AbstractString, mission::Mission, ip_range, ac; imission=1)

Write a [`Mission`](@ref) sweep to TOML format suitable for use as a
regression baseline fixture.

Each `[[points]]` table contains:
- Metadata: `ip`, `label`, `alt_m`, `Mach`.
- Engine performance: `Fe_N`, `TSFC_kg_Ns`, `BPR`, `mcore_kg_s`,
  `mfuel_kg_s`, `M0`, `T0_K`, `p0_Pa`, `a0_m_s`.
- Station subtables (e.g. `[points.stations.st0]`) with twelve scalar
  fields: `Tt`, `ht`, `pt`, `cpt`, `Rt`, `Ts`, `ps`, `cps`, `Rs`,
  `u`, `A`, `mdot`.

The TOML format is line-oriented and deterministic, making it suitable
for `git diff` review when baseline values change.
"""
function write_sweep_toml(io::IO, mission, ip_range, ac; imission::Int=1)
    ips = collect(Int, ip_range)
    n = length(ips)
    points = Vector{Dict{String,Any}}(undef, n)
    for (k, ip) in enumerate(ips)
        eng = mission.points[ip].engine
        lbl = (ip <= length(ip_labels)) ? ip_labels[ip] : string(ip)
        stations = Dict{String,Any}()
        for (_, _, stfld) in _TOML_STATION_ORDER
            stations[String(stfld)] = _station_to_dict(getfield(eng, stfld))
        end
        points[k] = Dict{String,Any}(
            "ip"          => ip,
            "label"       => lbl,
            "alt_m"       => Float64(ac.para[iaalt,  ip, imission]),
            "Mach"        => Float64(ac.para[iaMach, ip, imission]),
            "Fe_N"        => Float64(eng.Fe),
            "TSFC_kg_Ns"  => Float64(eng.TSFC),
            "BPR"         => Float64(eng.BPR),
            "mcore_kg_s"  => Float64(eng.st2.mdot),
            "mfuel_kg_s"  => Float64(eng.mfuel),
            "M0"          => Float64(eng.M0),
            "T0_K"        => Float64(eng.T0),
            "p0_Pa"       => Float64(eng.p0),
            "a0_m_s"      => Float64(eng.a0),
            "stations"    => stations,
        )
    end
    data = Dict{String,Any}(
        "metadata" => Dict{String,Any}(
            "n_points" => n,
            "aircraft" => "default",
            "description" => "Engine sweep baseline — default TASOPT aircraft.",
        ),
        "points" => points,
    )
    TOML.print(io, data)
    return nothing
end

function write_sweep_toml(path::AbstractString, mission, ip_range, ac; imission::Int=1)
    open(path, "w") do io
        write_sweep_toml(io, mission, ip_range, ac; imission=imission)
    end
    return nothing
end

"""
    write_sweep_toml(io::IO, result::SweepResult)
    write_sweep_toml(path::AbstractString, result::SweepResult)

!!! warning "Deprecated"
    This overload is deprecated.  Use `write_sweep_toml(io, mission, ip_range, ac)`
    with the `Mission{T}` returned by [`run_engine_sweep`](@ref).

Write a [`SweepResult`](@ref) to TOML format.
"""
function write_sweep_toml(io::IO, result::SweepResult)
    n = length(result.ip_indices)
    points = Vector{Dict{String,Any}}(undef, n)
    for k in 1:n
        eng = result.engines[k]
        stations = Dict{String,Any}()
        for (_, _, stfld) in _TOML_STATION_ORDER
            stations[String(stfld)] = _station_to_dict(getfield(eng, stfld))
        end
        points[k] = Dict{String,Any}(
            "ip"          => result.ip_indices[k],
            "label"       => result.ip_labels[k],
            "alt_m"       => Float64(result.alt[k]),
            "Mach"        => Float64(result.Mach[k]),
            "Fe_N"        => Float64(result.Fe[k]),
            "TSFC_kg_Ns"  => Float64(result.TSFC[k]),
            "BPR"         => Float64(result.BPR[k]),
            "mcore_kg_s"  => Float64(result.mcore[k]),
            "mfuel_kg_s"  => Float64(result.mdotf[k]),
            "M0"          => Float64(eng.M0),
            "T0_K"        => Float64(eng.T0),
            "p0_Pa"       => Float64(eng.p0),
            "a0_m_s"      => Float64(eng.a0),
            "stations"    => stations,
        )
    end
    data = Dict{String,Any}(
        "metadata" => Dict{String,Any}(
            "n_points" => n,
            "aircraft" => "default",
            "description" => "Engine sweep baseline — default TASOPT aircraft.",
        ),
        "points" => points,
    )
    TOML.print(io, data)
    return nothing
end

function write_sweep_toml(path::AbstractString, result::SweepResult)
    open(path, "w") do io
        write_sweep_toml(io, result)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# regenerate_engine_baseline
# ---------------------------------------------------------------------------

"""
    const ENGINE_BASELINE_PATH

Default path for the engine-sweep regression baseline fixture.
Points to `test/fixtures/engine_sweep_baseline.toml` relative to the
repository root.
"""
const ENGINE_BASELINE_PATH = joinpath(
    @__DIR__, "..", "..", "..", "test", "fixtures", "engine_sweep_baseline.toml")

"""
    regenerate_engine_baseline(ac; path=ENGINE_BASELINE_PATH)

**Auditable action** — regenerates the engine-sweep regression baseline
from a sized aircraft object and writes it to `path`.

Run the full off-design sweep on `ac` (which must already be sized with
`size_aircraft!`) and serialise the result to TOML via `write_sweep_toml`.

## ⚠  Warning

Regenerating the baseline changes the numerical reference used by the
engine regression test (`tasopt-j9l.5`).  Only call this function when a
physics or model change is intentional and the new values have been
verified.  Commit the updated baseline with an explanatory git message —
reviewers will see the numerical diff in the pull request.

## Usage
```julia
using TASOPT
ac = TASOPT.load_default_model()
TASOPT.size_aircraft!(ac; printiter=false)

# Regenerate to the default test-fixture path
TASOPT.engine.regenerate_engine_baseline(ac)

# Or write to a custom path for inspection first
TASOPT.engine.regenerate_engine_baseline(ac; path="/tmp/my_baseline.toml")
```
"""
function regenerate_engine_baseline(ac; path::AbstractString=ENGINE_BASELINE_PATH)
    @warn """
    *** REGENERATING ENGINE SWEEP BASELINE ***
    Output path: $(abspath(path))

    This rewrites the numerical reference for engine regression tests.
    Only proceed when the change is intentional and has been verified.
    Commit the updated baseline file with an explanatory git message so
    reviewers can audit the numerical diff.
    """
    mkpath(dirname(abspath(path)))
    mission = run_engine_sweep(ac)
    write_sweep_toml(path, mission, ipstatic:ipdescentn, ac)
    n = ipdescentn - ipstatic + 1
    @info "Baseline written: $(abspath(path))  ($n mission points)"
    return abspath(path)
end
