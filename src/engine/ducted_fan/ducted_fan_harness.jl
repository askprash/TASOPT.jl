"""
ducted_fan_harness.jl — Ducted-fan engine-standalone runners.

Provides:
- `DuctedFanState`: typed container for one ducted-fan operating point.
- `pare_to_ducted_fan_state!`: populate DuctedFanState from a typed EngineState.
- `run_ducted_fan_design_point`: single-point sizing runner.
- `run_ducted_fan_sweep`: multi-point off-design sweep.
- `write_ducted_fan_sweep_toml`: serialize sweep to TOML regression baseline.
- `DUCTED_FAN_BASELINE_PATH`: default fixture path.
- `regenerate_ducted_fan_baseline`: convenience regenerator.
"""

import TOML

# ---------------------------------------------------------------------------
# DuctedFanState
# ---------------------------------------------------------------------------

"""
    DuctedFanState{T<:AbstractFloat}

Mutable container for the complete aerothermodynamic state of a TASOPT
ducted-fan (electric or fuel-cell driven) at a single operating point.

## Ambient / flight condition scalars

| Field | Unit  | Description                         |
|:------|:------|:------------------------------------|
| `M0`  | —     | Freestream Mach number              |
| `T0`  | K     | Freestream static temperature       |
| `p0`  | Pa    | Freestream static pressure          |
| `a0`  | m/s   | Freestream speed of sound           |
| `u0`  | m/s   | Freestream flight speed             |

## Performance scalars

| Field  | Unit      | Description                              |
|:-------|:----------|:-----------------------------------------|
| `Fe`   | N         | Net thrust per engine                    |
| `Fsp`  | N/(kg/s)  | Specific thrust                          |
| `TSEC` | W/N       | Thrust-specific electric consumption     |
| `TSFC` | kg/N/s    | Thrust-specific fuel consumption         |
| `Pfan` | W         | Fan shaft power                          |
| `mfan` | kg/s      | Fan mass flow                            |
| `etaf` | —         | Fan adiabatic efficiency                 |
| `epf`  | —         | Fan polytropic efficiency                |

## Fan map operating point

| Field | Unit          | Description                    |
|:------|:--------------|:-------------------------------|
| `pif` | —             | Fan pressure ratio             |
| `mbf` | kg/s corrected| Fan corrected mass flow        |
| `Nbf` | —             | Fan corrected spool speed      |

## Design anchors (frozen at sizing point)

| Field  | Unit          | Description                    |
|:-------|:--------------|:-------------------------------|
| `pi_fan_des` | —             | Design fan pressure ratio      |
| `mb_fan_des` | kg/s corrected| Design fan corrected flow      |
| `Nb_fan_des` | —             | Design fan corrected speed     |
| `A2`   | m²            | Fan-face area                  |
| `A7`   | m²            | Fan-nozzle throat area         |

## Flow stations

| Field  | TASOPT # | Description                             |
|:-------|:--------:|:----------------------------------------|
| `st0`  | 0        | Freestream (total + u)                  |
| `st18` | 18       | Fan face, outer (total)                 |
| `st2`  | 2        | Fan face (total + static + A2 + mfan)  |
| `st21` | 21       | Fan exit (total)                        |
| `st7`  | 7        | Fan nozzle throat (total + static + A7) |
| `st8`  | 8        | Fan nozzle exit (static + A8 only)      |

## Constructors

```julia
DuctedFanState{T}()   # all fields zeroed, numeric type T
DuctedFanState()      # Float64 default
```
"""
mutable struct DuctedFanState{T<:AbstractFloat}
    # -----------------------------------------------------------------------
    # Ambient / flight condition
    # -----------------------------------------------------------------------
    M0 ::T   # freestream Mach number              [—]
    T0 ::T   # freestream static temperature       [K]
    p0 ::T   # freestream static pressure          [Pa]
    a0 ::T   # freestream speed of sound           [m/s]
    u0 ::T   # freestream flight speed             [m/s]

    # -----------------------------------------------------------------------
    # Performance scalars
    # -----------------------------------------------------------------------
    Fe   ::T   # net thrust per engine              [N]
    Fsp  ::T   # specific thrust                    [N/(kg/s)]
    TSEC ::T   # thrust-specific electric consumption [W/N]
    TSFC ::T   # thrust-specific fuel consumption   [kg/N/s]
    Pfan ::T   # fan shaft power                    [W]
    mfan ::T   # fan mass flow                      [kg/s]
    etaf ::T   # fan adiabatic efficiency           [—]
    epf  ::T   # fan polytropic efficiency          [—]

    # -----------------------------------------------------------------------
    # Fan map operating point
    # -----------------------------------------------------------------------
    pif  ::T   # fan pressure ratio                 [—]
    mbf  ::T   # fan corrected mass flow            [kg/s corrected]
    Nbf  ::T   # fan corrected spool speed          [—]

    # -----------------------------------------------------------------------
    # Design anchors (frozen at sizing point)
    # -----------------------------------------------------------------------
    pi_fan_des ::T   # design fan pressure ratio          [—]
    mb_fan_des ::T   # design fan corrected flow          [kg/s corrected]
    Nb_fan_des ::T   # design fan corrected speed         [—]
    A2   ::T   # fan-face area                      [m²]
    A7   ::T   # fan-nozzle throat area             [m²]

    # -----------------------------------------------------------------------
    # Flow stations
    # -----------------------------------------------------------------------
    st0  ::FlowStation{T}   # Freestream (total + u)
    st18 ::FlowStation{T}   # FanFaceOuter (total)
    st2  ::FlowStation{T}   # FanFace (total + static + A2 + mfan)
    st21 ::FlowStation{T}   # FanExit (total)
    st7  ::FlowStation{T}   # FanNozzle (total + static + A7)
    st8  ::FlowStation{T}   # FanNozzleExit (static + A8 only)
end

function DuctedFanState{T}() where {T<:AbstractFloat}
    DuctedFanState{T}(
        # ambient
        zero(T), zero(T), zero(T), zero(T), zero(T),
        # performance
        zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T),
        # fan map
        zero(T), zero(T), zero(T),
        # design anchors
        zero(T), zero(T), zero(T), zero(T), zero(T),
        # stations
        FlowStation{T}(), FlowStation{T}(), FlowStation{T}(),
        FlowStation{T}(), FlowStation{T}(), FlowStation{T}(),
    )
end

DuctedFanState() = DuctedFanState{Float64}()

# ---------------------------------------------------------------------------
# pare_to_ducted_fan_state!
# ---------------------------------------------------------------------------

"""
    pare_to_ducted_fan_state!(state, eng_ip::EngineState) -> state

Populate a `DuctedFanState` from a typed `EngineState`
(`ac.missions[imission].points[ip].engine`) after a successful `ductedfancalc!`
call.  Reads exclusively from the typed surface; no bare `pare` access.

Used internally by `run_ducted_fan_design_point` and `run_ducted_fan_sweep`.
"""
function pare_to_ducted_fan_state!(state::DuctedFanState, eng_ip::EngineState)
    # -----------------------------------------------------------------------
    # Ambient
    # -----------------------------------------------------------------------
    state.M0 = eng_ip.M0
    state.T0 = eng_ip.T0
    state.p0 = eng_ip.p0
    state.a0 = eng_ip.a0
    state.u0 = eng_ip.st0.u

    # -----------------------------------------------------------------------
    # Performance
    # -----------------------------------------------------------------------
    state.Fe   = eng_ip.Fe
    state.Fsp  = eng_ip.Fsp
    state.TSEC = eng_ip.TSEC
    state.TSFC = eng_ip.TSFC
    state.Pfan = eng_ip.Pfan
    state.mfan = eng_ip.mfan
    state.etaf = eng_ip.etaf
    state.epf  = eng_ip.epf

    # -----------------------------------------------------------------------
    # Fan map operating point
    # -----------------------------------------------------------------------
    state.pif  = eng_ip.pif
    state.mbf  = eng_ip.mbf
    state.Nbf  = eng_ip.Nbf

    # -----------------------------------------------------------------------
    # Design anchors
    # -----------------------------------------------------------------------
    state.pi_fan_des = eng_ip.design.pi_fan_des
    state.mb_fan_des = eng_ip.design.mb_fan_des
    state.Nb_fan_des = eng_ip.design.Nb_fan_des
    state.A2   = eng_ip.design.A2
    state.A7   = eng_ip.design.A7

    # -----------------------------------------------------------------------
    # Station 0 — freestream (total + u)
    # -----------------------------------------------------------------------
    let st = eng_ip.st0
        state.st0.Tt  = st.Tt
        state.st0.ht  = st.ht
        state.st0.pt  = st.pt
        state.st0.cpt = st.cpt
        state.st0.Rt  = st.Rt
        state.st0.u   = st.u
    end

    # -----------------------------------------------------------------------
    # Station 18 — FanFaceOuter (total only)
    # -----------------------------------------------------------------------
    let st = eng_ip.st18
        state.st18.Tt  = st.Tt
        state.st18.ht  = st.ht
        state.st18.pt  = st.pt
        state.st18.cpt = st.cpt
        state.st18.Rt  = st.Rt
    end

    # -----------------------------------------------------------------------
    # Station 2 — FanFace (total + static + A2 + mfan)
    # -----------------------------------------------------------------------
    let st = eng_ip.st2
        state.st2.Tt  = st.Tt
        state.st2.ht  = st.ht
        state.st2.pt  = st.pt
        state.st2.cpt = st.cpt
        state.st2.Rt  = st.Rt
        state.st2.ps  = st.ps
        state.st2.Ts  = st.Ts
        state.st2.Rs  = st.Rs
        state.st2.cps = st.cps
        state.st2.u   = st.u
    end
    state.st2.A    = eng_ip.design.A2
    state.st2.mdot = eng_ip.mfan

    # -----------------------------------------------------------------------
    # Station 21 — FanExit (total only)
    # -----------------------------------------------------------------------
    let st = eng_ip.st21
        state.st21.Tt  = st.Tt
        state.st21.ht  = st.ht
        state.st21.pt  = st.pt
        state.st21.cpt = st.cpt
        state.st21.Rt  = st.Rt
    end

    # -----------------------------------------------------------------------
    # Station 7 — FanNozzle (total + static + A7)
    # -----------------------------------------------------------------------
    let st = eng_ip.st7
        state.st7.Tt  = st.Tt
        state.st7.ht  = st.ht
        state.st7.pt  = st.pt
        state.st7.cpt = st.cpt
        state.st7.Rt  = st.Rt
        state.st7.ps  = st.ps
        state.st7.Ts  = st.Ts
        state.st7.Rs  = st.Rs
        state.st7.cps = st.cps
        state.st7.u   = st.u
    end
    state.st7.A = eng_ip.design.A7

    # -----------------------------------------------------------------------
    # Station 8 — FanNozzleExit (static + A8 only; no total in EngineState)
    # -----------------------------------------------------------------------
    let st = eng_ip.st8
        state.st8.ps  = st.ps
        state.st8.Ts  = st.Ts
        state.st8.Rs  = st.Rs
        state.st8.cps = st.cps
        state.st8.u   = st.u
        state.st8.A   = st.A
    end

    return state
end

# ---------------------------------------------------------------------------
# run_ducted_fan_design_point
# ---------------------------------------------------------------------------

"""
    run_ducted_fan_design_point(ac; imission=1, ip=ipcruise1) -> DuctedFanState

Run the ducted-fan design-point sizing for one mission point and return the
result as a populated `DuctedFanState`.

Only the engine thermodynamics (`ductedfancalc!` via `ac.engine.enginecalc!`)
are run — the full aircraft weight-convergence loop is not.  Ambient conditions
are computed from the ISA standard atmosphere at the design altitude and Mach
number stored in `para[iaalt, ip, imission]` and `para[iaMach, ip, imission]`.

## Arguments
- `ac`: an `aircraft` object whose `engine` is a ducted-fan model (e.g.
  `fuel_cell_with_ducted_fan`).  The design thrust must be stored in the typed
  `DuctedFanState` before calling this function.
- `imission::Int`: mission index (default `1`).
- `ip::Int`: mission-point index (default `ipcruise1`).

## Returns
A `DuctedFanState{Float64}` populated from the converged design-point typed state
at `ip`.

## Example
```julia
using TASOPT
ac = setup_ducted_fan_aircraft()   # aircraft with fuel_cell_with_ducted_fan model
df = TASOPT.engine.run_ducted_fan_design_point(ac)
@show df.Fe, df.Fsp, df.pif, df.etaf
```
"""
function run_ducted_fan_design_point(ac; imission::Int=1, ip::Int=ipcruise1)
    # -----------------------------------------------------------------------
    # Initialise atmosphere ISA offset (mirrors turbofan harness pattern).
    # -----------------------------------------------------------------------
    altTO = ac.parm[imaltTO]
    T_std = atmos(altTO).T
    ac.parm[imDeltaTatm] = ac.parm[imT0TO] - T_std

    # -----------------------------------------------------------------------
    # Set ReUnit from design altitude + Mach.
    # -----------------------------------------------------------------------
    ΔTatmos = ac.parm[imDeltaTatm]
    alt_m   = ac.para[iaalt, ip, imission]
    as      = atmos(alt_m, ΔTatmos)
    Mach    = ac.para[iaMach, ip, imission]

    ac.para[iaReunit, ip, imission] = Mach * as.a * as.ρ / as.μ

    # -----------------------------------------------------------------------
    # Dual-write ambient conditions directly to typed EngineState so that
    # ductedfancalc! can read them from the typed surface (tasopt-j9l.45.15).
    # Eliminates the pare_to_engine_state! blanket sync used previously.
    # -----------------------------------------------------------------------
    let eng_ip = ac.missions[imission].points[ip].engine
        eng_ip.p0    = as.p
        eng_ip.T0    = as.T
        eng_ip.a0    = as.a
        eng_ip.rho0  = as.ρ
        eng_ip.mu0   = as.μ
        eng_ip.M0    = Mach
        eng_ip.st0.u = Mach * as.a
    end

    # -----------------------------------------------------------------------
    # Run ducted-fan design-point sizing.
    # -----------------------------------------------------------------------
    ac.engine.enginecalc!(ac, "design", imission, ip, true)

    # -----------------------------------------------------------------------
    # Read converged typed EngineState → DuctedFanState (tasopt-j9l.45.15).
    # -----------------------------------------------------------------------
    state = DuctedFanState{Float64}()
    pare_to_ducted_fan_state!(state, ac.missions[imission].points[ip].engine)
    return state
end

# ---------------------------------------------------------------------------
# run_ducted_fan_sweep
# ---------------------------------------------------------------------------

"""
    run_ducted_fan_sweep(ac; imission=1, ip_range=ipstatic:ipdescentn,
                         initializes_engine=false) -> Dict{Int,DuctedFanState{Float64}}

Run the ducted-fan off-design routine across every mission point in `ip_range`,
and return a `Dict` mapping mission-point index → `DuctedFanState`.

The aircraft must already be sized (design point converged) so that:
- Required thrust at each point is stored in the typed `DuctedFanState`.
- Fan power limit at each point is stored in the typed `DuctedFanState`.
- Ambient conditions (T0, p0, a0, M0, …) are pre-populated in typed state.

Control logic mirrors `ductedfancalc!`: takeoff/rotation points
(`ipstatic:iprotate`) operate at fixed fan power (`iPspec=true`); all other
points operate at fixed thrust target.

## Arguments
- `ac`: a ducted-fan `aircraft` object with converged design-point typed state.
- `imission::Int`: mission index (default `1`).
- `ip_range`: iterable of mission-point indices (default all 16 points).
- `initializes_engine::Bool`: passed to `enginecalc!` for each point.

## Returns
`Dict{Int, DuctedFanState{Float64}}` — one entry per `ip` in `ip_range`.

## Example
```julia
using TASOPT
ac = setup_ducted_fan_aircraft()
TASOPT.engine.run_ducted_fan_design_point(ac)   # size at cruise
states = TASOPT.engine.run_ducted_fan_sweep(ac; ip_range=ipcruise1:ipcruise2)
@show states[ipcruise1].Fe
```
"""
function run_ducted_fan_sweep(ac;
                              imission::Int=1,
                              ip_range=ipstatic:ipdescentn,
                              initializes_engine::Bool=false)
    states = Dict{Int, DuctedFanState{Float64}}()
    for ip in ip_range
        ac.engine.enginecalc!(ac, "off_design", imission, ip, initializes_engine)
        # Read converged typed EngineState → DuctedFanState (tasopt-j9l.45.15).
        state = DuctedFanState{Float64}()
        pare_to_ducted_fan_state!(state, ac.missions[imission].points[ip].engine)
        states[ip] = state
    end
    return states
end

# ---------------------------------------------------------------------------
# write_ducted_fan_sweep_toml
# ---------------------------------------------------------------------------

# Stations serialised in flow-path order.
const _DF_TOML_STATION_ORDER = (
    ("0",  "Freestream",    :st0),
    ("18", "FanFaceOuter",  :st18),
    ("2",  "FanFace",       :st2),
    ("21", "FanExit",       :st21),
    ("7",  "FanNozzle",     :st7),
    ("8",  "FanNozzleExit", :st8),
)

"""
    write_ducted_fan_sweep_toml(io::IO, states, ip_range, ac; imission=1)
    write_ducted_fan_sweep_toml(path::AbstractString, states, ip_range, ac; imission=1)

Write a ducted-fan sweep (returned by [`run_ducted_fan_sweep`](@ref)) to TOML
format suitable for use as a regression baseline fixture.

Each `[[points]]` table contains:
- Metadata: `ip`, `alt_m`, `Mach`.
- Performance: `Fe_N`, `Fsp_N_kgs`, `TSEC_W_N`, `TSFC_kg_Ns`, `Pfan_W`,
  `mfan_kg_s`, `etaf`, `pif`, `mbf`.
- Station subtables (e.g. `[points.stations.st0]`) with twelve scalar
  fields: `Tt`, `ht`, `pt`, `cpt`, `Rt`, `Ts`, `ps`, `cps`, `Rs`, `u`,
  `A`, `mdot`.

The `path` overload opens, writes, and closes the file.
"""
function write_ducted_fan_sweep_toml(io::IO, states::Dict{Int,<:DuctedFanState},
                                     ip_range, ac; imission::Int=1)
    ips = collect(Int, ip_range)
    n   = length(ips)
    points = Vector{Dict{String,Any}}(undef, n)
    for (k, ip) in enumerate(ips)
        df = states[ip]
        station_dict = Dict{String,Any}()
        for (_, _, stfld) in _DF_TOML_STATION_ORDER
            station_dict[String(stfld)] = _station_to_dict(getfield(df, stfld))
        end
        points[k] = Dict{String,Any}(
            "ip"         => ip,
            "alt_m"      => Float64(ac.para[iaalt,  ip, imission]),
            "Mach"       => Float64(ac.para[iaMach, ip, imission]),
            "M0"         => Float64(df.M0),
            "T0_K"       => Float64(df.T0),
            "p0_Pa"      => Float64(df.p0),
            "a0_m_s"     => Float64(df.a0),
            "Fe_N"       => Float64(df.Fe),
            "Fsp_N_kgs"  => Float64(df.Fsp),
            "TSEC_W_N"   => Float64(df.TSEC),
            "TSFC_kg_Ns" => Float64(df.TSFC),
            "Pfan_W"     => Float64(df.Pfan),
            "mfan_kg_s"  => Float64(df.mfan),
            "etaf"       => Float64(df.etaf),
            "pif"        => Float64(df.pif),
            "mbf"        => Float64(df.mbf),
            "stations"   => station_dict,
        )
    end
    data = Dict{String,Any}(
        "metadata" => Dict{String,Any}(
            "n_points"    => n,
            "aircraft"    => "ducted_fan_test",
            "description" => "Ducted-fan sweep baseline — test aircraft.",
        ),
        "points" => points,
    )
    TOML.print(io, data)
    return nothing
end

function write_ducted_fan_sweep_toml(path::AbstractString,
                                     states::Dict{Int,<:DuctedFanState},
                                     ip_range, ac; imission::Int=1)
    open(path, "w") do io
        write_ducted_fan_sweep_toml(io, states, ip_range, ac; imission=imission)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# regenerate_ducted_fan_baseline
# ---------------------------------------------------------------------------

"""
    const DUCTED_FAN_BASELINE_PATH

Default path for the ducted-fan sweep regression baseline fixture.
Points to `test/fixtures/ducted_fan_sweep_baseline.toml` relative to the
repository root.
"""
const DUCTED_FAN_BASELINE_PATH = joinpath(
    @__DIR__, "..", "..", "..", "test", "fixtures", "ducted_fan_sweep_baseline.toml")

"""
    regenerate_ducted_fan_baseline(ac, ip_range; path=DUCTED_FAN_BASELINE_PATH)

**Auditable action** — regenerate the ducted-fan sweep regression baseline from
a configured aircraft object and write it to `path`.

`ac` must have its ducted-fan design point already converged (design thrust
and ambient conditions stored in typed `DuctedFanState` at each `ip` in `ip_range`).

## ⚠  Warning

Regenerating the baseline changes the numerical reference used by the ducted-fan
regression test.  Only call this when a physics or model change is intentional
and the new values have been verified.  Commit the updated baseline file with an
explanatory git message.

## Usage
```julia
using TASOPT
ac = setup_ducted_fan_aircraft()   # see test/generate_ducted_fan_baseline.jl
TASOPT.engine.regenerate_ducted_fan_baseline(ac, ipcruise1:ipcruise1)
```
"""
function regenerate_ducted_fan_baseline(ac, ip_range;
                                        path::AbstractString=DUCTED_FAN_BASELINE_PATH)
    @warn """
    *** REGENERATING DUCTED-FAN SWEEP BASELINE ***
    Output path: $(abspath(path))

    This rewrites the numerical reference for ducted-fan regression tests.
    Only proceed when the change is intentional and has been verified.
    Commit the updated baseline file with an explanatory git message.
    """
    mkpath(dirname(abspath(path)))
    states = run_ducted_fan_sweep(ac; ip_range=ip_range)
    write_ducted_fan_sweep_toml(path, states, ip_range, ac)
    n = length(collect(ip_range))
    @info "Baseline written: $(abspath(path))  ($n mission points)"
    return abspath(path)
end
